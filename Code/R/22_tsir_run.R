# 22_tsir_forecast.R
# At each origin the model is refitted on the preceding TSIR_WINDOW biweeks
# alone, then the recursion is iterated forward:
#
#   I(t+i) = beta[season] * S(t+i-1) * I(t+i-1)^alpha
#   S(t+i) = births(one year earlier, same season) + S(t+i-1) - I(t+i)

# 23_tsir_process.R merges these into the file the comparison reads.
#
# Input  : the converted biweekly files
# Output : tsir_forecast/raw/w<window>/<uid>.rds, one per district



source("/.../.../Code/R/00_config.R")
source("/.../.../Code/R/05_tsir_engine.R")

library(tsiR)
library(tidyverse)

# ---- settings ----------------------------------------------------------------

TSIR_WINDOW  <- as.integer(Sys.getenv("MS_TSIR_WINDOW", "130"))
TSIR_KS      <- as.integer(strsplit(Sys.getenv("MS_TSIR_KS", "1,4,6,12"), ",")[[1]])
TSIR_SIGMAX  <- 0.5
K_MAX        <- max(TSIR_KS)

# empty runs every district; set MS_TSIR_DISTRICTS to a comma separated list of
# uids for a trial
ONLY <- Sys.getenv("MS_TSIR_DISTRICTS", "")

RAW_DIR <- file.path(OUT_DIR, "tsir_forecast", "raw",
                     paste0("w", TSIR_WINDOW))
dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)

cat("window", TSIR_WINDOW, "biweeks |",
    round(TSIR_WINDOW / BIWEEKS_PER_YEAR, 1), "years |",
    round(TSIR_WINDOW / BIWEEKS_PER_YEAR), "observations per seasonal rate\n")
cat("horizons:", paste(TSIR_KS, collapse = ", "), "\n")

# ---- data --------------------------------------------------------------------
rd <- function(f) read_csv(file.path(OUT_DIR, f), show_col_types = FALSE)

cases  <- rd("cases_biweekly.csv")
births <- rd("births_biweekly.csv") %>% select(uid, t_global, births)
popn   <- rd("popn_biweekly.csv")   %>% select(uid, t_global, population)

dat <- cases %>%
  select(state, district, uid, t_global, time, year, biweek, cases) %>%
  left_join(births, by = c("uid", "t_global")) %>%
  left_join(popn,   by = c("uid", "t_global")) %>%
  arrange(uid, t_global) %>%
  mutate(cases = replace_na(cases, 0L), births = replace_na(births, 0L))

uids <- sort(unique(dat$uid))
if (nzchar(ONLY)) {
  want <- trimws(strsplit(ONLY, ",")[[1]])
  uids <- intersect(uids, want)
  if (!length(uids)) stop("none of the requested districts are in the data")
}
cat("districts:", length(uids), "\n")

# ---- one origin ---------------------------------------------------------------
# Refits on the window ending at t, then propagates K_MAX steps.
k_step_ahead <- function(t, window, k, d) {
  fit_data <- d[(t - window + 1):t, ]
  ahead    <- d[(t + 1):(t + k), c("t_global", "time", "year", "biweek")]
  ahead$k  <- seq_len(k)
  
  # TSIR cannot fit a window that is empty or all but empty
  if (sum(fit_data$cases == 0) >= nrow(fit_data) - 1) {
    ahead$predicted <- 0
    ahead$alpha <- NA_real_
    ahead$capped <- FALSE
    ahead$status <- "window empty"
    return(ahead)
  }
  
  # the package advises the break method above 30 percent zeros
  epi <- if (mean(fit_data$cases == 0) > 0.3) "break" else "cont"
  
  fit <- tryCatch(
    runtsir_updated(
      data = as.data.frame(fit_data[, c("time", "cases", "pop", "births")]),
      IP = 2, xreg = "cumcases", regtype = "lm",
      alpha = NULL, sbar = NULL,
      method = "deterministic", nsim = 1, pred = "step-ahead",
      family = "poisson", link = "log",
      sigmamax = TSIR_SIGMAX, epidemics = epi, printon = FALSE),
    error = function(e) e)
  
  if (inherits(fit, "error")) {
    ahead$predicted <- NA_real_
    ahead$alpha <- NA_real_
    ahead$capped <- FALSE
    ahead$status <- paste("fit failed:", conditionMessage(fit))
    return(ahead)
  }
  
  n   <- length(fit$S)
  rho <- fit$rho[length(fit$rho)]
  b   <- unname(fit$beta)
  a   <- as.numeric(fit$alpha)   # this window's own, not a fixed value
  
  S <- c(fit$S, rep(NA_real_, k))
  I <- c(fit$I, rep(NA_real_, k))
  
  # the seasonal position of each forecast step, continuing from the window
  seas <- fit$period[length(fit$period)]
  
  capped <- logical(k)
  
  for (i in seq_len(k)) {
    seas <- if (seas >= BIWEEKS_PER_YEAR) 1L else seas + 1L
    
    raw <- b[seas] * S[n + i - 1] * (max(I[n + i - 1], 0) ^ a)
    
    
    I[n + i] <- min(raw, S[n + i - 1])
    capped[i] <- raw > S[n + i - 1]
    
    # Births one year earlier at the same seasonal position. The reference
    b_idx <- n + i - BIWEEKS_PER_YEAR
    b_in  <- if (b_idx >= 1 && b_idx <= nrow(fit_data))
      fit_data$births[b_idx] else mean(fit_data$births)
    
    S[n + i] <- max(b_in + S[n + i - 1] - I[n + i], 0)
  }
  
  ahead$predicted <- I[(n + 1):(n + k)] / rho
  ahead$alpha  <- a
  ahead$capped <- capped
  ahead$status <- "ok"
  ahead
}


run_district <- function(u) {
  d <- dat %>% filter(uid == u) %>% arrange(t_global)
  d$pop <- d$population
  
  first <- TSIR_WINDOW + K_MAX
  last  <- nrow(d) - K_MAX
  if (last < first) {
    cat("  ", u, ": series too short for this window\n", sep = "")
    return(NULL)
  }
  
  out <- vector("list", last - first + 1)
  for (j in seq_along(out)) {
    t <- first + j - 1
    r <- k_step_ahead(t, TSIR_WINDOW, K_MAX, d)
    r$origin_t <- d$t_global[t]
    out[[j]] <- r
  }
  
  res <- bind_rows(out) %>%
    mutate(state = d$state[1], district = d$district[1], uid = u) %>%
    left_join(d %>% select(t_global, observed = cases), by = "t_global")
  
  res
}

safe_name <- function(u) gsub("[^A-Za-z0-9]+", "_", u)

t0 <- Sys.time()
done <- 0L

for (i in seq_along(uids)) {
  u <- uids[i]
  f <- file.path(RAW_DIR, paste0(safe_name(u), ".rds"))
  
  if (file.exists(f)) {
    cat(sprintf("  %d of %d | %s | already done\n", i, length(uids), u))
    next
  }
  
  r <- run_district(u)
  if (is.null(r)) next
  
  saveRDS(r, f)
  done <- done + 1L
  
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  %d of %d | %s | %d origins | %.1f min, %.1f projected\n",
              i, length(uids), u, n_distinct(r$origin_t), el,
              el / max(done, 1) * length(uids)))
}

cat("\nfitted this run:", done, "| files in", RAW_DIR, ":",
    length(list.files(RAW_DIR, pattern = "\\.rds$")), "\n")

# ---- a look at what was written ----------------------------------------------
# The merge is 23_tsir_process.R; this is only enough to see the run behaved.
fc <- bind_rows(lapply(
  list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE), readRDS))
if (!nrow(fc)) stop("nothing was written")
fc <- fc %>% filter(k %in% TSIR_KS)

# ---- report ------------------------------------------------------------------
cat("\nforecasts:", nrow(fc), "rows |", n_distinct(fc$uid), "districts\n")

cat("\nfit status:\n")
print(as.data.frame(fc %>% filter(k == min(TSIR_KS)) %>%
                      count(status, name = "origins")))

ok <- fc %>% filter(status == "ok", is.finite(predicted), is.finite(observed))

cat("\naccuracy over every origin:\n")
print(as.data.frame(ok %>% group_by(state, k) %>%
                      summarise(n = n(),
                                rmse = round(sqrt(mean((predicted - observed)^2)), 2),
                                mae  = round(mean(abs(predicted - observed)), 2),
                                corr = round(cor(predicted, observed), 3),
                                .groups = "drop")))

# The comparison against the SFNN must use the same period, so the test years
# are reported separately as well.
test_from <- (2017 - BASE_YEAR) * BIWEEKS_PER_YEAR + 7
te <- ok %>% filter(t_global >= test_from)

if (nrow(te)) {
  cat("\nover the SFNN test period only, from t_global", test_from, ":\n")
  print(as.data.frame(te %>% group_by(state, k) %>%
                        summarise(n = n(),
                                  rmse = round(sqrt(mean((predicted - observed)^2)), 2),
                                  corr = round(cor(predicted, observed), 3),
                                  .groups = "drop")))
}


cat("\norigins where the ceiling bound:\n")
print(as.data.frame(ok %>% group_by(k) %>%
                      summarise(n = n(), capped = sum(capped),
                                pct = round(100 * mean(capped), 1), .groups = "drop")))

cat("\nalpha, estimated per window:\n")
print(summary(ok$alpha))
cat("spread within a district, median of the per district sd:",
    round(median(ok %>% group_by(uid) %>%
                   summarise(s = sd(alpha, na.rm = TRUE),
                             .groups = "drop") %>% pull(s), na.rm = TRUE), 3), "\n")

cat("\nspread of the predictions against the data, by horizon:\n")
print(as.data.frame(ok %>% group_by(k) %>%
                      summarise(sd_ratio = round(sd(predicted) / sd(observed), 3),
                                median_pred = round(median(predicted), 2),
                                median_obs  = round(median(observed), 2),
                                max_pred = round(max(predicted), 1),
                                .groups = "drop")))

cat("\ntotal time:",
    sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))), "\n")
cat("next: 23_tsir_process.R merges these into the comparison file\n")