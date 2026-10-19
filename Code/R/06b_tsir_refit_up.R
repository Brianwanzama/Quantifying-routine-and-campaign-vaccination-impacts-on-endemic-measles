# 06_tsir_reconstruction.R
# The grid runs to SBAR_GRID_MAX rather than the package's 0.4, so a district
# whose level genuinely sits above 0.4 is distinguished from one with no
# minimum at all.
# Output : tsir_panel.csv, tsir_fit_log.csv, tsir_beta.csv

source("/.../.../00_config.R")
source("/.../.../05_tsir_engine.R")

library(tsiR)
library(tidyverse)

# ---- read --------------------------------------------------------------------
rd <- function(f) {
  p <- file.path(OUT_DIR, f)
  if (!file.exists(p)) stop("missing ", f, ", run the converters first")
  read_csv(p, show_col_types = FALSE)
}

cases <- rd("cases_biweekly.csv")
births <- rd("births_biweekly.csv") %>% select(uid, t_global, births)
popn   <- rd("popn_biweekly.csv")   %>% select(uid, t_global, population)
cov    <- rd("coverage_biweekly.csv") %>% select(uid, t_global, v1, v2)
camp   <- rd("campaign_factors.csv")  %>% select(uid, t_global, phi, p_age)

panel <- cases %>%
  select(state, district, uid, year, biweek, t_global, time, cases) %>%
  left_join(births, by = c("uid", "t_global")) %>%
  left_join(popn,   by = c("uid", "t_global")) %>%
  left_join(cov,    by = c("uid", "t_global")) %>%
  left_join(camp,   by = c("uid", "t_global")) %>%
  arrange(uid, t_global)

cat("panel:", nrow(panel), "rows |", n_distinct(panel$uid), "districts\n")
print(as.data.frame(panel %>% count(state, name = "rows")))

n_na <- sum(is.na(panel$cases))
panel <- panel %>%
  mutate(cases  = if (MISSING_CASES_ARE_ZERO) replace_na(cases, 0L) else cases,
         births = replace_na(births, 0L),
         phi    = replace_na(phi, 1),
         p_age  = replace_na(p_age, 1))
cat("unobserved case biweeks set to zero:", n_na,
    sprintf("(%.1f%%)", 100 * n_na / nrow(panel)), "\n")

stopifnot(!any(is.na(panel$population)),
          !any(is.na(panel$v1)), !any(is.na(panel$v2)))

# ---- effective births --------------------------------------------------------
# Coverage is shifted so each dose meets the birth cohort reaching that age.
# This belongs here rather than in the converter: it is a modelling choice.
lagv <- function(x, k) {
  y <- c(rep(NA_real_, k), head(x, -k))
  y[is.na(y)] <- x[1]
  y
}

panel <- panel %>%
  group_by(uid) %>% arrange(t_global, .by_group = TRUE) %>%
  mutate(v1_lag = lagv(v1, SHIFT_V1),
         v2_lag = pmin(lagv(v2, SHIFT_V2), lagv(v1, SHIFT_V1)),
         births_eff = births * (1 - E1 * v1_lag * (1 - v2_lag) -
                                    E12 * v1_lag * v2_lag)) %>%
  ungroup()

cat("effective births as a share of births:",
    sprintf("%.3f", sum(panel$births_eff) / sum(panel$births)), "\n")

# ---- has the deviance surface turned? ----------------------------------------
# Two ways it has not: the minimum sits at the grid edge, or the curve is still
# descending at the end of the range. Either means the minimum is outside what
# was searched, so the data has not determined the level.
surface_turned <- function(fit) {
  dev <- fit$loglik
  if (is.null(dev) || all(is.na(dev))) return(NA)

  n   <- length(dev)
  pos <- which.min(dev)
  if (pos >= n - SBAR_EDGE_TOL) return(FALSE)

  # still falling across the final tenth
  tail_n <- max(3L, n %/% 10L)
  tail_d <- dev[(n - tail_n + 1):n]
  if (all(diff(tail_d) < 0)) return(FALSE)

  TRUE
}

# ---- refit window: the only change from 06 -----------------------------------
# This script is 06_tsir_reconstruction.R with two changes and nothing else:
# Uttar Pradesh only, from April of CF_REFIT_YEARS_UP[1] (config), and outputs
# named with a suffix. It refits the reconstruction over years in which
# reporting was stable, for the counterfactual only; the main reconstruction
# from 06 is not read or changed, and the network does not use this one.
REFIT_FY <- CF_REFIT_YEARS_UP[1]
t0 <- (REFIT_FY - BASE_YEAR) * BIWEEKS_PER_YEAR + FIRST_BIWEEK_KEPT
panel <- panel %>% filter(state == "Uttar Pradesh", t_global >= t0)
RECON_TAG <- paste0("_up_from", REFIT_FY)
cat("refitting Uttar Pradesh from fiscal year", REFIT_FY, "(t_global", t0, ")\n")

# ---- fit ---------------------------------------------------------------------
uids <- sort(unique(panel$uid))
panel_out <- vector("list", length(uids))
log_out   <- vector("list", length(uids))
beta_out  <- vector("list", length(uids))
elapsed   <- rep(NA_real_, length(uids))

alpha_arg <- if (is.na(ALPHA_FIXED)) NULL else ALPHA_FIXED

fit_district <- function(fit_data, phi, p_age, sbar_arg, epi) {
  tryCatch(
    runtsir_campaign(
      data = fit_data, IP = 2, xreg = "cumcases", regtype = "lm",
      alpha = alpha_arg, sbar = sbar_arg, family = "gaussian",
      link = "identity", method = "negbin", nsim = 100, epidemics = epi,
      inits.fit = FALSE, phi_t = phi, p_age = p_age,
      sbar_grid_max = SBAR_GRID_MAX, printon = FALSE),
    error = function(e) e)
}

for (i in seq_along(uids)) {
  u  <- uids[i]
  df <- panel %>% filter(uid == u) %>% arrange(t_global)

  fit_data <- data.frame(time = df$time, cases = df$cases,
                         births = df$births_eff, pop = df$population)
  frac_zero <- mean(fit_data$cases == 0)
  epi <- if (frac_zero > 0.3) "break" else "cont"

  t0 <- Sys.time()

  # profile first
  fit <- fit_district(fit_data, df$phi, df$p_age, NULL, epi)
  turned <- if (inherits(fit, "error")) NA else surface_turned(fit)
  fell_back <- FALSE

  # refit with the level fixed only where the surface did not turn
  if (!inherits(fit, "error") && isFALSE(turned)) {
    refit <- fit_district(fit_data, df$phi, df$p_age, SBAR_FALLBACK, epi)
    if (!inherits(refit, "error")) { fit <- refit; fell_back <- TRUE }
  }

  elapsed[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  base_log <- tibble(state = df$state[1], district = df$district[1], uid = u,
                     frac_zero = frac_zero, total_cases = sum(fit_data$cases),
                     epidemics = epi, mean_pop = mean(fit_data$pop))

  if (inherits(fit, "error")) {
    log_out[[i]] <- base_log %>%
      mutate(fitted = FALSE, reason = conditionMessage(fit),
             sbar_identified = NA, sbar_source = "failed",
             rho = NA_real_, sbar = NA_real_, sbar_frac = NA_real_,
             alpha = NA_real_, rsquared = NA_real_, n_neg_S = NA_integer_,
             dev_min_pos = NA_integer_, dev_n_grid = NA_integer_,
             sigma_sq = NA_real_, S_start = NA_real_, I_start = NA_real_,
             campaign_removed = NA_real_)
    next
  }

  S <- fit$simS$mean
  D <- fit$simD$mean
  Z <- as.numeric(fit$Z); if (length(Z) != nrow(df)) Z <- rep(NA_real_, nrow(df))

  beta_vec <- as.numeric(fit$contact$beta)
  sigma_sq <- summary(fit$glmfit)$dispersion
  if (!is.finite(sigma_sq)) sigma_sq <- var(residuals(fit$glmfit), na.rm = TRUE)

  panel_out[[i]] <- df %>%
    mutate(susc = S, camp_removed = D, Z = Z,
           rho_t = as.numeric(1 / fit$rho),
           rho = mean(1 / fit$rho), sbar = fit$sbar,
           alpha = as.numeric(fit$alpha),
           sbar_identified = !fell_back) %>%
    select(state, district, uid, year, biweek, t_global, time,
           cases, births, births_eff, population, v1, v2, v1_lag, v2_lag,
           phi, susc, camp_removed, Z, rho_t, rho, sbar, alpha,
           sbar_identified)

  beta_out[[i]] <- tibble(uid = u, state = df$state[1],
                          district = df$district[1],
                          season = seq_along(beta_vec), beta = beta_vec)

  log_out[[i]] <- base_log %>%
    mutate(fitted = TRUE, reason = NA_character_,
           sbar_identified = !fell_back,
           sbar_source = if (fell_back) "fixed" else "profiled",
           rho = mean(1 / fit$rho), sbar = fit$sbar,
           sbar_frac = fit$sbar / mean(fit_data$pop),
           alpha = as.numeric(fit$alpha), rsquared = fit$rsquared,
           n_neg_S = sum(S <= 0),
           dev_min_pos = fit$dev_min_pos, dev_n_grid = fit$dev_n_grid,
           sigma_sq = sigma_sq,
           S_start = fit$inits[1], I_start = fit$inits[2],
           campaign_removed = sum(D, na.rm = TRUE))

  if (i %% 10 == 0)
    cat(sprintf("  %d of %d | %.1f s each\n", i, length(uids),
                mean(elapsed, na.rm = TRUE)))
}

tsir_panel <- bind_rows(panel_out)
fit_log    <- bind_rows(log_out)
tsir_beta  <- bind_rows(beta_out)

# ---- report ------------------------------------------------------------------
cat("\ndistricts fitted:", sum(fit_log$fitted), "of", nrow(fit_log), "\n")
cat("total time:", sprintf("%.1f min", sum(elapsed, na.rm = TRUE) / 60), "\n")

if (any(!fit_log$fitted)) {
  cat("\nfailures:\n")
  print(as.data.frame(fit_log %>% filter(!fitted) %>%
    select(state, district, frac_zero, total_cases, reason)))
}

ok <- fit_log %>% filter(fitted)

cat("\n---- the susceptible level ----\n")
cat("grid searched to", SBAR_GRID_MAX, "of mean population\n")
print(as.data.frame(ok %>% count(state, sbar_source)))
cat("\nidentified:", sum(ok$sbar_identified), "of", nrow(ok),
    sprintf("(%.0f%%)", 100 * mean(ok$sbar_identified)), "\n")
cat("fell back to", SBAR_FALLBACK, "of mean population:",
    sum(!ok$sbar_identified), "\n")

if (any(ok$sbar_identified)) {
  cat("\nwhere it was identified, as a share of mean population:\n")
  print(summary(ok$sbar_frac[ok$sbar_identified]))
  cat("position of the minimum on the grid:\n")
  print(summary(ok$dev_min_pos[ok$sbar_identified]))
}

# whether the two groups differ is worth knowing before any result separates
# them, since a difference in reporting or sparsity would explain the split
if (any(!ok$sbar_identified) && any(ok$sbar_identified)) {
  cat("\nhow the two groups differ:\n")
  print(as.data.frame(ok %>% group_by(sbar_identified) %>%
    summarise(n = n(),
              med_frac_zero = round(median(frac_zero), 3),
              med_cases = median(total_cases),
              med_rho = signif(median(rho), 3),
              med_r2 = round(median(rsquared, na.rm = TRUE), 3),
              med_pop = round(median(mean_pop)), .groups = "drop")))
}

cat("\n---- the fit ----\n")
cat("reporting rate:\n");   print(summary(ok$rho))
cat("mixing exponent:\n");  print(summary(ok$alpha))
cat("r squared:\n");        print(summary(ok$rsquared))
cat("break method:", sum(ok$epidemics == "break"), "districts\n")
cat("any non positive susceptible:", sum(ok$n_neg_S > 0), "\n")

if (nrow(tsir_beta)) {
  r0 <- tsir_beta %>% group_by(uid) %>%
    summarise(mean_beta = mean(beta), .groups = "drop") %>%
    left_join(ok %>% select(uid, sbar, sbar_identified), by = "uid") %>%
    mutate(R0 = mean_beta * sbar)
  cat("\nimplied reproduction number, mean beta times sbar:\n")
  print(summary(r0$R0))
  cat("measles requires 12 to 18\n")
}

# ---- write -------------------------------------------------------------------
write_csv(tsir_panel, file.path(OUT_DIR, paste0("tsir_panel", RECON_TAG, ".csv")))
write_csv(fit_log,    file.path(OUT_DIR, paste0("tsir_fit_log", RECON_TAG, ".csv")))
write_csv(tsir_beta,  file.path(OUT_DIR, paste0("tsir_beta", RECON_TAG, ".csv")))
cat("\nwrote tsir_panel.csv, tsir_fit_log.csv, tsir_beta.csv\n")
cat("seasonal rates:", nrow(tsir_beta), "rows,",
    n_distinct(tsir_beta$uid), "districts,",
    n_distinct(tsir_beta$season), "seasons each\n")
