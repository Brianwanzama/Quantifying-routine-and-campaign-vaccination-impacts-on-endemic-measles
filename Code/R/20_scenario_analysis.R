# 20_scenario_analysis.R

source("/.../.../.../Code/R/00_config.R")

library(tidyverse)
library(scales)

# In Uttar Pradesh. Reported cases there are near zero until the
# reporting expansion of late 2012, so a simulation started in 2008 begins from
# almost no infections. 
CF_START_FY_UP <- suppressWarnings(as.integer(Sys.getenv("MS_CF_START_FY_UP", "")))
CF_STATES      <- trimws(strsplit(Sys.getenv("MS_CF_STATES", ""), ",")[[1]])
SC_DIR <- file.path(OUT_DIR, if (is.na(CF_START_FY_UP)) "scenarios" else
                      paste0("scenarios_up_from", CF_START_FY_UP))
dir.create(SC_DIR, showWarnings = FALSE, recursive = TRUE)

N_BOOT       <- as.integer(Sys.getenv("MS_N_BOOT", CF_N_BOOT))
BLOCK_LENGTH <- CF_BLOCK_LENGTH
SEED         <- as.integer(Sys.getenv("MS_SEED", RNG_SEED))
ONLY         <- Sys.getenv("MS_BOOT_DISTRICTS", "")

set.seed(SEED)
cat("replicates:", N_BOOT, "| block:", BLOCK_LENGTH,
    "biweeks | seed:", SEED, "\n")

# MS_RECON_TAG reads a refitted reconstruction, for example "_up_from2013"
# from script 06; results then go to a folder named after it.
RECON_TAG <- Sys.getenv("MS_RECON_TAG", "")
panel <- read_csv(file.path(OUT_DIR, paste0("tsir_panel", RECON_TAG, ".csv")), show_col_types = FALSE)
beta  <- read_csv(file.path(OUT_DIR, paste0("tsir_beta", RECON_TAG, ".csv")),  show_col_types = FALSE)
flog  <- read_csv(file.path(OUT_DIR, paste0("tsir_fit_log", RECON_TAG, ".csv")), show_col_types = FALSE)
if (nzchar(RECON_TAG)) {
  SC_DIR <- file.path(OUT_DIR, paste0("scenarios_refit", RECON_TAG))
  dir.create(SC_DIR, showWarnings = FALSE, recursive = TRUE)
}

panel <- panel %>% filter(uid %in% (flog %>% filter(fitted) %>% pull(uid)))
if (length(CF_STATES) && nzchar(CF_STATES[1]))
  panel <- panel %>% filter(state %in% CF_STATES)
if (!is.na(CF_START_FY_UP)) {
  t_start <- (CF_START_FY_UP - BASE_YEAR) * BIWEEKS_PER_YEAR + FIRST_BIWEEK_KEPT
  panel <- panel %>% filter(!(state == "Uttar Pradesh" & t_global < t_start))
  cat("Uttar Pradesh starts at fiscal year", CF_START_FY_UP, "(t_global", t_start, ")\n")
}
cat("districts fitted:", n_distinct(panel$uid), "\n")

# ---- effective births, given which doses are in force --------------------------
# The Ferrari two dose correction. Setting a dose to zero removes its term.
eff_births <- function(B, v1, v2, use_v1, use_v2) {
  a <- if (use_v1) v1 else 0
  b <- if (use_v2) v2 else 0
  B * (1 - E1 * a * (1 - b) - E12 * a * b)
}

SCEN <- tribble(
  ~id,         ~label,           ~v1,   ~v2,   ~camp,
  "observed",  "observed",       TRUE,  TRUE,  TRUE,
  "no_mcv1",   "no MCV1",        FALSE, TRUE,  TRUE,
  "no_mcv2",   "no MCV2",        TRUE,  FALSE, TRUE,
  "no_camp",   "no campaign",    TRUE,  TRUE,  FALSE,
  "none",      "no vaccination", FALSE, FALSE, FALSE,
  "mcv1_only", "MCV1 only",      TRUE,  FALSE, FALSE,
  "routine",   "MCV1 and MCV2",  TRUE,  TRUE,  FALSE)

# ---- propagate one scenario ------------------------------------------------------
# Free running: incidence feeds on its own past and removals come from it.
# Capped at the susceptible pool, since infections cannot exceed the people
# they are drawn from, which the recursion does not otherwise enforce.
propagate <- function(n, b, alpha, S0, I0, be, surv, eps) {
  I <- numeric(n); S <- numeric(n)
  I[1] <- I0; S[1] <- S0
  hit <- FALSE
  for (t in 2:n) {
    raw <- b[t] * S[t - 1] * (max(I[t - 1], 0) ^ alpha) * eps[t]
    if (!is.finite(raw)) raw <- 0
    I[t] <- min(raw, S[t - 1])
    if (raw > S[t - 1]) hit <- TRUE
    S[t] <- surv[t] * max(S[t - 1] + be[t] - I[t], 0)
  }
  list(I = I, S = S, capped = hit)
}

block_resample <- function(x, block) {
  n <- length(x); nb <- ceiling(n / block)
  st <- sample.int(max(n - block + 1, 1), nb, replace = TRUE)
  unlist(lapply(st, function(s) x[s:min(s + block - 1, n)]))[seq_len(n)]
}

# ---- one district -----------------------------------------------------------------
run_district <- function(d, bd) {
  d <- d %>% arrange(t_global)
  n <- nrow(d)
  if (n < 3 * BLOCK_LENGTH) return(NULL)

  rho    <- d$rho[1]
  sbar   <- d$sbar[1]
  alpha0 <- d$alpha[1]

  I_obs  <- d$cases / rho
  # calendar biweek, used when the replicates refit the seasonal rates below
  season <- ((d$biweek - 1) %% BIWEEKS_PER_YEAR) + 1
  # The engine numbers its 26 rates by position, starting at the second biweek
  # of the series (the first transition). 
  season_engine <- ((d$biweek - FIRST_BIWEEK_KEPT - 1) %% BIWEEKS_PER_YEAR) + 1
  b0 <- bd$beta[match(season_engine, bd$season)]
  if (anyNA(b0)) b0[is.na(b0)] <- mean(bd$beta, na.rm = TRUE)

  S0 <- sbar + (if (is.finite(d$Z[1])) d$Z[1] else 0)
  I0 <- max(I_obs[1], 1)

  # the observed susceptible series, used to form the residuals that the
  # replicates resample and to offset the refit
  S_obs <- numeric(n); S_obs[1] <- S0
  for (t in 2:n)
    S_obs[t] <- d$phi[t] * max(S_obs[t - 1] +
      eff_births(d$births[t], d$v1_lag[t], d$v2_lag[t], TRUE, TRUE) -
      I_obs[t], 0)

  fitted0 <- c(I_obs[1],
               b0[2:n] * S_obs[1:(n - 1)] * pmax(I_obs[1:(n - 1)], 0) ^ alpha0)
  resid <- log(pmax(I_obs, 0.5)) - log(pmax(fitted0, 0.5))
  resid <- resid[is.finite(resid)]
  if (length(resid) < BLOCK_LENGTH) return(NULL)
  sigma <- sd(resid, na.rm = TRUE)

  inputs <- lapply(seq_len(nrow(SCEN)), function(i)
    list(be = eff_births(d$births, d$v1_lag, d$v2_lag, SCEN$v1[i], SCEN$v2[i]),
         surv = if (SCEN$camp[i]) d$phi else rep(1, n)))

  out <- vector("list", N_BOOT + 1)

  for (bnum in 0:N_BOOT) {

    if (bnum == 0) {
      # replicate zero is the point estimate: fitted parameters, no innovation
      b <- b0; alpha <- alpha0; eps <- rep(1, n)
    } else {
      # refit on a resampled series, so beta and alpha carry their uncertainty
      e <- block_resample(resid, BLOCK_LENGTH)
      I_star <- pmax(exp(log(pmax(fitted0, 0.5)) + e), 0.5)

      df <- data.frame(
        y  = log(I_star[2:n]),
        lS = log(pmax(S_obs[1:(n - 1)], 1)),
        lI = log(pmax(I_star[1:(n - 1)], 0.5)),
        s  = factor(season[2:n], levels = 1:BIWEEKS_PER_YEAR))

      fit <- tryCatch(lm(y ~ 0 + s + lI + offset(lS), data = df),
                      error = function(e) NULL)

      if (is.null(fit)) {
        b <- b0; alpha <- alpha0
      } else {
        cf <- coef(fit)
        alpha <- unname(cf[["lI"]])
        if (!is.finite(alpha) || alpha <= 0 || alpha > 1.2) alpha <- alpha0
        bs <- exp(cf[paste0("s", seq_len(BIWEEKS_PER_YEAR))])
        b <- unname(bs[season])
        b[!is.finite(b)] <- b0[!is.finite(b)]
      }

      # lognormal with mean one, so a replicate is neither inflated nor damped
      # CF_PROCESS_NOISE (config): FALSE keeps every replicate deterministic, so
      # replicates differ only in their refitted seasonal rates and alpha, and the
      # estimates, intervals and figures all come from the same deterministic model.
      # TRUE adds lognormal variation with mean one at each step.
      eps <- if (CF_PROCESS_NOISE) exp(rnorm(n, -sigma^2 / 2, sigma)) else rep(1, n)
    }

    sims <- lapply(inputs, function(z)
      propagate(n, b, alpha, S0, I0, z$be, z$surv, eps))

    I_mat <- vapply(sims, function(z) z$I, numeric(n))
    colnames(I_mat) <- SCEN$id

    row <- bind_cols(
      tibble(b = bnum, t_global = d$t_global, year = d$year,
             biweek = d$biweek,
             capped = any(vapply(sims, function(z) z$capped, logical(1)))),
      as_tibble(I_mat))

    # the susceptible series is kept for the point estimate only, so the
    # bounds check has something to test without carrying it 1000 times
    if (bnum == 0) {
      S_mat <- vapply(sims, function(z) z$S, numeric(n))
      row$S_min <- apply(S_mat, 1, min)
      row$S_max <- apply(S_mat, 1, max)
      row$population <- d$population
    } else {
      row$S_min <- NA_real_; row$S_max <- NA_real_
      row$population <- NA_real_
    }
    out[[bnum + 1]] <- row
  }

  bind_rows(out) %>%
    mutate(state = d$state[1], district = d$district[1], uid = d$uid[1])
}

# ---- run ---------------------------------------------------------------------------
uids <- sort(unique(panel$uid))
if (nzchar(ONLY)) {
  uids <- intersect(uids, trimws(strsplit(ONLY, ",")[[1]]))
  if (!length(uids)) stop("none of the requested districts are present")
}
cat("running", length(uids), "districts\n")

t0 <- Sys.time()
res <- vector("list", length(uids))

res_ts <- vector("list", length(uids))
# year totals per replicate, so a yearly figure can carry an interval rather
# than only a point
res_fy <- vector("list", length(uids))

res_md <- vector("list", length(uids))

for (i in seq_along(uids)) {
  bd <- beta %>% filter(uid == uids[i])
  if (!nrow(bd)) next
  r <- run_district(panel %>% filter(uid == uids[i]), bd)
  if (is.null(r)) { cat("  ", uids[i], ": too short\n", sep = ""); next }

  res[[i]] <- r %>%
    group_by(state, district, uid, b) %>%
    summarise(across(all_of(SCEN$id), sum), capped = any(capped),
              .groups = "drop")

  # every scenario, including the sequential ones, since the figures need
  # MCV1 only as well as the marginal set
  res_ts[[i]] <- r %>% filter(b == 0) %>%
    select(state, district, uid, t_global, year, biweek,
           all_of(SCEN$id), S_min, S_max, population)

  res_fy[[i]] <- r %>%
    mutate(fy = if_else(biweek < FIRST_BIWEEK_KEPT, year - 1L, year)) %>%
    group_by(state, district, uid, fy, b) %>%
    summarise(across(all_of(SCEN$id), sum), .groups = "drop")

  res_md[[i]] <- r %>% filter(b > 0) %>%
    mutate(seq_MCV1 = none - mcv1_only, seq_MCV2 = mcv1_only - routine,
           seq_camp = routine - observed) %>%
    group_by(state, district, uid, t_global, year, biweek) %>%
    summarise(across(c(all_of(SCEN$id), seq_MCV1, seq_MCV2, seq_camp), median),
              .groups = "drop")

  if (i %% 10 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  %d of %d | %.1f min, %.1f projected\n",
                i, length(uids), el, el / i * length(uids)))
  }
}

reps <- bind_rows(res)
if (!nrow(reps)) stop("nothing simulated")
cat("\nreplicate rows:", nrow(reps), "|", n_distinct(reps$uid), "districts\n")

# ---- validation ---------------------------------------------------------------------------

# 1. ordering. Incidence should fall as protection is added.
ord <- ts %>%
  group_by(state, district, uid) %>%
  summarise(bad_order = sum(none < mcv1_only - 1e-6 |
                            mcv1_only < routine - 1e-6 |
                            routine < observed - 1e-6),
            n = n(), .groups = "drop") %>%
  mutate(fail_order = bad_order > 0)

# 2. bounds. No susceptible pool below zero or above population.
bnd <- ts %>%
  group_by(state, district, uid) %>%
  summarise(bad_bounds = sum(S_min < -1e-6 | S_max > population * 1.001,
                             na.rm = TRUE), .groups = "drop") %>%
  mutate(fail_bounds = bad_bounds > 0)

# 3. sign. Removing a component should never reduce cases.
sgn <- reps %>% filter(b == 0) %>%
  group_by(state, district, uid) %>%
  summarise(bad_sign = sum(no_mcv1 < observed - 1e-6) +
                       sum(no_mcv2 < observed - 1e-6) +
                       sum(no_camp < observed - 1e-6), .groups = "drop") %>%
  mutate(fail_sign = bad_sign > 0)

# 4. plausibility. The counterfactual burden should be between one and twenty
# times observed; outside that the simulation has not behaved.
pls <- reps %>% filter(b == 0) %>%
  group_by(state, district, uid) %>%
  summarise(ratio = sum(none) / pmax(sum(observed), 1), .groups = "drop") %>%
  mutate(fail_plausible = !is.finite(ratio) | ratio < 1 | ratio > 20)

valid <- ord %>%
  left_join(bnd %>% select(uid, bad_bounds, fail_bounds), by = "uid") %>%
  left_join(sgn %>% select(uid, bad_sign, fail_sign), by = "uid") %>%
  left_join(pls %>% select(uid, ratio, fail_plausible), by = "uid") %>%
  mutate(any_fail = fail_order | fail_bounds | fail_sign | fail_plausible)

cat("\n---- validation ----\n")
print(as.data.frame(valid %>% group_by(state) %>%
  summarise(districts = n(),
            ordering = sum(fail_order),
            bounds = sum(fail_bounds),
            sign = sum(fail_sign),
            implausible = sum(fail_plausible),
            any = sum(any_fail), .groups = "drop")))

if (any(valid$any_fail)) {
  cat("\ndistricts failing a check:\n")
  print(as.data.frame(valid %>% filter(any_fail) %>%
    transmute(state, district, bad_order, bad_bounds, bad_sign,
              ratio = round(ratio, 1)) %>%
    arrange(desc(bad_order + bad_bounds + bad_sign)) %>% head(15)))
  cat("\nthese are flagged, not excluded. Decide whether to drop them and\n")
  cat("say so, rather than letting the count shrink silently.\n")
} else {
  cat("\nevery district passes every check\n")
}

cat("\ncounterfactual against observed, the ratio:\n")
print(summary(pls$ratio))

write_csv(valid, file.path(SC_DIR, "scenario_validation.csv"))

# ---- did the propagation behave? -------------------------------------------------------
cat("\nreplicates where the susceptible ceiling bound:", sum(reps$capped),
    "of", nrow(reps), sprintf("(%.1f%%)", 100 * mean(reps$capped)), "\n")
if (mean(reps$capped) > 0.1)
  cat("  a large share means the free running recursion is diverging and the\n  cap is doing the work; the counterfactuals should not be trusted\n")

cat("\ndoes the observed scenario reproduce the reconstructed epidemic?\n")
recon <- panel %>% group_by(state) %>% summarise(reconstructed = sum(cases / rho), .groups = "drop")
sim <- reps %>% group_by(state, b) %>% summarise(o = sum(observed), .groups = "drop")
print(as.data.frame(recon %>%
  left_join(sim %>% filter(b > 0) %>% group_by(state) %>%
              summarise(boot_md = median(o), boot_lo = quantile(o, .025),
                        boot_hi = quantile(o, .975), .groups = "drop"), by = "state") %>%
  left_join(sim %>% filter(b == 0) %>% select(state, replicate_zero = o), by = "state") %>%
  mutate(boot_pct = round(100 * boot_md / reconstructed),
         zero_pct = round(100 * replicate_zero / reconstructed),
         across(c(reconstructed, boot_md, boot_lo, boot_hi, replicate_zero), round))))
cat("a simulation that is anchored should recover most of the reconstructed total\n")

cat("\ncounterfactual burden against observed, at the point estimate:\n")
print(as.data.frame(reps %>% filter(b == 0) %>% group_by(state) %>%
  summarise(observed = round(sum(observed)), none = round(sum(none)),
            ratio = round(sum(none) / sum(observed), 2), .groups = "drop")))
cat("a ratio far above ten or below one would mean the no vaccination world\n")
cat("is not plausible and the simulation needs revisiting\n")

# ---- the two decompositions ---------------------------------------------------------------
d_rep <- reps %>%
  mutate(marg_MCV1 = no_mcv1 - observed,
         marg_MCV2 = no_mcv2 - observed,
         marg_camp = no_camp - observed,
         total_averted = none - observed,
         seq_MCV1 = none - mcv1_only,
         seq_MCV2 = mcv1_only - routine,
         seq_camp = routine - observed)

VARS <- c("marg_MCV1", "marg_MCV2", "marg_camp", "total_averted",
          "seq_MCV1", "seq_MCV2", "seq_camp")

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))

boot_d <- d_rep %>% filter(b > 0) %>%
  group_by(state, district, uid) %>%
  summarise(across(all_of(VARS),
                   list(md = ~ median(.x), lo = ~ q(.x, .025),
                        hi = ~ q(.x, .975))), .groups = "drop") %>%
  left_join(panel %>% group_by(uid) %>%
              summarise(mean_pop = mean(population), .groups = "drop"),
            by = "uid")

boot_s <- d_rep %>% filter(b > 0) %>%
  group_by(state, b) %>%
  summarise(across(all_of(VARS), sum), .groups = "drop") %>%
  group_by(state) %>%
  summarise(across(all_of(VARS),
                   list(md = ~ median(.x), lo = ~ q(.x, .025),
                        hi = ~ q(.x, .975))), .groups = "drop")

point <- d_rep %>% filter(b == 0) %>% group_by(state) %>%
  summarise(across(all_of(VARS), sum), .groups = "drop")

# ---- report -----------------------------------------------------------------------------------
fmt <- function(md, lo, hi)
  sprintf("%s (%s to %s)", comma(round(md)), comma(round(lo)), comma(round(hi)))

cat("\n---- marginal: what removing one component would have cost ----\n")
print(as.data.frame(boot_s %>%
  transmute(state,
            MCV1 = fmt(marg_MCV1_md, marg_MCV1_lo, marg_MCV1_hi),
            MCV2 = fmt(marg_MCV2_md, marg_MCV2_lo, marg_MCV2_hi),
            Campaign = fmt(marg_camp_md, marg_camp_lo, marg_camp_hi),
            all_three = fmt(total_averted_md, total_averted_lo,
                            total_averted_hi))))

cat("\nthese overlap and do not sum: removing one dose leaves the other still\n")
cat("protecting the same children\n")

cat("\noverlap, the three marginal effects against removing all:\n")
print(as.data.frame(boot_s %>%
  transmute(state,
            sum_marginal = round(marg_MCV1_md + marg_MCV2_md + marg_camp_md),
            all_three = round(total_averted_md),
            ratio = round((marg_MCV1_md + marg_MCV2_md + marg_camp_md) /
                            total_averted_md, 2))))

cat("\n---- sequential: the attribution, which sums exactly ----\n")
print(as.data.frame(boot_s %>%
  mutate(tot = seq_MCV1_md + seq_MCV2_md + seq_camp_md) %>%
  transmute(state,
            MCV1 = sprintf("%.1f%%", 100 * seq_MCV1_md / tot),
            MCV2 = sprintf("%.1f%%", 100 * seq_MCV2_md / tot),
            Campaign = sprintf("%.1f%%", 100 * seq_camp_md / tot),
            total = comma(round(tot)))))
cat("\nthese depend on the order the components are introduced\n")

cat("\ndistricts where a marginal effect's interval reaches zero:\n")
print(as.data.frame(boot_d %>% group_by(state) %>%
  summarise(districts = n(),
            MCV1 = sum(marg_MCV1_lo <= 0),
            MCV2 = sum(marg_MCV2_lo <= 0),
            Campaign = sum(marg_camp_lo <= 0), .groups = "drop")))
cat("those components are not established there and should be reported as\n")
cat("such rather than as a small positive number\n")

write_csv(reps,   file.path(SC_DIR, "scenario_replicates.csv"))
write_csv(ts, file.path(SC_DIR, "scenario_series_point.csv"))
write_csv(bind_rows(res_md), file.path(SC_DIR, "scenario_series_median.csv"))
write_csv(bind_rows(res_fy),
          file.path(SC_DIR, "scenario_year_replicates.csv"))
write_csv(boot_d, file.path(SC_DIR, "bootstrap_district.csv"))
write_csv(boot_s, file.path(SC_DIR, "bootstrap_state.csv"))
write_csv(point,  file.path(SC_DIR, "scenario_point.csv"))

