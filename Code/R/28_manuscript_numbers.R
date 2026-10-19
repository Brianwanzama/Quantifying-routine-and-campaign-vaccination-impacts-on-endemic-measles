# 28_manuscript_numbers.R
# Prints every number the manuscript quotes, grouped by the section it
# appears in, and saves the same output to Derived/tables/manuscript_numbers.txt.
# Run after scripts 06, 15, 19, 20, 13, 18 and 25.

suppressPackageStartupMessages({ library(tidyverse); library(scales) })
source("/.../.../.../Code/R/00_config.R")

OUTF <- file.path(TAB_DIR, "manuscript_numbers.txt")
rd   <- function(p) read_csv(p, show_col_types = FALSE)
q    <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
CF_LABEL <- ""
blk  <- function(title, expr) {
  cat("\n==============================================================\n")
  cat(if (nzchar(CF_LABEL)) paste0("[", CF_LABEL, "] ") else "", title, "\n", sep = "")
  cat("==============================================================\n")
  tryCatch(expr, error = function(e) cat("  BLOCK FAILED:", conditionMessage(e), "\n"))
}
fmt <- function(md, lo, hi) sprintf("%s (%s to %s)", comma(round(md)), comma(round(lo)), comma(round(hi)))

sink(OUTF, split = TRUE)

# ---------------- Data: coverage at the ends of the period -----------------
blk("DATA | median district coverage, 2008 and 2019 (%)", {
  m <- rd(file.path(OUT_DIR, "sfnn_master.csv"))
  print(as.data.frame(m %>% filter(year %in% c(2008, 2019)) %>%
    group_by(state, year) %>%
    summarise(MCV1 = round(100 * median(v1), 1), MCV2 = round(100 * median(v2), 1),
              .groups = "drop")))
})

# ---------------- Reconstruction ----------------
blk("RESULTS: reconstruction | median reporting rate (%)", {
  m <- rd(file.path(OUT_DIR, "sfnn_master.csv"))
  print(as.data.frame(m %>% distinct(uid, .keep_all = TRUE) %>% group_by(state) %>%
    summarise(median_rho_pct = round(100 * median(rho), 2), .groups = "drop")))
})

blk("RESULTS: reconstruction | Sbar identified or fixed, and median alpha", {
  fl <- rd(file.path(OUT_DIR, "tsir_fit_log.csv"))
  cat("  columns:", paste(names(fl), collapse = ", "), "\n")
  src <- intersect(c("sbar_source", "sbar_identified", "identified"), names(fl))
  if (length(src)) print(table(fl$state, fl[[src[1]]]))
  if ("alpha" %in% names(fl)) cat("  median alpha:", round(median(fl$alpha, na.rm = TRUE), 2), "\n")
  if ("sbar_frac" %in% names(fl)) {
    idf <- fl %>% filter(sbar_source == "profiled")
    cat("  Sbar fraction, identified districts only: median", round(median(idf$sbar_frac), 3),
        "| range", round(min(idf$sbar_frac), 3), "to", round(max(idf$sbar_frac), 3), "\n")
  }
})

blk("RESULTS: reconstruction | susceptible fraction at end of series (%)", {
  sm <- rd(file.path(OUT_DIR, "sfnn", "forecast", "susceptibility_map_data.csv"))
  print(as.data.frame(sm %>% group_by(state) %>%
    summarise(median = round(median(susc_pct), 1), above_10 = sum(susc_pct >= 10),
              .groups = "drop")))
  print(as.data.frame(sm %>% filter(susc_pct >= 10) %>% arrange(desc(susc_pct)) %>%
    transmute(state, district, susc_pct = round(susc_pct, 1))))
})

# MANUSCRIPT: what vaccination averted, from the runs in Table 2

t2 <- function(f) bind_rows(lapply(names(CF_TABLE_RUN), function(st)
  rd(file.path(OUT_DIR, CF_TABLE_RUN[[st]], f)) %>% filter(state == st)))
cell <- function(md, lo, hi) {
  if (abs(md) >= 1e6) sprintf("%.2f million (%.2f to %.2f)", md / 1e6, lo / 1e6, hi / 1e6)
  else sprintf("%s (%s to %s)", comma(round(md, -3)), comma(round(lo, -3)), comma(round(hi, -3)))
}

blk("MANUSCRIPT | Table 2, infections averted (median, 95% interval)", {
  bs <- t2("bootstrap_state.csv")
  print(as.data.frame(bs %>% rowwise() %>% transmute(state,
    period = CF_TABLE_PERIOD[[state]],
    first_and_second = cell(marg_MCV1_md, marg_MCV1_lo, marg_MCV1_hi),
    second_only      = cell(marg_MCV2_md, marg_MCV2_lo, marg_MCV2_hi),
    campaign         = if (abs(marg_camp_md) < 0.5) "no campaign" else
                         cell(marg_camp_md, marg_camp_lo, marg_camp_hi),
    all_vaccination  = cell(total_averted_md, total_averted_lo, total_averted_hi))))
  print(as.data.frame(bs %>% mutate(tot = seq_MCV1_md + seq_MCV2_md + seq_camp_md) %>%
    transmute(state,
              second_as_pct_of_routine = round(100 * marg_MCV2_md / marg_MCV1_md, 1),
              share_first = round(100 * seq_MCV1_md / tot, 1),
              share_second = round(100 * seq_MCV2_md / tot, 1),
              share_campaign = round(100 * seq_camp_md / tot, 1))))
})

blk("MANUSCRIPT | calibration of the runs in Table 2 (% of reconstructed infections recovered)", {
  rows <- lapply(names(CF_TABLE_RUN), function(st) {
    r <- CF_TABLE_RUN[[st]]
    tag <- if (grepl("^scenarios_refit", r)) sub("^scenarios_refit", "", r) else ""
    fy  <- suppressWarnings(as.integer(sub(".*from", "", r)))
    m <- rd(file.path(OUT_DIR, paste0("tsir_panel", tag, ".csv"))) %>% filter(state == st)
    if (!is.na(fy)) m <- m %>% filter(t_global >= (fy - BASE_YEAR) * BIWEEKS_PER_YEAR + FIRST_BIWEEK_KEPT)
    rr <- rd(file.path(OUT_DIR, r, "scenario_replicates.csv")) %>% filter(state == st)
    m <- m %>% filter(uid %in% unique(rr$uid))
    sim <- rr %>% group_by(b) %>% summarise(o = sum(observed), .groups = "drop")
    tibble(state = st, run = r,
           reconstructed = comma(round(sum(m$cases / m$rho))),
           recovered_pct = round(100 * median(sim$o[sim$b > 0]) / sum(m$cases / m$rho)),
           at_fitted_parameters_pct = round(100 * sim$o[sim$b == 0] / sum(m$cases / m$rho)))
  })
  print(as.data.frame(bind_rows(rows)))
})

blk("MANUSCRIPT | validation of the runs in Table 2", {
  v <- t2("scenario_validation.csv")
  print(as.data.frame(v %>% group_by(state) %>%
    summarise(districts = n(), failing = sum(any_fail), .groups = "drop")))
  print(as.data.frame(v %>% filter(any_fail) %>%
    transmute(state, district, ordering = fail_order, ratio = round(ratio, 1))))
})

blk("MANUSCRIPT | infections averted per fiscal year, by component (medians), and the peak year", {
  yr <- t2("scenario_year_replicates.csv") %>% filter(b > 0) %>%
    mutate(first = none - mcv1_only, second = mcv1_only - routine, campaign = routine - observed) %>%
    group_by(state, fy, b) %>% summarise(across(c(first, second, campaign), sum), .groups = "drop") %>%
    group_by(state, fy) %>% summarise(across(c(first, second, campaign), median), .groups = "drop") %>%
    mutate(total = first + second + campaign)
  print(as.data.frame(yr %>% mutate(across(c(first, second, campaign, total), ~ comma(round(.x))))))
  print(as.data.frame(yr %>% group_by(state) %>% slice_max(total, n = 1) %>%
    transmute(state, peak_fy = fy, total = comma(round(total)))))
})

blk("MANUSCRIPT | the campaign by district, per 1,000 (top 5, Uttar Pradesh)", {
  bd <- t2("bootstrap_district.csv") %>% filter(state == "Uttar Pradesh")
  print(as.data.frame(bd %>%
    mutate(across(starts_with("marg_camp"), ~ round(1000 * .x / mean_pop, 2))) %>%
    arrange(desc(marg_camp_md)) %>% select(district, starts_with("marg_camp")) %>% head(5)))
})

blk("MANUSCRIPT | infections averted per 1,000 by district, for the map categories", {
  bd <- t2("bootstrap_district.csv")
  print(as.data.frame(bd %>% mutate(v = 1000 * total_averted_md / mean_pop) %>%
    group_by(state) %>%
    summarise(min = round(min(v), 1), p25 = round(quantile(v, .25), 1),
              median = round(median(v), 1), p75 = round(quantile(v, .75), 1),
              p90 = round(quantile(v, .9), 1), max = round(max(v), 1), .groups = "drop")))
  if (exists("CF_MAP_BREAKS")) {
    print(as.data.frame(bd %>% mutate(v = 1000 * total_averted_md / mean_pop,
                                      band = cut(v, CF_MAP_BREAKS, right = FALSE)) %>%
      count(state, band, .drop = FALSE)))
  }
})

# ---------------- Counterfactuals, one pass per scenario run ----------------
# The main run is Derived/scenarios (both states from April 2008). Runs with a
# later start for Uttar Pradesh are Derived/scenarios_up_from<year>.
cf_report <- function(SC, label, start_fy = NA, recon_tag = "") {
  CF_LABEL <<- label
  if (nzchar(recon_tag)) blk("RECONSTRUCTION of this refit | reporting rate, Sbar and alpha", {
    fl <- rd(file.path(OUT_DIR, paste0("tsir_fit_log", recon_tag, ".csv")))
    m0 <- rd(file.path(OUT_DIR, "tsir_fit_log.csv")) %>% filter(state == "Uttar Pradesh")
    cat("  median reporting rate (%): refit", round(100 * median(fl$rho, na.rm = TRUE), 2),
        "| main reconstruction", round(100 * median(m0$rho, na.rm = TRUE), 2), "\n")
    print(table(fl$sbar_source))
    idf <- fl %>% filter(sbar_source == "profiled")
    if (nrow(idf)) cat("  Sbar fraction, identified: median", round(median(idf$sbar_frac), 3),
                       "| range", round(min(idf$sbar_frac), 3), "to", round(max(idf$sbar_frac), 3), "\n")
    cat("  median alpha:", round(median(fl$alpha, na.rm = TRUE), 2), "\n")
  })
  blk("RESULTS: averted | calibration: does the observed scenario reproduce the reconstructed epidemic?", {
    m <- rd(file.path(OUT_DIR, paste0("tsir_panel", recon_tag, ".csv")))
    if (!is.na(start_fy) && "t_global" %in% names(m)) {
      t0 <- (start_fy - BASE_YEAR) * BIWEEKS_PER_YEAR + FIRST_BIWEEK_KEPT
      m <- m %>% filter(!(state == "Uttar Pradesh" & t_global < t0))
    }
    reps <- rd(file.path(SC, "scenario_replicates.csv"))
    m <- m %>% filter(uid %in% unique(reps$uid))
    recon <- m %>% group_by(state) %>% summarise(reconstructed = sum(cases / rho), .groups = "drop")
    sim <- reps %>% group_by(state, b) %>% summarise(o = sum(observed), .groups = "drop")
    print(as.data.frame(recon %>%
      left_join(sim %>% filter(b > 0) %>% group_by(state) %>%
                  summarise(boot_median = median(o), .groups = "drop"), by = "state") %>%
      left_join(sim %>% filter(b == 0) %>% select(state, replicate_zero = o), by = "state") %>%
      transmute(state,
                boot_pct = round(100 * boot_median / reconstructed),
                zero_pct = round(100 * replicate_zero / reconstructed),
                reconstructed = comma(round(reconstructed)),
                boot_median = comma(round(boot_median)),
                replicate_zero = comma(round(replicate_zero)))))
    cat("  anchored if the bootstrap median recovers most of the reconstructed total\n")
  })

  blk("RESULTS: averted | Table 2, marginal effects", {
    bs <- rd(file.path(SC, "bootstrap_state.csv"))
    print(as.data.frame(bs %>% transmute(state,
      first_and_second = fmt(marg_MCV1_md, marg_MCV1_lo, marg_MCV1_hi),
      second_only      = fmt(marg_MCV2_md, marg_MCV2_lo, marg_MCV2_hi),
      campaign         = fmt(marg_camp_md, marg_camp_lo, marg_camp_hi),
      all_vaccination  = fmt(total_averted_md, total_averted_lo, total_averted_hi))))
    print(as.data.frame(bs %>% transmute(state,
      second_as_pct_of_routine = round(100 * marg_MCV2_md / marg_MCV1_md, 1),
      all_minus_routine = comma(round(total_averted_md - marg_MCV1_md)))))
  })

  blk("RESULTS: averted | sequential shares (%)", {
    bs <- rd(file.path(SC, "bootstrap_state.csv"))
    print(as.data.frame(bs %>% mutate(tot = seq_MCV1_md + seq_MCV2_md + seq_camp_md) %>%
      transmute(state, MCV1 = round(100 * seq_MCV1_md / tot, 1),
                MCV2 = round(100 * seq_MCV2_md / tot, 1),
                campaign = round(100 * seq_camp_md / tot, 1))))
  })

  blk("RESULTS: averted | no-vaccination burden against the observed scenario (bootstrap)", {
    reps <- rd(file.path(SC, "scenario_replicates.csv"))
    print(as.data.frame(reps %>% filter(b > 0) %>% group_by(state, b) %>%
      summarise(r = sum(none) / sum(observed), .groups = "drop") %>% group_by(state) %>%
      summarise(ratio = round(median(r), 2), lo = round(q(r, .025), 2), hi = round(q(r, .975), 2))))
    r <- reps %>% filter(b > 0) %>% group_by(uid, b) %>%
      summarise(r = sum(none) / sum(observed), .groups = "drop") %>%
      group_by(uid) %>% summarise(r = median(r))
    cat("  district ratio (median over replicates): median", round(median(r$r), 1),
        "| range", round(min(r$r), 1), "to", round(max(r$r), 1), "\n")
    if ("capped" %in% names(reps))
      cat("  ceiling reached:", sum(reps$capped), "of", nrow(reps),
          sprintf("(%.1f%%)", 100 * mean(reps$capped)), "\n")
    cat("\n  for comparison only, replicate zero (no noise):\n")
    print(as.data.frame(reps %>% filter(b == 0) %>% group_by(state) %>%
      summarise(averted = comma(round(sum(none - observed))), observed = comma(round(sum(observed))),
                none = comma(round(sum(none))), .groups = "drop")))
  })

  blk("RESULTS: averted | validation checks", {
    v <- rd(file.path(SC, "scenario_validation.csv"))
    print(as.data.frame(v %>% group_by(state) %>%
      summarise(districts = n(), ordering = sum(fail_order), implausible = sum(fail_plausible),
                any = sum(any_fail), .groups = "drop")))
    print(as.data.frame(v %>% filter(any_fail) %>%
      transmute(state, district, ordering = fail_order, ratio = round(ratio, 1))))
  })

  blk("RESULTS: averted | peak fiscal year of infections averted by routine doses (bootstrap)", {
    yr <- rd(file.path(SC, "scenario_year_replicates.csv")) %>% filter(b > 0)
    print(as.data.frame(yr %>% group_by(state, fy, b) %>%
      summarise(a = sum(no_mcv1 - observed), .groups = "drop") %>%
      group_by(state, fy) %>% summarise(averted = round(median(a)), .groups = "drop") %>%
      group_by(state) %>% slice_max(averted, n = 1)))
  })

  blk("RESULTS: averted | campaign by district, per 1,000 (top 5)", {
    bd <- rd(file.path(SC, "bootstrap_district.csv"))
    print(as.data.frame(bd %>% filter(state == "Uttar Pradesh") %>%
      mutate(across(starts_with("marg_camp"), ~ round(1000 * .x / mean_pop, 2))) %>%
      arrange(desc(marg_camp_md)) %>% select(district, starts_with("marg_camp")) %>% head(5)))
  })

  blk("RESULTS: averted | remaining burden and share prevented (bootstrap); * = failed a check", {
    reps <- rd(file.path(SC, "scenario_replicates.csv")) %>% filter(b > 0)
    bd   <- rd(file.path(SC, "bootstrap_district.csv"))
    flag <- rd(file.path(SC, "scenario_validation.csv")) %>% select(uid, any_fail)
    d <- reps %>% group_by(state, district, uid, b) %>%
      summarise(none = sum(none), observed = sum(observed), .groups = "drop") %>%
      group_by(state, district, uid) %>%
      summarise(none = median(none), prevented = median(100 * (1 - observed / none)),
                .groups = "drop") %>%
      left_join(bd %>% select(uid, mean_pop), by = "uid") %>%
      left_join(flag, by = "uid") %>%
      mutate(cf_per_1000 = round(1000 * none / mean_pop, 1), prevented = round(prevented),
             district = if_else(any_fail, paste0(district, " *"), district))
    cat("  largest counterfactual burden\n")
    print(as.data.frame(d %>% group_by(state) %>% slice_max(cf_per_1000, n = 4) %>%
      select(state, district, cf_per_1000, prevented)))
    cat("  smallest share prevented\n")
    print(as.data.frame(d %>% group_by(state) %>% slice_min(prevented, n = 2) %>%
      select(state, district, prevented)))
    cat("  largest share prevented\n")
    print(as.data.frame(d %>% group_by(state) %>% slice_max(prevented, n = 3) %>%
      select(state, district, prevented)))
  })

  CF_LABEL <<- ""
}

runs <- c("scenarios",
          sort(list.files(OUT_DIR, pattern = "^scenarios_up_from[0-9]{4}$")),
          sort(list.files(OUT_DIR, pattern = "^scenarios_refit_up_from[0-9]{4}$")))
for (r in runs) {
  fy    <- suppressWarnings(as.integer(sub(".*from", "", r)))
  refit <- grepl("^scenarios_refit", r)
  tag   <- if (refit) sub("^scenarios_refit", "", r) else ""
  lab   <- if (is.na(fy)) "main: both states from 2008" else
             if (refit) paste0("Uttar Pradesh refitted and run from ", fy) else
             paste0("Uttar Pradesh run from ", fy, ", main reconstruction")
  cf_report(file.path(OUT_DIR, r), lab, fy, tag)
}

blk("SUMMARY | all counterfactual runs side by side", {
  rows <- lapply(runs, function(r) {
    fy    <- suppressWarnings(as.integer(sub(".*from", "", r)))
    refit <- grepl("^scenarios_refit", r)
    tag   <- if (refit) sub("^scenarios_refit", "", r) else ""
    SCr   <- file.path(OUT_DIR, r)
    m <- rd(file.path(OUT_DIR, paste0("tsir_panel", tag, ".csv")))
    if (!is.na(fy)) {
      t0 <- (fy - BASE_YEAR) * BIWEEKS_PER_YEAR + FIRST_BIWEEK_KEPT
      m <- m %>% filter(!(state == "Uttar Pradesh" & t_global < t0))
    }
    reps <- rd(file.path(SCr, "scenario_replicates.csv"))
    bs   <- rd(file.path(SCr, "bootstrap_state.csv"))
    recon <- m %>% filter(uid %in% unique(reps$uid)) %>%
      group_by(state) %>% summarise(recon = sum(cases / rho), .groups = "drop")
    sim <- reps %>% filter(b > 0) %>% group_by(state, b) %>%
      summarise(o = sum(observed), .groups = "drop") %>%
      group_by(state) %>% summarise(boot = median(o), .groups = "drop")
    bs %>% mutate(tot = seq_MCV1_md + seq_MCV2_md + seq_camp_md) %>%
      transmute(state, all_vacc = round(total_averted_md), campaign = round(marg_camp_md),
                seq_MCV1 = round(100 * seq_MCV1_md / tot, 1),
                seq_camp = round(100 * seq_camp_md / tot, 1)) %>%
      left_join(recon, by = "state") %>% left_join(sim, by = "state") %>%
      transmute(run = r, state, recovered_pct = round(100 * boot / recon),
                all_vacc = comma(all_vacc), campaign = comma(campaign), seq_MCV1, seq_camp)
  })
  print(as.data.frame(bind_rows(rows)))
  cat("  recovered_pct: bootstrap median of the observed scenario as % of reconstructed infections\n")
})

# ---------------- Forecasting ----------------
blk("RESULTS: forecasting | median district correlation and RMSE (Table 3)", {
  dp <- rd(file.path(TAB_DIR, "district_performance.csv"))
  print(as.data.frame(dp %>% group_by(state, k) %>%
    summarise(median_r = round(median(correlation, na.rm = TRUE), 2),
              median_rmse = round(median(rmse), 1), .groups = "drop")))
})

blk("RESULTS: forecasting | SFNN against TSIR", {
  hit <- list.files(OUT_DIR, pattern = "^compare_rmse_by_district\\.csv$",
                    recursive = TRUE, full.names = TRUE)
  if (!length(hit)) stop("compare_rmse_by_district.csv not found under Derived; rerun script 25")
  cmp <- rd(hit[1])
  print(as.data.frame(cmp %>% group_by(k) %>%
    summarise(n = n(), median_tsir = round(median(tsir_rmse), 2),
              median_sfnn = round(median(nn_rmse), 2),
              median_gain = round(median(tsir_rmse - nn_rmse), 2),
              sfnn_wins_pct = round(100 * mean(nn_rmse < tsir_rmse), 1),
              .groups = "drop")))
})

# ---------------- Covariates ----------------
blk("RESULTS: static covariates | with against without (Table 4)", {
  print(as.data.frame(rd(file.path(TAB_DIR, "nfhs_comparison.csv"))))
})

# ---------------- Attribution ----------------
blk("RESULTS: attribution | share of attribution by group (%), model with covariates", {
  imp <- rd(file.path(TAB_DIR, "shap_importance_all_horizons.csv"))
  print(as.data.frame(imp %>% filter(variant == "with_nfhs") %>%
    select(k, group, share) %>%
    pivot_wider(names_from = k, values_from = share) %>%
    mutate(across(-group, ~ round(.x, 1)))))
})

sink()
cat("\nsaved to", OUTF, "\n")
