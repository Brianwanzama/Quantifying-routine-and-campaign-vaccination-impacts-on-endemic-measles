# 15_district_performance.R
# How well the model does district by district, across every horizon.
# Input  : district_errors_k<K>.csv, sfnn_master.csv
# Output : tables/district_performance.csv, and a report


source("/.../.../Code/R/00_config.R")

library(tidyverse)

FIN <- file.path(OUT_DIR, "sfnn", "final")

KS  <- as.integer(strsplit(Sys.getenv("MS_KS", "1,4,6,12"), ",")[[1]])
VAR <- Sys.getenv("MS_VARIANT", "no_nfhs")

cat("variant", VAR, "| horizons", paste(KS, collapse = ", "), "\n")

# ---- the errors, every horizon ----------------------------------------------
read_k <- function(k) {
  f <- file.path(FIN, paste0("district_errors_k", k, ".csv"))
  if (!file.exists(f)) { cat("  no errors for k =", k, "\n"); return(NULL) }
  read_csv(f, show_col_types = FALSE) %>% filter(variant == VAR) %>%
    mutate(k = k)
}

e <- bind_rows(lapply(KS, read_k))
if (!nrow(e)) stop("no district errors found for variant ", VAR)
KS <- sort(unique(e$k))

# ---- what each district looks like, over the test period ---------------------
# Judged on the test period rather than the whole series, since that is what
# the correlations describe.
test_from <- (2017 - BASE_YEAR) * BIWEEKS_PER_YEAR + 7

m <- read_csv(file.path(OUT_DIR, "sfnn_master.csv"), show_col_types = FALSE)

ctx <- m %>%
  group_by(state, district) %>%
  summarise(pct_nonzero_all  = 100 * mean(cases > 0),
            pct_nonzero_test = 100 * mean(cases[t_global >= test_from] > 0),
            total            = sum(cases),
            mean_test        = mean(cases[t_global >= test_from]),
            pop              = mean(population),
            .groups = "drop")

d <- e %>% left_join(ctx, by = c("state", "district"))

# ---- correlation by horizon --------------------------------------------------
cat("\n---- correlation between predicted and observed, by horizon ----\n")
print(as.data.frame(d %>% group_by(state, k) %>%
  summarise(districts = n(),
            min = round(min(correlation, na.rm = TRUE), 3),
            median = round(median(correlation, na.rm = TRUE), 3),
            max = round(max(correlation, na.rm = TRUE), 3),
            .groups = "drop") %>%
  pivot_wider(names_from = state, values_from = c(min, median, max))))

# ---- the best at each horizon ------------------------------------------------
for (s in STATES) {
  cat("\n===", s, ", the five best correlations at each horizon ===\n")
  for (k in KS) {
    sub <- d %>% filter(state == s, k == !!k) %>%
      arrange(desc(correlation)) %>% head(5)
    cat("\n  k =", k, "\n")
    print(as.data.frame(sub %>%
      transmute(district, r = round(correlation, 3), rmse = round(rmse, 1),
                mean_obs = round(mean_observed, 1),
                pct_nonzero_test = round(pct_nonzero_test, 1))))
  }
}

# ---- consistent across horizons ---------------------------------------------
# A district that does well at one horizon may be a fluke. Ranking on the mean
# across horizons, and reporting the worst as well, finds those that hold up.
avg <- d %>% group_by(state, district) %>%
  summarise(mean_r = mean(correlation, na.rm = TRUE),
            min_r  = min(correlation, na.rm = TRUE),
            n_k    = n(),
            mean_rmse = mean(rmse, na.rm = TRUE),
            mean_obs  = first(mean_test),
            pct_nonzero_test = first(pct_nonzero_test),
            pop = first(pop), .groups = "drop") %>%
  arrange(state, desc(mean_r))

cat("\n---- averaged across every horizon ----\n")
for (s in STATES) {
  cat("\n===", s, ", the ten best on average ===\n")
  print(as.data.frame(avg %>% filter(state == s) %>% head(10) %>%
    transmute(district,
              mean_r = round(mean_r, 3), worst_r = round(min_r, 3),
              mean_rmse = round(mean_rmse, 1),
              mean_obs = round(mean_obs, 1),
              pct_nonzero_test = round(pct_nonzero_test, 1),
              pop = round(pop / 1000))))
}


picks <- avg %>% group_by(state) %>%
  arrange(desc(mean_r), .by_group = TRUE) %>%
  mutate(rank = row_number(), n = n()) %>%
  filter(rank %in% c(1, ceiling(n / 2), n)) %>%
  mutate(role = case_when(rank == 1 ~ "strong",
                          rank == n ~ "weak",
                          TRUE ~ "typical")) %>%
  ungroup()

print(as.data.frame(picks %>%
  transmute(state, role, district, mean_r = round(mean_r, 3),
            mean_obs = round(mean_obs, 1),
            pct_nonzero_test = round(pct_nonzero_test, 1))))

cat("\nto draw these:\n")
for (s in STATES) {
  p <- picks %>% filter(state == s, role == "strong")
  if (nrow(p)) cat("  ", s, ": ", p$district[1], "\n", sep = "")
}

# If correlation follows population or how often a district reports, that
# belongs in the results rather than only in the choice of figure.
cat("\n---- what performance tracks ----\n")
for (k in KS) {
  sub <- d %>% filter(k == !!k, is.finite(correlation))
  cat(sprintf("  k=%-3d  with population %+.3f   with how often it reports %+.3f   with mean cases %+.3f\n",
              k,
              cor(sub$correlation, sub$pop, method = "spearman", use = "complete.obs"),
              cor(sub$correlation, sub$pct_nonzero_test, method = "spearman", use = "complete.obs"),
              cor(sub$correlation, sub$mean_test, method = "spearman", use = "complete.obs")))
}

write_csv(d, file.path(TAB_DIR, "district_performance.csv"))
write_csv(avg, file.path(TAB_DIR, "district_performance_mean.csv"))
cat("\nwrote district_performance.csv and district_performance_mean.csv\n")

