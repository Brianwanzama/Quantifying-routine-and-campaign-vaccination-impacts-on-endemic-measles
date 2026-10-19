# 23_tsir_process.R
# Merges the per-district TSIR forecasts into one file, and puts it on the same
# footing as the SFNN predictions so the two can be compared.
#
# 22_tsir_run.R saves every horizon per district. This selects the horizons of
# interest, restricts to the SFNN test period, and writes a single file keyed
# the same way.
#TSIR refits at every forecast origin on a window of recent
# data; the SFNN trains once on 2008 to 2016 and never updates. So TSIR sees
# fresher data at each origin while the SFNN sees more of it in total. Neither
# is the fairer arrangement, but the difference belongs in the methods.
#
# Input  : tsir_forecast/raw/w<window>/*.rds
# Output : tsir_forecast/tsir_predictions_w<window>.csv
#          tsir_forecast/tsir_accuracy_w<window>.csv


source("/.../.../Code/R/00_config.R")

library(tidyverse)

TSIR_WINDOW <- as.integer(Sys.getenv("MS_TSIR_WINDOW", "130"))
TSIR_KS     <- as.integer(strsplit(Sys.getenv("MS_TSIR_KS", "1,4,6,12"), ",")[[1]])

OUT_F   <- file.path(OUT_DIR, "tsir_forecast")
RAW_DIR <- file.path(OUT_F, "raw", paste0("w", TSIR_WINDOW))

if (!dir.exists(RAW_DIR)) stop("no raw forecasts at ", RAW_DIR)

files <- list.files(RAW_DIR, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("no district files in ", RAW_DIR)

cat("window", TSIR_WINDOW, "| district files:", length(files), "\n")
cat("horizons kept:", paste(TSIR_KS, collapse = ", "), "\n")

# ---- merge -------------------------------------------------------------------
fc <- bind_rows(lapply(files, readRDS))
cat("\nrows read:", nrow(fc), "| districts:", n_distinct(fc$uid),
    "| horizons available:", paste(sort(unique(fc$k)), collapse = ", "), "\n")

missing_k <- setdiff(TSIR_KS, unique(fc$k))
if (length(missing_k))
  cat("requested but not present:", paste(missing_k, collapse = ", "), "\n")

fc <- fc %>% filter(k %in% TSIR_KS)

# ---- how the fitting went ----------------------------------------------------
# Reported before anything is filtered, since a district that failed at most
cat("\nfit status, per origin:\n")
print(as.data.frame(fc %>% filter(k == min(TSIR_KS)) %>%
                      count(status, name = "origins") %>%
                      mutate(pct = round(100 * origins / sum(origins), 1))))

by_dist <- fc %>% filter(k == min(TSIR_KS)) %>%
  group_by(state, district, uid) %>%
  summarise(origins = n(),
            ok = sum(status == "ok"),
            pct_ok = 100 * mean(status == "ok"), .groups = "drop")

poor <- by_dist %>% filter(pct_ok < 50)
if (nrow(poor)) {
  cat("\ndistricts where the fit succeeded at under half the origins:",
      nrow(poor), "\n")
  print(as.data.frame(poor %>% arrange(pct_ok) %>% head(10) %>%
                        transmute(state, district, origins, ok, pct_ok = round(pct_ok, 1))))
  cat("these are usually the sparsest districts; consider whether to exclude\n")
  cat("them from the comparison rather than let them weight it\n")
}

# ---- the test period ---------------------------------------------------------
# The comparison must use the period the SFNN was held out on, or the two are
# not measuring the same thing.
test_from <- (2017 - BASE_YEAR) * BIWEEKS_PER_YEAR + 7

pred <- fc %>%
  filter(status == "ok", is.finite(predicted), is.finite(observed)) %>%
  mutate(period = if_else(t_global >= test_from, "test", "pre-test")) %>%
  select(any_of(c("state", "district", "uid", "k", "origin_t", "t_global",
                  "time", "year", "biweek", "observed", "predicted",
                  "alpha", "capped", "period")))

cat("\nusable forecasts:", nrow(pred), "of", nrow(fc), "\n")
print(as.data.frame(pred %>% count(period, k) %>%
                      pivot_wider(names_from = period, values_from = n)))

# ---- accuracy ----------------------------------------------------------------
acc <- pred %>%
  group_by(state, k, period) %>%
  summarise(n = n(),
            rmse = sqrt(mean((predicted - observed)^2)),
            mae  = mean(abs(predicted - observed)),
            corr = if (sd(predicted) > 0 && sd(observed) > 0)
              cor(predicted, observed) else NA_real_,
            sd_ratio = sd(predicted) / sd(observed),
            .groups = "drop")

cat("\naccuracy by state and horizon:\n")
print(as.data.frame(acc %>%
                      transmute(state, k, period, n,
                                rmse = round(rmse, 2), corr = round(corr, 3),
                                sd_ratio = round(sd_ratio, 3))))

# The spread ratio matters as much as the error. Iterated propagation can
# collapse towards a constant or run away, and RMSE alone does not say which.
flat <- acc %>% filter(period == "test", sd_ratio < 0.3)
if (nrow(flat)) {
  cat("\nwhere the forecast varies far less than the data:\n")
  print(as.data.frame(flat %>% transmute(state, k,
                                         sd_ratio = round(sd_ratio, 3))))
  cat("the error there is not informative about forecasting skill\n")
}
big <- acc %>% filter(period == "test", sd_ratio > 3)
if (nrow(big)) {
  cat("\nwhere it varies far more, which is propagation running away:\n")
  print(as.data.frame(big %>% transmute(state, k,
                                        sd_ratio = round(sd_ratio, 3))))
}

# ---- the ceiling and the alpha ------------------------------------------------
.
if ("capped" %in% names(pred)) {
  cat("\norigins where the ceiling bound:\n")
  print(as.data.frame(pred %>% group_by(state, k, period) %>%
                        summarise(n = n(), capped = sum(capped),
                                  pct = round(100 * mean(capped), 1), .groups = "drop") %>%
                        filter(capped > 0)))
}

if ("alpha" %in% names(pred)) {
  cat("\nalpha, estimated in each window:\n")
  print(as.data.frame(pred %>% filter(k == min(TSIR_KS)) %>%
                        group_by(state) %>%
                        summarise(median = round(median(alpha, na.rm = TRUE), 3),
                                  min = round(min(alpha, na.rm = TRUE), 3),
                                  max = round(max(alpha, na.rm = TRUE), 3),
                                  within_district_sd = round(median(
                                    (pred %>% filter(k == min(TSIR_KS)) %>% group_by(uid) %>%
                                       summarise(s = sd(alpha, na.rm = TRUE),
                                                 .groups = "drop"))$s, na.rm = TRUE), 3),
                                  .groups = "drop")))
  cat("a wide spread within a district means the estimate is unstable between\n")
  cat("adjacent windows, which would argue for fixing it instead\n")
}

# ---- per district, for the comparison ----------------------------------------
by_d <- pred %>%
  filter(period == "test") %>%
  group_by(state, district, uid, k) %>%
  summarise(n = n(),
            rmse = sqrt(mean((predicted - observed)^2)),
            mae  = mean(abs(predicted - observed)),
            corr = if (sd(predicted) > 0 && sd(observed) > 0)
              cor(predicted, observed) else NA_real_,
            mean_observed = mean(observed),
            .groups = "drop") %>%
  mutate(method = "TSIR")

cat("\nper district over the test period, correlation by horizon:\n")
print(as.data.frame(by_d %>% group_by(state, k) %>%
                      summarise(districts = n(),
                                median_rmse = round(median(rmse), 2),
                                median_corr = round(median(corr, na.rm = TRUE), 3),
                                negative_corr = sum(corr < 0, na.rm = TRUE),
                                .groups = "drop")))

# ---- write -------------------------------------------------------------------
f1 <- file.path(OUT_F, paste0("tsir_predictions_w", TSIR_WINDOW, ".csv"))
f2 <- file.path(OUT_F, paste0("tsir_accuracy_w", TSIR_WINDOW, ".csv"))
f3 <- file.path(OUT_F, paste0("tsir_district_errors_w", TSIR_WINDOW, ".csv"))

write_csv(pred,  f1)
write_csv(acc,   f2)
write_csv(by_d,  f3)

cat("\nnext: 24_compare_methods.R sets these against the SFNN\n")

