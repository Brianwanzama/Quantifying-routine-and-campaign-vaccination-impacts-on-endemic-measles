# 07_build_master.R
# The campaign is carried as mu_mr, the share of the susceptible pool removed,
# which is one minus the survival factor. West Bengal has no campaign, so its
# districts carry zero throughout.
#
# Output : sfnn_master.csv


source("/.../.../Code/R/00_config.R")
source("/.../.../Code/R/01_district_names.R")

library(tidyverse)

FILE_NFHS <- Sys.getenv("MS_NFHS_FILE",
                        file.path(OUT_DIR, "all_states_district_covariates.csv"))

refs <- load_all_references(STATES, file_of)

# ---- the reconstruction ------------------------------------------------------
panel <- read_csv(file.path(OUT_DIR, "tsir_panel.csv"), show_col_types = FALSE)
cat("panel:", nrow(panel), "rows |", n_distinct(panel$uid), "districts\n")
cat("t_global", min(panel$t_global), "to", max(panel$t_global), "\n")

# ---- coordinates -------------------------------------------------------------
# Two rows per file, longitude then latitude, one column per district.
read_coords <- function(state) {
  co  <- read_csv(file_of(state, "coords"), show_col_types = FALSE)
  lab <- tolower(trimws(co[[1]]))
  dc  <- names(co)[-1]
  ref <- refs[[state]]
  st  <- state

  tibble(district_raw = dc,
         lon = as.numeric(co[grep("^lon", lab)[1], dc]),
         lat = as.numeric(co[grep("^lat", lab)[1], dc])) %>%
    mutate(state = st,
           district = standardise_districts(district_raw, ref,
                                            source = paste(st, "coords"))) %>%
    filter(!is.na(district), !is.na(lat), !is.na(lon)) %>%
    mutate(uid = make_uid(state, district)) %>%
    select(uid, lat, lon)
}

coords <- bind_rows(lapply(STATES, read_coords))
cat("coordinates:", nrow(coords), "districts\n")

miss_co <- setdiff(unique(panel$uid), coords$uid)
if (length(miss_co)) {
  cat("no coordinates for", length(miss_co), "districts:\n")
  for (u in head(miss_co, 10)) cat("  ", u, "\n")
}

# ---- static covariates -------------------------------------------------------
# n_children_vax records how many children underlie each district's vaccination
# covariates. It is metadata rather than a covariate.
if (!file.exists(FILE_NFHS))
  stop("no covariate file at ", FILE_NFHS,
       "\nset MS_NFHS_FILE to its location")

nfhs_raw <- read_csv(FILE_NFHS, show_col_types = FALSE)

# the file names states in lower case, so they are matched case insensitively
state_of <- setNames(STATES, tolower(STATES))
nfhs <- nfhs_raw %>%
  mutate(state_std = unname(state_of[tolower(trimws(state))]))

unknown <- unique(nfhs$state[is.na(nfhs$state_std)])
if (length(unknown))
  cat("covariate rows with an unrecognised state:",
      paste(unknown, collapse = ", "), "\n")

nfhs <- nfhs %>%
  filter(!is.na(state_std)) %>%
  group_by(state_std) %>%
  group_modify(~ {
    ref <- refs[[.y$state_std]]
    .x %>% mutate(district_std = standardise_districts(
      district, ref, source = paste(.y$state_std, "NFHS")))
  }) %>%
  ungroup() %>%
  filter(!is.na(district_std)) %>%
  mutate(uid = make_uid(state_std, district_std)) %>%
  select(-any_of(c("state", "district", "state_std", "district_std",
                   "n_children_vax")))

cov_cols <- setdiff(names(nfhs), "uid")
cat("static covariates:", length(cov_cols), "|",
    paste(cov_cols, collapse = ", "), "\n")

miss_cov <- setdiff(unique(panel$uid), nfhs$uid)
if (length(miss_cov)) {
  cat("no covariates for", length(miss_cov), "districts:\n")
  for (u in head(miss_cov, 10)) cat("  ", u, "\n")
}

# ---- join --------------------------------------------------------------------
master <- panel %>%
  left_join(coords, by = "uid") %>%
  left_join(nfhs,   by = "uid") %>%
  arrange(uid, t_global)

# ---- campaign intensity ------------------------------------------------------
# The survival factor is the share of the pool remaining, so the share removed
# is one minus it. Zero outside a window, and zero throughout for West Bengal.
master <- master %>%
  mutate(phi       = replace_na(phi, 1),
         mu_mr     = 1 - phi,
         mr_active = as.integer(mu_mr > 0))

# ---- checks ------------------------------------------------------------------
cat("\nrows:", nrow(master), "| districts:", n_distinct(master$uid), "\n")
print(as.data.frame(master %>%
  group_by(state) %>%
  summarise(districts = n_distinct(uid), rows = n(),
            cases = sum(cases, na.rm = TRUE),
            campaign_biweeks = sum(mr_active), .groups = "drop")))

# The loader shifts within a district by position, so every district must share
# the same complete grid or a lag would pull the wrong time point.
g <- master %>% group_by(uid) %>%
  summarise(n = n(), first = min(t_global), last = max(t_global), .groups = "drop")
if (n_distinct(g$n) != 1 || n_distinct(g$first) != 1 || n_distinct(g$last) != 1) {
  cat("\ndistricts do not share one grid, the loader would misalign lags:\n")
  print(as.data.frame(g %>% filter(n != max(n) | first != min(first) |
                                     last != max(last)) %>% head(20)))
  stop("grid mismatch")
}
cat("every district shares the same", g$n[1], "step grid\n")

na_report <- master %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_na") %>%
  filter(n_na > 0)
if (nrow(na_report)) {
  cat("\nmissing values:\n"); print(as.data.frame(na_report))
} else {
  cat("no missing values\n")
}

cat("\nbiweeks with no cases:",
    sprintf("%.0f%%", 100 * mean(master$cases == 0, na.rm = TRUE)), "\n")
cat("susceptible fraction:\n")
print(summary(master$susc / master$population))

cat("\ncases by year:\n")
print(as.data.frame(master %>% group_by(state, year) %>%
  summarise(cases = sum(cases, na.rm = TRUE),
            zero = round(mean(cases == 0), 2), .groups = "drop") %>%
  pivot_wider(names_from = state, values_from = c(cases, zero))))

# ---- write -------------------------------------------------------------------
out <- master %>%
  select(state, district, uid, t_global, time, year, biweek,
         cases, births, births_eff, population, v1, v2, v1_lag, v2_lag,
         susc, camp_removed, Z, phi, mu_mr, mr_active,
         rho_t, rho, sbar, alpha, sbar_identified,
         lat, lon, all_of(cov_cols))

write_csv(out, file.path(OUT_DIR, "sfnn_master.csv"))
cat("\nwrote sfnn_master.csv:", nrow(out), "rows,", ncol(out), "columns\n")
