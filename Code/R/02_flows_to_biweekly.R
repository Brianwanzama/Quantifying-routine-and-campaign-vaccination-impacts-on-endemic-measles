# 02_flows_to_biweekly.R
# Cases and births, monthly to biweekly, both states.
# Output : cases_biweekly.csv, births_biweekly.csv


source("/.../.../Code/R//00_config.R")
source("/.../...//Code/R/01_district_names.R")

library(tidyverse)
library(lubridate)

refs <- load_all_references(STATES, file_of)

stochastic_round <- function(x) {
  keep <- !is.na(x)
  out  <- rep(NA_integer_, length(x))
  fl   <- floor(x[keep])
  out[keep] <- as.integer(fl + rbinom(sum(keep), 1, x[keep] - fl))
  out
}

biweek_grid <- function(years) {
  bind_rows(lapply(years, function(yr) {
    s <- as.Date(paste0(yr, "-01-01")) + (0:(BIWEEKS_PER_YEAR - 1)) * BIWEEK_DAYS
    tibble(year = yr, biweek = seq_len(BIWEEKS_PER_YEAR),
           bw_start = s, bw_end = s + BIWEEK_DAYS - 1)
  }))
}

# ---- read one state's monthly file, long ------------------------------------
read_monthly <- function(state, what, value) {
  path <- file_of(state, what)
  
  if (is.na(path) || !file.exists(path)) {
    stop("missing ", what, " for ", state)
  }
  
  d <- read_csv(path, show_col_types = FALSE)
  
  names(d)[tolower(names(d)) == "year"]  <- "year"
  names(d)[tolower(names(d)) == "month"] <- "month"
  
  dc  <- setdiff(names(d), c("year", "month"))
  ref <- refs[[state]]
  src <- paste(state, what)
  st  <- state
  
  d %>%
    pivot_longer(
      all_of(dc),
      names_to = "district_raw",
      values_to = value
    ) %>%
    mutate(
      state = st,
      district = standardise_districts(
        district_raw,
        ref,
        source = src
      )
    )
}

# ---- allocate a monthly series onto the biweek clock ------------------------
to_biweekly <- function(long, value) {
  long <- long %>%
    mutate(month_start   = as.Date(sprintf("%04d-%02d-01", year, month)),
           month_end     = ceiling_date(month_start, "month") - days(1),
           days_in_month = as.integer(month_end - month_start) + 1L) %>%
    filter(month_start >= STUDY_START_DATE, month_start <= STUDY_END_DATE)

  years <- sort(unique(long$year))
  grid  <- biweek_grid(years)

  # the share of each month falling in each biweek, computed once
  overlap <- long %>%
    distinct(year, month, month_start, month_end, days_in_month) %>%
    crossing(grid %>% select(g_year = year, biweek, bw_start, bw_end)) %>%
    mutate(days = as.integer(pmin(month_end, bw_end) -
                             pmax(month_start, bw_start)) + 1L) %>%
    filter(days > 0) %>%
    mutate(w = days / days_in_month) %>%
    select(year, month, g_year, biweek, w)

  alloc <- long %>%
    select(state, district, year, month, all_of(value)) %>%
    left_join(overlap, by = c("year", "month")) %>%
    mutate(part = .data[[value]] * w) %>%
    group_by(state, district, year = g_year, biweek) %>%
    summarise(frac = sum(part), .groups = "drop")

  # the full grid, so a month absent from the input leaves NA rather than
  # vanishing and shortening that district's series
  full <- alloc %>% distinct(state, district) %>%
    crossing(grid) %>%
    left_join(alloc, by = c("state", "district", "year", "biweek"))

  full %>%
    mutate(!!value := stochastic_round(frac),
           uid      = make_uid(state, district),
           t_global = (year - BASE_YEAR) * BIWEEKS_PER_YEAR + biweek,
           time     = year + (biweek - 1) / BIWEEKS_PER_YEAR) %>%
    filter(!(year == min(years) & biweek < FIRST_BIWEEK_KEPT),
           bw_start <= STUDY_END_DATE) %>%
    arrange(state, district, t_global) %>%
    select(state, district, uid, year, biweek, t_global, time,
           bw_start, bw_end, all_of(value))
}

# ---- run for both states and both quantities --------------------------------
set.seed(RNG_SEED)

for (spec in list(list(what = "cases",  value = "cases",
                       out = "cases_biweekly.csv"),
                  list(what = "births", value = "births",
                       out = "births_biweekly.csv"))) {

  cat("\n=====", spec$what, "=====\n")

  long <- bind_rows(lapply(STATES, read_monthly,
                           what = spec$what, value = spec$value))
  cat("monthly rows:", nrow(long),
      "| districts:", n_distinct(paste(long$state, long$district)), "\n")

  out <- to_biweekly(long, spec$value)

  in_tot  <- sum(long[[spec$value]], na.rm = TRUE)
  out_tot <- sum(out[[spec$value]], na.rm = TRUE)
  cat("monthly total  ", format(round(in_tot), big.mark = ","), "\n")
  cat("biweekly total ", format(round(out_tot), big.mark = ","), "\n")
  cat("difference     ", sprintf("%.2f%%", 100 * abs(in_tot - out_tot) / in_tot),
      "  (the pre-April 2008 lead in, and the last two days of each December,\n",
      "   which fall outside 26 blocks of 14 days)\n", sep = "")

  cat("\nrows:", nrow(out), "| districts:", n_distinct(out$uid),
      "| t_global", min(out$t_global), "to", max(out$t_global), "\n")
  cat("per district:\n")
  print(summary(as.numeric(table(out$uid))))
  cat("unobserved biweeks:", sum(is.na(out[[spec$value]])),
      sprintf(" (%.1f%%)", 100 * mean(is.na(out[[spec$value]]))), "\n")

  by_state <- out %>% group_by(state) %>%
    summarise(districts = n_distinct(district),
              total = sum(.data[[spec$value]], na.rm = TRUE),
              na = sum(is.na(.data[[spec$value]])), .groups = "drop")
  print(as.data.frame(by_state))

  write_csv(out, file.path(OUT_DIR, spec$out))
  cat("wrote", spec$out, "\n")
}
