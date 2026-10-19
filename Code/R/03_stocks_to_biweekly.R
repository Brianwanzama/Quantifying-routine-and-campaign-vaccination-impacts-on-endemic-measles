
#03_stocks_to_biweekly.R
#Population and coverage, annual to biweekly, both states.
#Output : popn_biweekly.csv, coverage_biweekly.csv

source("/.../.../Code/R/00_config.R")
source("/.../.../Code/R/01_district_names.R")

library(tidyverse)

refs <- load_all_references(STATES, file_of)

clean_number <- function(x) {
  if (is.numeric(x)) x else
    suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "", as.character(x))))
}

biweek_grid <- function(years) {
  bind_rows(lapply(years, function(yr) {
    s <- as.Date(paste0(yr, "-01-01")) + (0:(BIWEEKS_PER_YEAR - 1)) * BIWEEK_DAYS
    tibble(year = yr, biweek = seq_len(BIWEEKS_PER_YEAR),
           bw_start = s, bw_end = s + BIWEEK_DAYS - 1,
           bw_mid = s + (BIWEEK_DAYS %/% 2))
  }))
}

read_annual <- function(state, what, value) {
  path <- file_of(state, what)
  
  if (is.na(path) || !file.exists(path)) {
    stop("missing ", what, " for ", state)
  }
  
  d <- read_csv(path, show_col_types = FALSE)
  names(d)[tolower(names(d)) == "year"] <- "year"
  
  dc  <- setdiff(names(d), "year")
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
      ),
      !!value := clean_number(.data[[value]])
    )
}
# ---- interpolate an annual series onto biweek midpoints ---------------------
interp_vec <- function(yrs, vals, xout) {
  ok <- !is.na(vals)
  if (!any(ok)) return(rep(NA_real_, length(xout)))
  if (sum(ok) == 1) return(rep(vals[ok][1], length(xout)))
  anchor <- as.Date(sprintf("%04d-%02d-%02d", yrs[ok],
                            POP_ANCHOR_MONTH, POP_ANCHOR_DAY))
  o <- order(anchor)
  approx(as.numeric(anchor)[o], vals[ok][o], xout = as.numeric(xout), rule = 2)$y
}

to_biweekly <- function(long, values) {
  years <- min(long$year):max(long$year)
  grid  <- biweek_grid(years)

  long %>%
    group_by(state, district) %>%
    group_modify(~ {
      out <- grid
      for (v in values) out[[v]] <- interp_vec(.x$year, .x[[v]], grid$bw_mid)
      out
    }) %>%
    ungroup() %>%
    mutate(uid      = make_uid(state, district),
           t_global = (year - BASE_YEAR) * BIWEEKS_PER_YEAR + biweek,
           time     = year + (biweek - 1) / BIWEEKS_PER_YEAR) %>%
    filter(!(year == min(years) & biweek < FIRST_BIWEEK_KEPT),
           bw_start <= STUDY_END_DATE) %>%
    arrange(state, district, t_global)
}

# ---- population --------------------------------------------------------------
cat("===== population =====\n")

popn <- bind_rows(lapply(STATES, read_annual, what = "popn", value = "population"))
cat("annual rows:", nrow(popn),
    "| districts:", n_distinct(paste(popn$state, popn$district)),
    "| years", min(popn$year), "to", max(popn$year), "\n")

pop_bw <- to_biweekly(popn, "population") %>%
  mutate(population = round(population)) %>%
  select(state, district, uid, year, biweek, t_global, time,
         bw_start, bw_end, population)

cat("rows:", nrow(pop_bw), "| districts:", n_distinct(pop_bw$uid),
    "| t_global", min(pop_bw$t_global), "to", max(pop_bw$t_global), "\n")
cat("missing:", sum(is.na(pop_bw$population)), "\n")
cat("range:", format(min(pop_bw$population, na.rm = TRUE), big.mark = ","),
    "to", format(max(pop_bw$population, na.rm = TRUE), big.mark = ","), "\n")
n_zero <- sum(pop_bw$population == 0, na.rm = TRUE)
if (n_zero)
  cat("zero population in", n_zero, "biweeks across",
      n_distinct(pop_bw$uid[pop_bw$population == 0]),
      "districts, which are those created part way through\n")

write_csv(pop_bw, file.path(OUT_DIR, "popn_biweekly.csv"))
cat("wrote popn_biweekly.csv\n")

# ---- coverage ----------------------------------------------------------------
cat("\n===== coverage =====\n")

mcv1 <- bind_rows(lapply(STATES, read_annual, what = "mcv1", value = "v1"))
mcv2 <- bind_rows(lapply(STATES, read_annual, what = "mcv2", value = "v2"))
cov <- full_join(mcv1, mcv2, by = c("state", "district", "year"))

cat("annual rows:", nrow(cov),
    "| districts:", n_distinct(paste(cov$state, cov$district)),
    "| years", min(cov$year), "to", max(cov$year), "\n")

if (COVERAGE_IS_PERCENT) cov <- cov %>% mutate(v1 = v1 / 100, v2 = v2 / 100)

n_gt <- sum(cov$v2 > cov$v1, na.rm = TRUE)
if (n_gt) {
  cat("MCV2 exceeds MCV1 in", n_gt, "district years, capped at MCV1\n")
  cov <- cov %>% mutate(v2 = pmin(v2, v1))
}

cat("\nmean coverage by year, for the MCV2 introduction:\n")
print(as.data.frame(cov %>% group_by(state, year) %>%
  summarise(v1 = round(100 * mean(v1, na.rm = TRUE), 1),
            v2 = round(100 * mean(v2, na.rm = TRUE), 1),
            v2_zero = sum(v2 == 0, na.rm = TRUE), .groups = "drop")))

cov_bw <- to_biweekly(cov, c("v1", "v2")) %>%
  mutate(v1 = pmin(pmax(v1, 0), 1),
         v2 = pmin(pmin(pmax(v2, 0), 1), v1)) %>%
  select(state, district, uid, year, biweek, t_global, time,
         bw_start, bw_end, v1, v2)

cat("\nrows:", nrow(cov_bw), "| districts:", n_distinct(cov_bw$uid), "\n")
cat("missing  v1:", sum(is.na(cov_bw$v1)), " v2:", sum(is.na(cov_bw$v2)), "\n")
cat("v1 range:", sprintf("%.3f to %.3f",
    min(cov_bw$v1, na.rm = TRUE), max(cov_bw$v1, na.rm = TRUE)), "\n")
cat("v2 range:", sprintf("%.3f to %.3f",
    min(cov_bw$v2, na.rm = TRUE), max(cov_bw$v2, na.rm = TRUE)), "\n")
cat("v2 above v1 anywhere:", sum(cov_bw$v2 > cov_bw$v1, na.rm = TRUE), "\n")

write_csv(cov_bw, file.path(OUT_DIR, "coverage_biweekly.csv"))
cat("wrote coverage_biweekly.csv\n")

# ---- do the four files agree? ------------------------------------------------
# A mismatch here would surface much later as a silent join failure.
cat("\n===== alignment =====\n")
ca <- read_csv(file.path(OUT_DIR, "cases_biweekly.csv"), show_col_types = FALSE)
bi <- read_csv(file.path(OUT_DIR, "births_biweekly.csv"), show_col_types = FALSE)

for (nm in c("cases", "births", "population", "coverage")) {
  d <- switch(nm, cases = ca, births = bi, population = pop_bw, coverage = cov_bw)
  cat(sprintf("%-11s rows %7d  districts %4d  t_global %4d to %4d\n",
              nm, nrow(d), n_distinct(d$uid), min(d$t_global), max(d$t_global)))
}

sets <- list(ca$uid, bi$uid, pop_bw$uid, cov_bw$uid)
cat("districts in all four:", length(Reduce(intersect, sets)), "\n")
for (nm in c("births", "population", "coverage")) {
  d <- switch(nm, births = bi, population = pop_bw, coverage = cov_bw)
  miss <- setdiff(unique(ca$uid), unique(d$uid))
  if (length(miss)) {
    cat("in cases but not", nm, ":", length(miss), "\n")
    for (u in head(miss, 10)) cat("   ", u, "\n")
  }
}

