# 04_campaign_factors.R
# The campaign removes mu = e * kappa of the susceptible children it reaches.
# That total is spread across the biweeks of the delivery window by the share
# of campaign days in each, then delayed one biweek for the time to protection.
#
#     Phi_t = exp(-lambda * w_{t-l}),   lambda = -log(1 - mu)

# Output : campaign_factors.csv, campaign_log.csv

source("/.../.../00_config.R")
source("/.../.../Code/R/01_district_names.R")

library(tidyverse)
library(lubridate)

refs <- load_all_references(STATES, file_of)

clean_number <- function(x) {
  if (is.numeric(x)) x else
    suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "", as.character(x))))
}

parse_date_dmy <- function(x) {
  x <- trimws(as.character(x))
  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
  for (f in c("%d-%m-%Y", "%d/%m/%Y", "%Y-%m-%d", "%Y/%m/%d")) {
    bad <- is.na(out)
    if (!any(bad)) break
    out[bad] <- suppressWarnings(as.Date(x[bad], format = f))
  }
  out
}

# ---- the grid comes from the cases file, so the output aligns ----------------
cases <- read_csv(file.path(OUT_DIR, "cases_biweekly.csv"), show_col_types = FALSE)
grid  <- cases %>% distinct(t_global, year, biweek, bw_start, bw_end) %>%
  arrange(t_global)
allu  <- cases %>% distinct(state, district, uid)

cat("grid:", nrow(grid), "biweeks | districts:", nrow(allu), "\n")

# ---- read whichever states have a campaign ----------------------------------
read_campaign <- function(state) {
  path <- file_of(state, "campaign")
  if (is.na(path) || !file.exists(path)) {
    cat(state, ": no campaign file, Phi held at 1\n", sep = "")
    return(NULL)
  }
  d <- read_csv(path, show_col_types = FALSE)

  pick <- function(...) {
    for (nm in c(...)) if (nm %in% names(d)) return(nm)
    stop("none of ", paste(c(...), collapse = ", "), " in ", basename(path))
  }
  dc <- pick("District", "district")
  sc <- pick("Start Date", "start_date", "StartDate")
  ec <- pick("End Date", "end_date", "EndDate")
  pc <- pick("Total Population (9months -15 years)",
             "Total Population (9months-15 years)", "pop_9m15y")
  rc <- pick("Target Achieved #", "Target Achieved", "reached")

  tibble(state      = state,
         district   = standardise_districts(trimws(d[[dc]]), refs[[state]],
                                            source = paste(state, "campaign")),
         start_date = parse_date_dmy(d[[sc]]),
         end_date   = parse_date_dmy(d[[ec]]),
         pop_9m15y  = clean_number(d[[pc]]),
         reached    = clean_number(d[[rc]])) %>%
    filter(!is.na(district))
}

mr <- bind_rows(lapply(STATES, read_campaign))

if (is.null(mr) || nrow(mr) == 0) {
  cat("no campaign anywhere, writing Phi = 1 for every district and biweek\n")
  out <- allu %>% crossing(grid %>% select(t_global, year, biweek)) %>%
    mutate(phi = 1, p_age = 1) %>%
    arrange(state, district, t_global)
  write_csv(out, file.path(OUT_DIR, "campaign_factors.csv"))
  quit(save = "no")
}

mr <- mr %>% mutate(uid = make_uid(state, district))

cat("\ncampaign rows:", nrow(mr), "| districts:", n_distinct(mr$uid), "\n")
cat("window:", format(min(mr$start_date, na.rm = TRUE)), "to",
    format(max(mr$end_date, na.rm = TRUE)), "\n")

bad <- mr %>% filter(is.na(start_date) | is.na(end_date) | end_date < start_date)
if (nrow(bad)) {
  cat("\nunusable dates:\n")
  print(as.data.frame(bad %>% select(state, district, start_date, end_date)))
}

# ---- coverage and the removal fraction --------------------------------------
mr <- mr %>%
  mutate(kappa_raw = reached / pop_9m15y,
         capped    = !is.na(kappa_raw) & kappa_raw > KAPPA_CAP,
         kappa     = pmin(pmax(kappa_raw, 0), KAPPA_CAP),
         mu        = E_CAMPAIGN * kappa,
         lambda    = -log(1 - pmin(mu, 0.999)))

cat("\nkappa before capping:\n"); print(summary(mr$kappa_raw))
cat("capped at", KAPPA_CAP, ":", sum(mr$capped, na.rm = TRUE), "of", nrow(mr), "\n")
cat("mu = e * kappa, e =", E_CAMPAIGN, ":\n"); print(summary(mr$mu))

# ---- day overlap weights across the delivery window -------------------------
weights <- mr %>%
  filter(!is.na(start_date), !is.na(end_date), end_date >= start_date) %>%
  select(uid, start_date, end_date) %>%
  crossing(grid %>% select(t_global, bw_start, bw_end)) %>%
  mutate(days = as.integer(pmin(end_date, bw_end) -
                           pmax(start_date, bw_start)) + 1L) %>%
  filter(days > 0) %>%
  group_by(uid) %>%
  mutate(w = days / sum(days)) %>%     # shares sum to one within a campaign
  ungroup() %>%
  select(uid, t_global, w)

cat("\nbiweeks touched per district:\n")
print(summary(weights %>% count(uid) %>% pull(n)))

# ---- the lag, then the survival factor --------------------------------------
phi_long <- weights %>%
  mutate(t_global = t_global + CAMPAIGN_LAG) %>%
  left_join(mr %>% select(uid, lambda), by = "uid") %>%
  mutate(phi = exp(-lambda * w)) %>%
  select(uid, t_global, phi)

# ---- every district and biweek, Phi = 1 where no campaign -------------------
out <- allu %>%
  crossing(grid %>% select(t_global, year, biweek)) %>%
  left_join(phi_long, by = c("uid", "t_global")) %>%
  mutate(phi   = replace_na(phi, 1),
         p_age = 1) %>%              # MR covers the whole 9m to 15y band
  arrange(state, district, t_global)

# ---- conservation ------------------------------------------------------------
# The product of Phi across a district should equal 1 - mu. This is the one
# check that would catch an error in the weights, the lag, or the join.
prod_phi <- out %>% group_by(uid) %>%
  summarise(prod_phi = prod(phi), .groups = "drop")

chk <- mr %>% select(uid, mu) %>%
  left_join(prod_phi, by = "uid") %>%
  mutate(expected = 1 - mu, err = abs(prod_phi - expected))

cat("\nconservation, product of Phi against 1 - mu:\n")
cat("  largest error:", format(max(chk$err, na.rm = TRUE), scientific = TRUE), "\n")
if (max(chk$err, na.rm = TRUE) > 1e-8) {
  cat("  failing:\n")
  print(as.data.frame(chk %>% filter(err > 1e-8) %>% head(10)))
}

cat("\nbiweeks with an effect:", sum(out$phi < 1), "\n")
by_state <- out %>% group_by(state) %>%
  summarise(districts = n_distinct(uid),
            with_campaign = n_distinct(uid[phi < 1]), .groups = "drop")
print(as.data.frame(by_state))

none <- setdiff(allu$uid, unique(out$uid[out$phi < 1]))
if (length(none)) {
  in_file <- none %in% mr$uid
  cat("\ndistricts with no campaign effect:", length(none), "\n")
  cat("  of which had a campaign row but no effect:", sum(in_file), "\n")
  cat("  and had no row at all:", sum(!in_file), "\n")
}

write_csv(out, file.path(OUT_DIR, "campaign_factors.csv"))
write_csv(mr %>% select(state, district, uid, start_date, end_date,
                        kappa_raw, kappa, mu),
          file.path(OUT_DIR, "campaign_log.csv"))
cat("\nwrote campaign_factors.csv and campaign_log.csv\n")

