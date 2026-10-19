# 29_supplementary_tables.R
# Run after 06, 06b, 20_run_all.R and the network tuning (11).

source("/.../.../.../Code/R/00_config.R")
source("/.../.../.../Code/R/01_district_names.R")

suppressPackageStartupMessages({ library(tidyverse); library(scales) })

SUP_DIR  <- file.path(TAB_DIR, "supplement")
dir.create(SUP_DIR, showWarnings = FALSE, recursive = TRUE)
TUNE_DIR <- file.path(OUT_DIR, "sfnn", "tuning")
SUP_WRITE_CSV <- FALSE   # TRUE also writes each table's data as CSV
COV_FILE <- file.path(OUT_DIR, "all_states_district_covariates.csv")

`%||%` <- function(a, b) if (is.null(a)) b else a
rd  <- function(p) read_csv(p, show_col_types = FALSE)
f1  <- function(x, d = 1) formatC(x, format = "f", digits = d)
cell <- function(md, lo, hi) {
  if (abs(md) >= 1e6) sprintf("%.2f million (%.2f to %.2f)", md / 1e6, lo / 1e6, hi / 1e6)
  else sprintf("%s (%s to %s)", comma(round(md, -3)), comma(round(lo, -3)), comma(round(hi, -3)))
}

if (!requireNamespace("gridExtra", quietly = TRUE) ||
    !requireNamespace("gtable", quietly = TRUE))
  stop("install the packages gridExtra and gtable: install.packages(c(\"gridExtra\", \"gtable\"))")

caps <- list()   # captions, by table

# A table drawn on A4 pages

draw_table <- function(df, file, title, text, landscape = FALSE, base = 7.5) {
  W <- if (landscape) 11.69 else 8.27; H <- if (landscape) 8.27 else 11.69
  margin <- 0.6
  df <- as.data.frame(lapply(df, as.character), check.names = FALSE)
  grDevices::cairo_pdf(file, width = W, height = H, onefile = TRUE)
  on.exit(grDevices::dev.off())

  theme_at <- function(sz) gridExtra::ttheme_minimal(
    base_size = sz, padding = grid::unit(c(3, 2.2), "mm"),
    core    = list(fg_params = list(hjust = 0, x = 0.04)),
    colhead = list(fg_params = list(hjust = 0, x = 0.04, fontface = "bold")))
  grob_of <- function(rows, sz) {
    g <- gridExtra::tableGrob(df[rows, , drop = FALSE], rows = NULL, theme = theme_at(sz))
    n <- ncol(g); rule <- function(y, lwd) grid::segmentsGrob(0, y, 1, y, gp = grid::gpar(lwd = lwd))
    g <- gtable::gtable_add_grob(g, rule(1, 1.1), t = 1, l = 1, r = n, name = "top")
    g <- gtable::gtable_add_grob(g, rule(0, 0.6), t = 1, l = 1, r = n, name = "mid")
    gtable::gtable_add_grob(g, rule(0, 1.1), t = nrow(g), l = 1, r = n, name = "bottom")
  }
  inch <- function(u) grid::convertUnit(u, "in", valueOnly = TRUE)

  # largest text size, from base down to 5 pt, at which the table fits the width
  grid::grid.newpage()
  sz <- base
  while (sz > 5 && inch(sum(grob_of(seq_len(nrow(df)), sz)$widths)) > W - 2 * margin)
    sz <- sz - 0.5

  cap_lines <- strwrap(text, width = floor((W - 2 * margin) * 15.5))
  row_h <- inch(grob_of(1, sz)$heights[2])
  first_room <- H - 2 * margin - 0.35 - 0.16 * length(cap_lines) - 0.5
  later_room <- H - 2 * margin - 0.6
  n1 <- max(1, floor(first_room / row_h) - 1)
  nk <- max(1, floor(later_room / row_h) - 1)
  idx <- seq_len(nrow(df)); pages <- list(head(idx, n1))
  rest <- idx[-seq_len(min(n1, length(idx)))]
  while (length(rest)) { pages[[length(pages) + 1]] <- head(rest, nk); rest <- rest[-seq_len(min(nk, length(rest)))] }

  for (p in seq_along(pages)) {
    if (p > 1) grid::grid.newpage()
    y <- 1 - margin / H
    if (p == 1) {
      grid::grid.text(title, x = margin / W, y = y, just = c("left", "top"),
                      gp = grid::gpar(fontsize = 9.5, fontface = "bold"))
      y <- y - 0.3 / H
      for (ln in cap_lines) {
        grid::grid.text(ln, x = margin / W, y = y, just = c("left", "top"), gp = grid::gpar(fontsize = 8.5))
        y <- y - 0.16 / H
      }
      y <- y - 0.12 / H
    } else {
      grid::grid.text(paste(title, "(continued)"), x = margin / W, y = y, just = c("left", "top"),
                      gp = grid::gpar(fontsize = 9, fontface = "bold"))
      y <- y - 0.35 / H
    }
    g <- grob_of(pages[[p]], sz)
    g$vp <- grid::viewport(x = margin / W, y = y, just = c("left", "top"),
                           width = sum(g$widths), height = sum(g$heights))
    grid::grid.draw(g)
  }
  invisible(sz)
}

# one table: draws its PDF, keeps its caption, and optionally writes its data
make_table <- function(df, id, title, text, landscape = FALSE) {
  stem <- file.path(SUP_DIR, gsub(" ", "_", id))
  if (SUP_WRITE_CSV) write_csv(df, paste0(stem, ".csv"))
  sz <- draw_table(df, paste0(stem, ".pdf"), paste(id, "\u2014", title), text, landscape)
  caps[[id]] <<- paste0(title, " ", text)
  cat(sprintf("  %s: %d rows, %s, text %.1f pt\n", id, nrow(df),
              if (landscape) "landscape" else "portrait", sz))
}

# ============================================================
# S1  Reconstruction by district
# ============================================================
fl <- rd(file.path(OUT_DIR, "tsir_fit_log.csv")) %>% filter(fitted)
s1 <- fl %>% arrange(state, district) %>%
  transmute(State = state, District = district,
            `Mean susceptible level` = if_else(sbar_source == "profiled", "identified",
                                               paste0("fixed at ", 100 * SBAR_FALLBACK, "%")),
            `S, % of population` = f1(100 * sbar_frac),
            `Reporting, %` = f1(100 * rho, 2),
            alpha = f1(alpha, 2),
            `R squared` = f1(rsquared, 2),
            `Biweeks with no cases, %` = f1(100 * frac_zero, 0))
make_table(s1, "S1 Table", "Reconstruction by district.",
  paste0("Main reconstruction, both states, April 2008 to March 2020. The mean ",
         "susceptible level was identified by profile likelihood where the fit reached a ",
         "minimum within the range searched, up to ", 100 * SBAR_GRID_MAX, "% of mean ",
         "population, and fixed at ", 100 * SBAR_FALLBACK, "% elsewhere. Reporting is the ",
         "estimated fraction of infections reported; alpha is the mixing exponent; R squared ",
         "is for the regression of cumulative effective births on cumulative reported cases."))

# ============================================================
# S2  Every counterfactual run
# ============================================================
runs <- c("scenarios",
          sort(list.files(OUT_DIR, pattern = "^scenarios_up_from[0-9]{4}$")),
          sort(list.files(OUT_DIR, pattern = "^scenarios_refit_up_from[0-9]{4}$")))

s2 <- map_dfr(runs, function(r) {
  fy    <- suppressWarnings(as.integer(sub(".*from", "", r)))
  refit <- grepl("^scenarios_refit", r)
  tag   <- if (refit) sub("^scenarios_refit", "", r) else ""
  m  <- rd(file.path(OUT_DIR, paste0("tsir_panel", tag, ".csv")))
  if (!is.na(fy)) {
    t0 <- (fy - BASE_YEAR) * BIWEEKS_PER_YEAR + FIRST_BIWEEK_KEPT
    m  <- m %>% filter(!(state == "Uttar Pradesh" & t_global < t0))
  }
  rr <- rd(file.path(OUT_DIR, r, "scenario_replicates.csv"))
  bs <- rd(file.path(OUT_DIR, r, "bootstrap_state.csv"))
  recon <- m %>% filter(uid %in% unique(rr$uid)) %>% group_by(state) %>%
    summarise(recon = sum(cases / rho), .groups = "drop")
  sim <- rr %>% group_by(state, b) %>% summarise(o = sum(observed), .groups = "drop")
  bs %>%
    mutate(tot = seq_MCV1_md + seq_MCV2_md + seq_camp_md) %>%
    left_join(recon, by = "state") %>%
    left_join(sim %>% filter(b > 0) %>% group_by(state) %>%
                summarise(boot = median(o), .groups = "drop"), by = "state") %>%
    left_join(sim %>% filter(b == 0) %>% select(state, zero = o), by = "state") %>%
    rowwise() %>%
    transmute(Run = case_when(is.na(fy) ~ "Main",
                              refit     ~ paste0("Refitted from ", fy),
                              TRUE      ~ paste0("From ", fy)),
              State = state,
              Period = if (is.na(fy)) "2008-2020" else paste0(fy, "-2020"),
              `In Table 2` = if_else(CF_TABLE_RUN[[state]] == r, "yes", ""),
              `Recovered, %` = f1(100 * boot / recon, 0),
              `Recovered at fitted parameters, %` = f1(100 * zero / recon, 0),
              `All vaccination` = cell(total_averted_md, total_averted_lo, total_averted_hi),
              Campaign = if (abs(marg_camp_md) < 0.5) "no campaign" else
                cell(marg_camp_md, marg_camp_lo, marg_camp_hi),
              `First / second dose / campaign, %` = paste(f1(100 * seq_MCV1_md / tot),
                                                        f1(100 * seq_MCV2_md / tot),
                                                        f1(100 * seq_camp_md / tot), sep = " / ")) %>%
    ungroup()
})
make_table(s2, "S2 Table", "Every counterfactual run: calibration and estimates.",
  paste0("Recovered is the simulated infections in the observed scenario as a percentage of ",
         "the reconstructed infections (reported cases divided by the reporting fraction): the ",
         "median across bootstrap replicates, and at the fitted parameters. Infections averted ",
         "are medians with 95% intervals across 1,000 bootstrap replicates. The last column ",
         "divides the infections averted between the first dose, the second dose and the ",
         "campaign, added in that order. Uttar Pradesh reported few cases before late 2012, so ",
         "its counterfactual was also run from April 2012 and April 2013, and with the ",
         "reconstruction refitted from April 2013. The runs reported in Table 2 are marked."), landscape = TRUE)

# ============================================================
# S3  Correlation and variance inflation of the static covariates
# ============================================================
COVS <- c(measles_zero = "Children never given MCV",
          penta_zero = "Children never given pentavalent",
          full_imm = "Full immunisation",
          women10plus = "Women with 10 or more years of schooling",
          poorest = "Women in poorest wealth quintile",
          institutional_birth = "Institutional birth",
          maternal_immunisation = "Last birth protected against tetanus",
          sc_st = "Scheduled Caste and Tribe share",
          hhsize = "Mean household size",
          urbanization = "Women living in urban areas",
          stunting = "Child stunting",
          vitamin_a = "Vitamin A supplementation",
          early_breastfeeding = "Early initiation of breastfeeding")
PATH <- c(rep("Supply", 3), rep("Access", 5), rep("Contact", 2), rep("Nutrition", 3))

if (file.exists(COV_FILE)) {
  cv <- rd(COV_FILE) %>% filter(tolower(state) %in% tolower(STATES)) %>%
    select(all_of(names(COVS)))
  cm <- cor(cv, use = "pairwise.complete.obs")
  # VIF of each indicator, from the regression of it on the other twelve
  vif <- sapply(names(COVS), function(v) {
    r2 <- summary(lm(reformulate(setdiff(names(COVS), v), v), data = cv))$r.squared
    1 / (1 - r2)
  })
  lower <- sapply(seq_along(COVS), function(j)
    sapply(seq_along(COVS), function(i)
      if (i > j) f1(cm[i, j], 2) else if (i == j) "1" else ""))
  s3 <- tibble(` ` = seq_along(COVS), Indicator = unname(COVS), Pathway = PATH) %>%
    bind_cols(as_tibble(lower, .name_repair = ~ as.character(seq_along(COVS)))) %>%
    mutate(VIF = f1(vif, 1))
  make_table(s3, "S3 Table", "Correlation and variance inflation of the static covariates.",
    paste0("Pearson correlations between the thirteen NFHS-5 indicators across the ",
           nrow(cv), " districts (lower triangle; column numbers refer to the numbered rows), ",
           "and each indicator's variance inflation factor (VIF). Pathways: Supply, ",
           "susceptible supply; Access, service access; Contact, contact structure; ",
           "Nutrition, nutrition and outreach."), landscape = TRUE)
} else {
  cat("  S3 Table skipped: no covariate file at", COV_FILE, "\n")
}

# ============================================================
# S4  Selected hyperparameters by horizon and variant
# ============================================================
VARIANTS <- c(no_nfhs = "Without covariates", with_nfhs = "With covariates")
s4 <- map_dfr(names(VARIANTS), function(v) map_dfr(c(1, 4, 6, 12), function(k) {
  p <- file.path(TUNE_DIR, v, paste0("tuning_k", k, ".csv"))
  if (!file.exists(p)) return(NULL)
  b <- rd(p) %>% arrange(val_mse) %>% slice(1)
  tibble(Model = VARIANTS[[v]], `Horizon, biweeks` = k,
         `Lags, biweeks` = b$t_lag, `Hidden layers` = b$num_hidden_layers,
         `Units per layer` = b$hidden_dim, `Learning rate` = formatC(b$lr, format = "e", digits = 1),
         `Weight decay` = formatC(b$weight_decay, format = "e", digits = 1),
         Epochs = b$epoch, `Validation MSE` = f1(b$val_mse, 3))
}))
if (nrow(s4)) {
  make_table(s4, "S4 Table", "Selected hyperparameters by horizon and model variant.",
    paste0("For each horizon and variant, the configuration with the lowest validation error ",
           "(April 2016 to March 2017) among the trials of the hyperparameter search."))
} else {
  cat("  S4 Table skipped: no tuning results in", TUNE_DIR, "\n")
}

# ============================================================
# S5  Validation checks on the runs reported in Table 2
# ============================================================
s5 <- map_dfr(names(CF_TABLE_RUN), function(st) {
  v <- rd(file.path(OUT_DIR, CF_TABLE_RUN[[st]], "scenario_validation.csv")) %>%
    filter(state == st)
  bad <- v %>% filter(any_fail)
  if (!nrow(bad))
    return(tibble(State = st, District = paste0("none of ", nrow(v), " districts"),
                  `Check failed` = "", `No-vaccination ratio` = ""))
  bad %>% arrange(district) %>%
    transmute(State = st, District = district,
              `Check failed` = case_when(fail_order & fail_plausible ~ "ordering; ratio above 20",
                                         fail_order ~ "ordering of scenarios",
                                         fail_plausible ~ "ratio above 20",
                                         TRUE ~ "other"),
              `No-vaccination ratio` = f1(ratio, 1))
})
make_table(s5, "S5 Table", "Validation checks on the runs reported in Table 2.",
  paste0("Districts failing at least one check: the ordering of scenarios (more vaccination ",
         "never producing more infections) inverted in some biweeks, or the no-vaccination ",
         "burden exceeding twenty times the observed scenario. Flagged districts are included ",
         "in all totals. West Bengal, ", CF_TABLE_PERIOD[["West Bengal"]], "; Uttar Pradesh, ",
         CF_TABLE_PERIOD[["Uttar Pradesh"]], "."))

# ============================================================
# S6  District name harmonisation and boundary changes
# ============================================================
refs <- load_all_references(STATES, file_of)
state_of <- function(nm) {
  hit <- names(refs)[sapply(refs, function(r) nm %in% r)]
  if (length(hit)) hit[1] else ""
}
s6 <- tibble(`Name in a source file` = names(DISTRICT_ALIASES),
             `Name used` = unname(DISTRICT_ALIASES)) %>%
  mutate(State = map_chr(`Name used`, state_of), Note = "") %>%
  select(State, `Name in a source file`, `Name used`, Note) %>%
  arrange(State, `Name used`) %>%
  bind_rows(tibble(
    State = "West Bengal",
    `Name in a source file` = c("Purba Barddhaman, Paschim Barddhaman", "Alipurduar"),
    `Name used` = c("Barddhaman (maps only)", "Jalpaiguri (maps only)"),
    Note = c("divided in 2017; one boundary in the shapefile, so the two are pooled on the maps",
             "created from Jalpaiguri in 2014; no separate data, shown with Jalpaiguri's estimate")))
make_table(s6, "S6 Table", "District name harmonisation and boundary changes.",
  paste0("Names differing between the case, population, coverage, survey and boundary files ",
         "were matched to one reference list per state. Names are shown as found after ",
         "removing spaces, punctuation and capitals. District definitions follow the case ",
         "data; the last two rows record where the boundary data differ from them."),
  landscape = TRUE)
# ============================================================
# the document and the captions
# ============================================================
writeLines(paste0(names(caps), ": ", unlist(caps), "\n"),
           file.path(SUP_DIR, "S_tables_captions.txt"))
cat("\n", length(caps), " tables written to ", SUP_DIR, "\n", sep = "")
