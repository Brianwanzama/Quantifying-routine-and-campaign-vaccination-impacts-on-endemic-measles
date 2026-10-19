# 21_fig_scenarios.R
# Tables and figures from the counterfactual analysis.
#Each state is read from the run set in 00_config.R (CF_TABLE_RUN):
# West Bengal from the main run, April 2008 to March 2020; Uttar Pradesh from
# the reconstruction refitted from April 2013, because its reporting before late
# 2012 was too sparse to anchor a simulation (Fig 1). The two states therefore
# cover different periods, and every caption says so.

source("/.../.../.../Code/R/00_config.R")
source("/.../.../.../Code/R/01_district_names.R")

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(scales); library(patchwork)
})

SHAPE <- Sys.getenv("MS_SHAPEFILE",
                    file.path(dirname(DATA_DIR), "shapefiles", "gadm36_IND_2_sp.rds"))

# ---- reading each state from its own run ------------------------------------------
rd_run <- function(st, f, req = TRUE) {
  p <- file.path(OUT_DIR, CF_TABLE_RUN[[st]], f)
  if (!file.exists(p)) {
    if (req) stop("no ", f, " for ", st, " in ", CF_TABLE_RUN[[st]],
                  ": run 20_run_all.R first, or check CF_TABLE_RUN in 00_config.R")
    return(NULL)
  }
  read_csv(p, show_col_types = FALSE) %>% filter(state == st)
}
both <- function(f, req = TRUE) bind_rows(lapply(STATES, rd_run, f = f, req = req))

bd    <- both("bootstrap_district.csv")
bs    <- both("bootstrap_state.csv")
reps  <- both("scenario_replicates.csv") %>% filter(b > 0)
fyrep <- both("scenario_year_replicates.csv", req = FALSE)
# The series has 313 biweeks, so its last is the first biweek of fiscal 2020
# (April 2020) on its own; one biweek is not a year, so it is left out of the
# yearly summaries. Fiscal years are numbered by the year they start in.
if (!is.null(fyrep) && nrow(fyrep)) fyrep <- fyrep %>% filter(b > 0, fy < 2020)
valid <- both("scenario_validation.csv", req = FALSE)

period_note <- paste0(names(CF_TABLE_PERIOD), ", ", CF_TABLE_PERIOD, collapse = "; ")

# ---- presentation -------------------------------------------------------------------
COMP <- c("MCV1", "MCV2", "Campaign")
pal  <- setNames(c("#1F78B4", "#33A02C", "#E31A1C"), COMP)
pal_state <- setNames(c("#2166AC", "#B2182B"), STATES)

# fiscal years written out, since 2018 alone reads as a calendar year
fy_label <- function(x) paste0(x, "\u2013", substr(x + 1, 3, 4))

theme_pol <- theme_classic(base_size = 10) +
  theme(axis.title = element_text(size = 9.5, colour = "grey15"),
        axis.text  = element_text(size = 8.5, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line  = element_line(colour = "grey20", linewidth = 0.4),
        axis.ticks = element_line(colour = "grey20", linewidth = 0.35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 9),
        plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0),
        plot.margin = margin(5, 9, 5, 6))
theme_set(theme_pol)

# Every figure is saved as a vector PDF for LaTeX, and as PNG and LZW-compressed
# TIFF at FIG_DPI for slides and the journal upload, written with ragg when
# available for cleaner anti-aliasing.
FIG_DPI  <- as.integer(Sys.getenv("MS_FIG_DPI", "600"))
has_ragg <- requireNamespace("ragg", quietly = TRUE)
sv <- function(p, nm, w, h) {
  f <- file.path(FIG_DIR, paste0("fig_cf_", nm))
  ggsave(paste0(f, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(f, ".png"), p, width = w, height = h, dpi = FIG_DPI, bg = "white",
         device = if (has_ragg) ragg::agg_png else "png")
  if (has_ragg) {
    ggsave(paste0(f, ".tif"), p, width = w, height = h, dpi = FIG_DPI, bg = "white",
           device = ragg::agg_tiff, compression = "lzw")
  } else {
    ggsave(paste0(f, ".tif"), p, width = w, height = h, dpi = FIG_DPI, bg = "white",
           device = "tiff", compression = "lzw")
  }
  cat("  saved fig_cf_", nm, " as pdf, png and tif\n", sep = "")
}

# numbers as the manuscript writes them: millions to two decimals, otherwise to
# the nearest thousand; the interval takes the estimate's unit
cell <- function(md, lo, hi) {
  if (abs(md) >= 1e6) sprintf("%.2f million (%.2f to %.2f)", md / 1e6, lo / 1e6, hi / 1e6)
  else sprintf("%s (%s to %s)", comma(round(md, -3)), comma(round(lo, -3)), comma(round(hi, -3)))
}

# ============================================================
# Table 2: infections averted
# ============================================================
rows <- c(marg_MCV1     = "First and second doses",
          marg_MCV2     = "Second dose alone",
          marg_camp     = "Campaign",
          total_averted = "All vaccination")

tab_states <- names(CF_TABLE_RUN)
tab <- map_dfr(tab_states, function(st) {
  b1 <- bs %>% filter(state == st)
  stopifnot(nrow(b1) == 1)
  tibble(state = st, row = unname(rows),
         value = map_chr(names(rows), function(v) {
           md <- b1[[paste0(v, "_md")]]
           if (v == "marg_camp" && abs(md) < 0.5) return("no campaign")
           cell(md, b1[[paste0(v, "_lo")]], b1[[paste0(v, "_hi")]])
         }))
}) %>% pivot_wider(names_from = state, values_from = value)

write_csv(tab, file.path(TAB_DIR, "table2_averted.csv"))
cat("\n---- Table 2: infections averted (median, 95% interval) ----\n")
print(as.data.frame(tab))

# Open table: rules above and below the header, before the total, and at the end.
body <- apply(tab, 1, function(r) paste0(paste(gsub(" to ", "--", r), collapse = " & "), " \\\\"))
tex <- c("\\begin{tabular}{lll}", "\\hline",
         paste0("Removed & ", paste0(tab_states, ", ", CF_TABLE_PERIOD[tab_states], collapse = " & "),
                " \\\\"), "\\hline",
         head(body, -1), "\\hline", tail(body, 1), "\\hline",
         "\\end{tabular}")
writeLines(tex, file.path(TAB_DIR, "table2_averted.tex"))

# ============================================================
# S Table: calibration of every counterfactual run
# ============================================================
runs <- c("scenarios",
          sort(list.files(OUT_DIR, pattern = "^scenarios_up_from[0-9]{4}$")),
          sort(list.files(OUT_DIR, pattern = "^scenarios_refit_up_from[0-9]{4}$")))

calib <- map_dfr(runs, function(r) {
  fy    <- suppressWarnings(as.integer(sub(".*from", "", r)))
  refit <- grepl("^scenarios_refit", r)
  tag   <- if (refit) sub("^scenarios_refit", "", r) else ""
  m <- read_csv(file.path(OUT_DIR, paste0("tsir_panel", tag, ".csv")), show_col_types = FALSE)
  if (!is.na(fy)) {
    t0 <- (fy - BASE_YEAR) * BIWEEKS_PER_YEAR + FIRST_BIWEEK_KEPT
    m <- m %>% filter(!(state == "Uttar Pradesh" & t_global < t0))
  }
  rr <- read_csv(file.path(OUT_DIR, r, "scenario_replicates.csv"), show_col_types = FALSE)
  m %>% filter(uid %in% unique(rr$uid)) %>%
    group_by(state) %>% summarise(reconstructed = sum(cases / rho), .groups = "drop") %>%
    left_join(rr %>% filter(b > 0) %>% group_by(state, b) %>%
                summarise(o = sum(observed), .groups = "drop") %>%
                group_by(state) %>% summarise(simulated = median(o), .groups = "drop"),
              by = "state") %>%
    left_join(rr %>% filter(b == 0) %>% group_by(state) %>%
                summarise(no_noise = sum(observed), .groups = "drop"), by = "state") %>%
    transmute(run = case_when(is.na(fy) ~ "Main, from 2008",
                              refit     ~ paste0("Refitted, from ", fy),
                              TRUE      ~ paste0("From ", fy, ", main reconstruction")),
              state, reconstructed = round(reconstructed), simulated = round(simulated),
              recovered_pct = round(100 * simulated / reconstructed),
              no_noise_pct  = round(100 * no_noise / reconstructed))
})
write_csv(calib, file.path(TAB_DIR, "S_table_calibration.csv"))
cat("\n---- S Table: calibration of the counterfactual runs ----\n")
print(as.data.frame(calib %>% mutate(across(c(reconstructed, simulated), comma))))

# ============================================================
# district medians, used by the map and the district figures
# ============================================================
dist <- reps %>%
  group_by(state, district, uid, b) %>%
  summarise(none = sum(none), observed = sum(observed), .groups = "drop") %>%
  group_by(state, district, uid) %>%
  summarise(none_md = median(none), observed_md = median(observed),
            prevented_pct = median(100 * (1 - observed / pmax(none, 1))),
            .groups = "drop") %>%
  left_join(bd %>% select(uid, mean_pop, starts_with("total_averted"),
                          starts_with("marg_camp"), seq_MCV1_md, seq_MCV2_md, seq_camp_md),
            by = "uid")

# ============================================================
# 1. the map: infections averted per 1,000
# ============================================================
refs <- load_all_references(STATES, file_of)

shp <- st_as_sf(readRDS(SHAPE)) %>%
  filter(NAME_1 %in% STATES) %>%
  group_by(NAME_1) %>%
  group_modify(~ {
    ref <- refs[[.y$NAME_1]]
    .x %>% mutate(district = standardise_districts(
      NAME_2, ref, source = paste(.y$NAME_1, "shapefile")))
  }) %>%
  ungroup() %>% st_as_sf() %>%
  mutate(state = NAME_1,
         district = if_else(is.na(district), trimws(NAME_2), district),
         uid = make_uid(state, district)) %>%
  select(state, uid, geometry)

rate <- dist %>%
  transmute(state, district, uid, mean_pop,
            md = 1000 * total_averted_md / mean_pop,
            lo = 1000 * total_averted_lo / mean_pop,
            hi = 1000 * total_averted_hi / mean_pop,
            camp_md = 1000 * marg_camp_md / mean_pop,
            camp_lo = 1000 * marg_camp_lo / mean_pop,
            camp_hi = 1000 * marg_camp_hi / mean_pop)

POOL <- list(`West Bengal__Barddhaman` = c("West Bengal__Purba Barddhaman",
                                           "West Bengal__Paschim Barddhaman"))
INHERIT <- c(`West Bengal__Alipurduar` = "West Bengal__Jalpaiguri")

md_map <- rate
for (poly in names(POOL)) {
  have <- intersect(POOL[[poly]], md_map$uid)
  if (!poly %in% shp$uid || !length(have)) next
  bits <- strsplit(poly, UID_SEP, fixed = TRUE)[[1]]
  p <- md_map %>% filter(uid %in% have)
  w <- sum(p$mean_pop)
  md_map <- bind_rows(md_map %>% filter(!uid %in% have),
                      tibble(state = bits[1], district = bits[2], uid = poly, mean_pop = w,
                             md = sum(p$md * p$mean_pop) / w, lo = sum(p$lo * p$mean_pop) / w,
                             hi = sum(p$hi * p$mean_pop) / w,
                             camp_md = sum(p$camp_md * p$mean_pop) / w,
                             camp_lo = sum(p$camp_lo * p$mean_pop) / w,
                             camp_hi = sum(p$camp_hi * p$mean_pop) / w))
}
for (child in names(INHERIT)) {
  parent <- INHERIT[[child]]
  if (!child %in% shp$uid || !parent %in% md_map$uid || child %in% md_map$uid) next
  bits <- strsplit(child, UID_SEP, fixed = TRUE)[[1]]
  md_map <- bind_rows(md_map, md_map %>% filter(uid == parent) %>%
                        mutate(state = bits[1], district = bits[2], uid = child))
}

# Infections averted per 1,000, both states stacked as a and b at the same scale,
# in the style of the susceptibility map in script 18: four categories set in the
# config, traffic-light colours, every district shaded and named, one legend
# drawn as a stacked bar.
MAP_W  <- as.numeric(Sys.getenv("MS_FIG_W", "5.5"))
MAP_H  <- as.numeric(Sys.getenv("MS_FIG_H", "6.5"))
BORDER <- as.numeric(Sys.getenv("MS_BORDER", "0.25"))

mapdat <- shp %>% left_join(md_map %>% select(uid, md), by = "uid") %>%
  mutate(district = sub(paste0("^.*", UID_SEP), "", uid))

# Both states, Uttar Pradesh above West Bengal, sharing one legend. The four
# categories are set by hand in the config, CF_MAP_BREAKS, as five values from
# -Inf to Inf at equal intervals, for example
#   CF_MAP_BREAKS <- c(-Inf, 25, 50, 75, Inf)
# The district values of each state are summarised first, to help choose them.
cat("\ninfections averted per 1,000, district values by state:\n")
print(as.data.frame(st_drop_geometry(mapdat) %>% filter(!is.na(md)) %>%
                      group_by(state) %>%
                      summarise(min = round(min(md), 1), p25 = round(quantile(md, .25), 1),
                                median = round(median(md), 1), p75 = round(quantile(md, .75), 1),
                                p90 = round(quantile(md, .9), 1), max = round(max(md), 1), .groups = "drop")))

if (!exists("CF_MAP_BREAKS")) {
  cat("  map not drawn: set CF_MAP_BREAKS in 00_config.R, five values from -Inf to Inf,",
      "for example c(-Inf, 25, 50, 75, Inf); the summary above shows the range\n")
} else {
  brk <- CF_MAP_BREAKS
  if (length(brk) != 5 || !is.infinite(brk[1]) || !is.infinite(brk[5]) || is.unsorted(brk))
    stop("CF_MAP_BREAKS must be five increasing values from -Inf to Inf, giving four categories")
  nb    <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
  inner <- brk[is.finite(brk)]
  map_labels <- c(paste0("<", nb(inner[1])),
                  paste0(nb(head(inner, -1)), "\u2013", nb(inner[-1])),
                  paste0("\u2265", nb(tail(inner, 1))))
  mapdat <- mapdat %>% mutate(band = cut(md, brk, labels = map_labels, right = FALSE))
  cat("\ncategories:", paste(map_labels, collapse = " | "), "\n")
  print(as.data.frame(mapdat %>% st_drop_geometry() %>%
                        count(state, band, .drop = FALSE) %>% pivot_wider(names_from = band, values_from = n)))
  

  MAP_PAL <- c("#D55E00", "#E69F00", "#F0E442", "#009E73")
  MAP_TXT <- c("white", "black", "black", "white")
  
  # Uniform scale: both panels share the same east-west extent, centred on each
  # state, and their heights follow each state's north-south extent, so a
  # kilometre is the same length in both panels.
  bb_of  <- lapply(setNames(STATES, STATES), function(st) st_bbox(mapdat %>% filter(state == st)))
  x_span <- max(sapply(bb_of, function(b) unname(b["xmax"] - b["xmin"])))
  lim_of <- function(st) {
    b <- bb_of[[st]]; cx <- unname(b["xmin"] + b["xmax"]) / 2
    list(x = c(cx - x_span / 2, cx + x_span / 2), y = unname(c(b["ymin"], b["ymax"])))
  }
  panel_heights <- sapply(STATES, function(st) unname(diff(lim_of(st)$y)))
  # legend blocks touching, as one stacked bar (the setting exists from ggplot2 3.5)
  key_gap <- if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.key.spacing.y = unit(0, "pt")) else theme()
  
  map_one <- function(st, tag) {
    d <- mapdat %>% filter(state == st) %>%
      mutate(txt = coalesce(MAP_TXT[as.integer(band)], "black"))
    ggplot(d) +
      geom_sf(aes(fill = band), colour = "white", linewidth = BORDER) +
      geom_sf_text(aes(label = district, colour = txt), size = MAP_LABEL_SIZE,
                   fontface = "bold", check_overlap = TRUE) +
      scale_fill_manual(values = setNames(MAP_PAL, map_labels), drop = FALSE,
                        na.value = "grey88", name = "Infections averted\nper 1,000",
                        guide = guide_legend(ncol = 1, reverse = TRUE,
                                             title.position = "top", title.hjust = 0,
                                             keyheight = unit(16, "pt"),
                                             keywidth = unit(14, "pt"),
                                             byrow = TRUE,
                                             override.aes = list(colour = NA))) +
      scale_colour_identity() +
      labs(tag = tag, title = st) +
      coord_sf(xlim = lim_of(st)$x, ylim = lim_of(st)$y, expand = FALSE) +
      theme_void(base_size = 8) +
      theme(plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
            plot.tag          = element_text(face = "bold", size = 10),
            plot.tag.position = c(0.02, 0.98),
            legend.position   = "right",
            legend.title      = element_text(size = 8, face = "bold"),
            legend.text       = element_text(size = 8),
            plot.margin       = margin(2, 2, 2, 2))
  }
  
  # One legend, Uttar Pradesh's, for both panels: the categories and colours are
  # the same in each, so it describes West Bengal too. West Bengal's own legend is
  # switched off, since the two would otherwise both be drawn whenever they differ
  # in any detail, such as one state having a district with no data.
  p1 <- (map_one(STATES[1], "a") /
           (map_one(STATES[2], "b") + guides(fill = "none"))) +
    plot_layout(nrow = 2, ncol = 1, heights = panel_heights, guides = "collect") &
    theme(legend.position = "right") & key_gap
  sv(p1, "map", MAP_W, MAP_H)
}

# ============================================================
# 2. the scenarios for two districts, biweekly
# ============================================================

stack_file <- if (CF_STACK_SOURCE == "point") "scenario_series_point.csv" else
  "scenario_series_median.csv"
rd_stack <- function(st) {
  p <- file.path(OUT_DIR, CF_STACK_RUN[[st]], stack_file)
  if (!file.exists(p)) return(NULL)
  x <- read_csv(p, show_col_types = FALSE) %>% filter(state == st)
  if (CF_STACK_SOURCE == "point")
    x <- x %>% mutate(seq_MCV1 = none - mcv1_only, seq_MCV2 = mcv1_only - routine,
                      seq_camp = routine - observed)
  x
}
smed <- bind_rows(lapply(STATES, rd_stack))
stack_period <- setNames(ifelse(CF_STACK_RUN[STATES] == "scenarios", "April 2008 to March 2020",
                                CF_TABLE_PERIOD[STATES]), STATES)
sel <- tibble(state = STATES, district = unname(FIG_DISTRICTS[STATES])) %>%
  mutate(uid = make_uid(state, district))
missing <- setdiff(sel$uid, unique(smed$uid))
if (length(missing)) {
  cat("  scenario figure not drawn: no", stack_file, "for",
      paste(sub(UID_SEP, ", ", missing, fixed = TRUE), collapse = " and "),
      "\n  rerun script 20 for that state's run (CF_STACK_RUN), or check FIG_DISTRICTS\n")
} else {
  
  stack <- smed %>% filter(uid %in% sel$uid) %>%
    transmute(panel = paste0(district, ", ", state),
              time = year + (biweek - 1) / BIWEEKS_PER_YEAR,
              `still occurred`     = pmax(observed, 0),
              `MCV1 prevented`     = pmax(seq_MCV1, 0),
              `MCV2 prevented`     = pmax(seq_MCV2, 0),
              `campaign prevented` = pmax(seq_camp, 0)) %>%
    pivot_longer(-c(panel, time), names_to = "band", values_to = "I") %>%
    mutate(band = factor(band, levels = c("campaign prevented", "MCV2 prevented",
                                          "MCV1 prevented", "still occurred")),
           # Uttar Pradesh on top, as in Fig 3, rather than alphabetical by district
           panel = factor(panel, levels = paste0(sel$district, ", ", sel$state)))
  
  pal_band <- c(`still occurred` = "grey55", `MCV1 prevented` = "#1F78B4",
                `MCV2 prevented` = "#33A02C", `campaign prevented` = "#E31A1C")
  
  if (nrow(stack)) {
    p2 <- ggplot(stack, aes(time, I, fill = band)) +
      geom_area(position = "stack", colour = NA, alpha = 0.92) +
      # one district above the other, each panel with its own axes, as in Fig 3
      facet_wrap(~ panel, ncol = 1, scales = "free") +
      scale_fill_manual(values = pal_band,
                        breaks = c("still occurred", "MCV1 prevented",
                                   "MCV2 prevented", "campaign prevented")) +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.03))) +
      scale_x_continuous(breaks = seq(2008, 2020, 2)) +
      guides(fill = guide_legend(nrow = 1)) +
      labs(x = NULL, y = "Infections per biweek") +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
            strip.background = element_blank(),                 # no box behind the district name
            strip.text  = element_text(face = "bold", size = 9, hjust = 0),
            legend.key.size = unit(8, "pt"))
    sv(p2, "scenarios_stack", 7.5, 5.5)   # PLOS: width <= 7.5 in
  }
}

if (!is.null(fyrep) && nrow(fyrep)) {
  # ============================================================
  # 3. infections averted per fiscal year, by component
  # ============================================================
  seqyr <- fyrep %>%
    mutate(MCV1 = none - mcv1_only, MCV2 = mcv1_only - routine,
           Campaign = routine - observed) %>%
    group_by(state, fy, b) %>%
    summarise(across(all_of(COMP), sum), .groups = "drop") %>%
    group_by(state, fy) %>%
    summarise(across(all_of(COMP), median), .groups = "drop") %>%
    pivot_longer(all_of(COMP), names_to = "component", values_to = "averted") %>%
    mutate(averted = pmax(averted, 0),
           component = factor(component, levels = rev(COMP)))
  
  # totals written above each column, so the yearly figure can be read directly
  totals <- seqyr %>%
    group_by(state, fy) %>%
    summarise(total = sum(averted), .groups = "drop") %>%
    mutate(lab = ifelse(total >= 1e6,
                        label_number(scale = 1e-6, accuracy = 0.01, suffix = "M")(total),
                        label_number(scale = 1e-3, accuracy = 1, suffix = "k")(total)))
  fy_all <- sort(unique(seqyr$fy))   # shared year axis, so the two states line up
  
  p3 <- ggplot(seqyr, aes(factor(fy, levels = fy_all), averted, fill = component)) +
    geom_col(width = 0.6) +
    geom_text(data = totals,
              aes(factor(fy, levels = fy_all), total, label = lab),
              inherit.aes = FALSE, vjust = -0.4, size = 8 / .pt, colour = "grey15") +
    # each panel draws its own year axis; drop = FALSE keeps every year on both,
    # so the columns of the two states still line up
    facet_wrap(~ state, ncol = 1, scales = "free", drop = FALSE) +
    scale_fill_manual(values = pal, breaks = COMP,
                      labels = c(MCV1 = "First dose", MCV2 = "Second dose",
                                 Campaign = "Campaign")) +
    scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_x_discrete(drop = FALSE, labels = function(x) fy_label(as.integer(x))) +
    labs(x = "Fiscal year", y = "Infections averted") +
    theme(strip.background = element_blank(),                 # no box behind the state name
          strip.text  = element_text(face = "bold", size = 9, hjust = 0),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.key.size = unit(8, "pt"))
  sv(p3, "by_year", 7.5, 5.5)   # PLOS: width <= 7.5 in
  
  cat("\ninfections averted per fiscal year, by component (medians):\n")
  print(as.data.frame(seqyr %>%
                        pivot_wider(names_from = component, values_from = averted) %>%
                        transmute(state, fy = fy_label(fy), across(all_of(COMP), ~ comma(round(.x))))))
}

# ============================================================
# 4. the campaign, where it acted
# ============================================================
camp <- md_map %>% filter(state == "Uttar Pradesh", is.finite(camp_md)) %>%
  arrange(desc(camp_md)) %>% slice_head(n = 20) %>%
  mutate(district = fct_reorder(district, camp_md))

if (nrow(camp)) {
  p4 <- ggplot(camp, aes(camp_md, district)) +
    geom_linerange(aes(xmin = camp_lo, xmax = camp_hi), linewidth = 0.8,
                   colour = "#E31A1C", alpha = 0.45) +
    geom_point(size = 2, colour = "#E31A1C") +
    geom_text(aes(x = camp_hi, label = sprintf("%.1f", camp_md)),
              hjust = -0.3, size = 2.6) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    labs(x = "Infections averted per 1,000 by the campaign", y = NULL) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  sv(p4, "campaign", 6.5, 5)
}

# ============================================================
# 5. the burden that remains, and what removed the rest
# ============================================================

remain <- dist %>%
  mutate(remaining = 1000 * observed_md / mean_pop,
         counterfactual = 1000 * none_md / mean_pop)

build_top <- function(st, n_top = 20) {
  remain %>% filter(state == st) %>%
    arrange(desc(remaining)) %>% slice_head(n = n_top) %>%
    mutate(seq_sum = pmax(seq_MCV1_md + seq_MCV2_md + seq_camp_md, 1),
           averted_1000 = pmax(counterfactual - remaining, 0),
           MCV1 = averted_1000 * pmax(seq_MCV1_md, 0) / seq_sum,
           MCV2 = averted_1000 * pmax(seq_MCV2_md, 0) / seq_sum,
           Campaign = averted_1000 * pmax(seq_camp_md, 0) / seq_sum,
           district = fct_reorder(district, counterfactual))
}

pal_stack <- c(MCV1 = "#1F78B4", MCV2 = "#33A02C", Campaign = "#E31A1C",
               `still occurred` = "grey70")

for (st in STATES) {
  top <- build_top(st)
  if (!nrow(top)) next
  stk <- top %>%
    select(district, MCV1, MCV2, Campaign, `still occurred` = remaining) %>%
    pivot_longer(-district, names_to = "part", values_to = "per1000") %>%
    mutate(part = factor(part, levels = c("MCV1", "MCV2", "Campaign", "still occurred")))
  p <- ggplot(stk, aes(per1000, district, fill = part)) +
    geom_col(width = 0.72) +
    geom_text(data = top, aes(x = counterfactual, y = district,
                              label = sprintf("%.0f%% prevented", prevented_pct)),
              inherit.aes = FALSE, hjust = -0.08, size = 2.4) +
    scale_fill_manual(values = pal_stack, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.24))) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(x = "Infections per 1,000 without vaccination", y = NULL) +
    theme(axis.text.y = element_text(size = 7.5),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
  sv(p, paste0("remaining_top20_", gsub("[^a-z]+", "_", tolower(st))), 7.5, 5.5)
}

# ============================================================
# 6. what it bought against what was at stake
# ============================================================
burden <- dist %>%
  mutate(cf = 1000 * none_md / mean_pop, av = 1000 * total_averted_md / mean_pop) %>%
  filter(cf > 0, av > 0)

p5 <- ggplot(burden, aes(cf, av, colour = state, size = mean_pop)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = pal_state) +
  scale_size_continuous(range = c(0.9, 4), guide = "none") +
  scale_x_log10(labels = label_number(accuracy = 1)) +
  scale_y_log10(labels = label_number(accuracy = 1)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Infections without vaccination, per 1,000, log scale",
       y = "Infections averted, per 1,000, log scale") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
sv(p5, "averted_vs_burden", 6.5, 4.5)

if (!is.null(valid) && nrow(valid) && any(valid$any_fail)) {
  cat("\n", sum(valid$any_fail), "districts in these runs failed a validation check and",
      "are included; see scenario_validation.csv in each run\n")
}
# ============================================================
# captions, for the manuscript rather than the figures
# ============================================================
p_up <- CF_TABLE_PERIOD[["Uttar Pradesh"]]; p_wb <- CF_TABLE_PERIOD[["West Bengal"]]
med  <- if (CF_PROCESS_NOISE) "Medians across 1,000 bootstrap replicates" else
  paste("Medians across 1,000 bootstrap replicates, each a deterministic simulation",
        "with refitted seasonal transmission rates and alpha")
caps <- c(
  fig_cf_map = paste0("Infections averted per 1,000 population by vaccination, by district. ",
                      "(a) Uttar Pradesh, ", p_up, "; (b) West Bengal, ", p_wb, ". ", med,
                      ", in four categories of equal width shared by both panels; both states are drawn at ",
                      "the same scale. Barddhaman appears undivided because the boundary data predate its ",
                      "2017 division; Alipurduar, created from Jalpaiguri in 2014, is shown with ",
                      "Jalpaiguri's estimate."),
  fig_cf_scenarios_stack = paste0("Simulated epidemics in ", FIG_DISTRICTS[["Uttar Pradesh"]],
                                  ", Uttar Pradesh (", stack_period[["Uttar Pradesh"]], ") and ", FIG_DISTRICTS[["West Bengal"]],
                                  ", West Bengal (", stack_period[["West Bengal"]],
                                  "). The full height is the epidemic without vaccination; the grey band is what occurred, ",
                                  "and the coloured bands are what the first dose, the second dose and the campaign each ",
                                  "prevented, added in that order. ",
                                  if (CF_STACK_SOURCE == "point") paste0("The deterministic simulation from the fitted ",
                                                                         "model, without random variation; estimates in the text are medians across 1,000 ",
                                                                         "bootstrap replicates.") else paste0(med, ".")),
  fig_cf_by_year = paste0("Infections averted per fiscal year by the first dose, the ",
                          "second dose and the campaign, added in that order so that the three parts sum to ",
                          "the total. ", med, ". Totals with 95% intervals are in Table 2. Uttar Pradesh, ",
                          p_up, "; West Bengal, ", p_wb, "."),
  fig_cf_campaign = paste0("Infections averted per 1,000 population by the measles-rubella ",
                           "campaign in the twenty Uttar Pradesh districts with the largest effect, ", p_up,
                           ". Points are medians and bars 95% intervals across 1,000 bootstrap replicates."),
  fig_cf_remaining_top20 = paste0("Infections per 1,000 without vaccination in the twenty ",
                                  "districts with the largest remaining burden, in each state. The coloured part of each ",
                                  "bar is what was prevented, divided between the components in proportion to their ",
                                  "sequential contributions; the grey part is what occurred. The percentage is the share ",
                                  "prevented. ", med, ". Uttar Pradesh, ", p_up, "; West Bengal, ", p_wb, "."),
  fig_cf_averted_vs_burden = paste0("Infections averted against infections without ",
                                    "vaccination, per 1,000, by district, on logarithmic axes. Points are sized by ",
                                    "population; the dashed line marks every infection prevented. ", med,
                                    ". Uttar Pradesh, ", p_up, "; West Bengal, ", p_wb, ".")
)
cap_file <- file.path(FIG_DIR, "captions_counterfactual.txt")
writeLines(paste0(names(caps), ":\n", caps, "\n"), cap_file)
cat("\ncaptions written to", cap_file, "\n")

cat("\ntables in", TAB_DIR, "| figures in", FIG_DIR, "\n")