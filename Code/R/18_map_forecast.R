# 18_map_forecast.R

source("/.../.../.../Code/R/00_config.R")
source("/.../.../.../Code/R/01_district_names.R")


library(tidyverse)
library(sf)
library(scales)
library(patchwork)

FC_DIR <- file.path(OUT_DIR, "sfnn", "forecast")

K       <- as.integer(Sys.getenv("MS_K", "26"))
VARIANT <- Sys.getenv("MS_VARIANT", "no_nfhs")
SHAPE   <- Sys.getenv("MS_SHAPEFILE",
                      file.path(dirname(DATA_DIR), "shapefiles",
                                "gadm36_IND_2_sp.rds"))
LABEL_ON <- tolower(Sys.getenv("MS_MAP_LABELS", "on")) %in% c("on", "true", "1")

# ---- print size --------------------------------------------------------------
# Saved at the size it will occupy on the page, so LaTeX places it 1:1 and the
# text and boundary lines are not shrunk. 
FIG_W   <- as.numeric(Sys.getenv("MS_FIG_W", "5.5"))
FIG_H   <- as.numeric(Sys.getenv("MS_FIG_H", "6.5"))
FIG_DPI <- as.integer(Sys.getenv("MS_FIG_DPI", "600"))
BORDER  <- as.numeric(Sys.getenv("MS_BORDER", "0.25"))   # boundary linewidth

cat("k =", K, "| variant", VARIANT, "| district labels",
    if (LABEL_ON) "on" else "off",
    "| figure", FIG_W, "x", FIG_H, "in at", FIG_DPI, "dpi\n")
if (!file.exists(SHAPE)) stop("no shapefile at ", SHAPE)

refs <- load_all_references(STATES, file_of)

# ---- the susceptible fraction ------------------------------------------------

m <- read_csv(file.path(OUT_DIR, "sfnn_master.csv"), show_col_types = FALSE)
last_t <- max(m$t_global)
from_t <- last_t - SUSC_AVG_BIWEEKS + 1L

susc <- m %>%
  filter(t_global >= from_t) %>%
  group_by(state, district, uid) %>%
  summarise(susc_pct   = 100 * mean(susc / population),
            population = mean(population), .groups = "drop")

cat("susceptible fraction averaged over t_global", from_t, "to", last_t, "\n")
print(as.data.frame(susc %>% group_by(state) %>%
                      summarise(districts = n(), min = round(min(susc_pct), 1),
                                median = round(median(susc_pct), 1),
                                max = round(max(susc_pct), 1), .groups = "drop")))

# ---- the forecast ------------------------------------------------------------
f <- file.path(FC_DIR, paste0("forecast_summary_k", K, "_", VARIANT, ".csv"))
HAVE_FC <- file.exists(f)

if (HAVE_FC) {
  fc <- read_csv(f, show_col_types = FALSE) %>%
    left_join(susc %>% select(uid, population), by = "uid") %>%
    mutate(rate = 1e5 * predicted_total / population)

  cat("\nforecast incidence per 100,000 over the year:\n")
  print(as.data.frame(fc %>% group_by(state) %>%
                        summarise(districts = n(), min = round(min(rate), 1),
                                  median = round(median(rate), 1),
                                  max = round(max(rate), 1), .groups = "drop")))
} else {
  cat("\nno forecast file; the map is susceptibility alone\n")
  fc <- tibble(uid = character(), rate = numeric(), predicted_total = numeric())
}

# ---- boundaries --------------------------------------------------------------
g <- readRDS(SHAPE)
shp <- st_as_sf(g) %>% filter(NAME_1 %in% STATES) %>%
  group_by(NAME_1) %>%
  group_modify(~ {
    ref <- refs[[.y$NAME_1]]
    .x %>% mutate(district = standardise_districts(
      NAME_2, ref, source = paste(.y$NAME_1, "shapefile")))
  }) %>%
  ungroup() %>% st_as_sf() %>%
  # Where standardisation fails the shapefile keeps its own name, so a polygon
  # whose district was later split still has a key the pooling can match.
  # Without this, Barddhaman resolves to NA and the pooling below never fires.
  mutate(state = NAME_1,
         district = if_else(is.na(district), trimws(NAME_2), district),
         uid = make_uid(state, district)) %>%
  select(state, district, uid, geometry)

# ---- districts the shapefile draws as one ------------------------------------
# A district created after the study began has no estimate of its own. Rather
# than leave a hole, it takes its parent's value, which is the estimate that
# covered its territory while the cases were reported.
INHERIT <- c(`West Bengal__Alipurduar` = "West Bengal__Jalpaiguri")

POOL <- list(
  `West Bengal__Barddhaman` = c("West Bengal__Purba Barddhaman",
                                "West Bengal__Paschim Barddhaman"))

for (poly in names(POOL)) {
  have <- intersect(POOL[[poly]], susc$uid)
  if (!poly %in% shp$uid || !length(have)) next

  bits    <- strsplit(poly, UID_SEP, fixed = TRUE)[[1]]
  st_name <- bits[1]; d_name <- bits[2]

  ps <- susc %>% filter(uid %in% have)
  susc <- bind_rows(
    susc %>% filter(!uid %in% have),
    tibble(state = st_name, district = d_name, uid = poly,
           susc_pct = weighted.mean(ps$susc_pct, ps$population),
           population = sum(ps$population)))

  pf <- fc %>% filter(uid %in% have)
  if (nrow(pf)) {
    fc <- bind_rows(
      fc %>% filter(!uid %in% have),
      tibble(state = st_name, district = d_name, uid = poly,
             biweeks = max(pf$biweeks),
             predicted_total = sum(pf$predicted_total),
             predicted_mean = sum(pf$predicted_mean),
             population = sum(ps$population)) %>%
        mutate(rate = 1e5 * predicted_total / population))
  }

  cat("\npooled onto one polygon: ",
      paste(sub(paste0("^.*", UID_SEP), "", have), collapse = " + "),
      " -> ", d_name, "\n", sep = "")
}

for (child in names(INHERIT)) {
  parent <- INHERIT[[child]]
  if (!child %in% shp$uid || child %in% susc$uid) next
  if (!parent %in% susc$uid) next
  bits <- strsplit(child, UID_SEP, fixed = TRUE)[[1]]

  susc <- bind_rows(susc,
                    susc %>% filter(uid == parent) %>%
                      mutate(state = bits[1], district = bits[2], uid = child))

  if (HAVE_FC && parent %in% fc$uid)
    fc <- bind_rows(fc,
                    fc %>% filter(uid == parent) %>%
                      mutate(state = bits[1], district = bits[2], uid = child))

  cat("\n", bits[2], " takes the estimate of ",
      sub(paste0("^.*", UID_SEP), "", parent), "\n", sep = "")
}

miss <- setdiff(susc$uid, shp$uid)
if (length(miss)) {
  cat("\ndistricts with data but no boundary, not drawn:", length(miss), "\n")
  for (u in miss) cat("  ", u, "\n")
}
blank <- setdiff(shp$uid, susc$uid)
if (length(blank)) {
  cat("boundaries with no data, drawn grey:", length(blank), "\n")
  for (u in blank) cat("  ", u, "\n")
}

mapdat <- shp %>%
  left_join(susc %>% select(uid, susc_pct), by = "uid") %>%
  left_join(fc %>% select(uid, rate, predicted_total), by = "uid") %>%
  mutate(band_susc = cut(susc_pct, SUSC_BREAKS, labels = SUSC_LABELS,
                         right = FALSE),
         band_rate = if (HAVE_FC)
           cut(rate, RATE_BREAKS, labels = RATE_LABELS, right = FALSE)
         else factor(NA, levels = RATE_LABELS))

cat("\ndistricts in each susceptibility band:\n")
print(as.data.frame(mapdat %>% st_drop_geometry() %>%
                      count(state, band_susc, .drop = FALSE) %>%
                      pivot_wider(names_from = band_susc, values_from = n)))
if (HAVE_FC) {
  cat("\ndistricts in each forecast band:\n")
  print(as.data.frame(mapdat %>% st_drop_geometry() %>%
                        count(state, band_rate, .drop = FALSE) %>%
                        pivot_wider(names_from = band_rate, values_from = n)))
}

# ---- one panel ---------------------------------------------------------------
# Four categories, from SUSC_BREAKS and SUSC_LABELS in the config, the same in
# both panels. drop = FALSE keeps all four in the legend even where a category
# occurs in only one state, so the one collected legend is complete.

SUSC_MAP_PAL <- c("#009E73", "#F0E442", "#E69F00", "#D55E00")
SUSC_MAP_TXT <- c("white", "black", "black", "white")

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


panel <- function(st, tag, title) {
  d <- mapdat %>% filter(state == st) %>%
    mutate(txt = coalesce(SUSC_MAP_TXT[as.integer(band_susc)], "black"))

  p <- ggplot(d) +
    geom_sf(aes(fill = band_susc), colour = "white", linewidth = BORDER)

  # every district named; check_overlap drops any that would collide, and the
  # full lists of the highest districts are printed below
  if (LABEL_ON)
    p <- p + geom_sf_text(aes(label = district, colour = txt),
                          size = MAP_LABEL_SIZE, fontface = "bold",
                          check_overlap = TRUE)

  p +
    scale_fill_manual(values = setNames(SUSC_MAP_PAL, SUSC_LABELS), drop = FALSE,
                      na.value = "grey88",
                      name = "Susceptible (S/N),\n% of population",
                      guide = guide_legend(ncol = 1, reverse = TRUE,
                                           title.position = "top", title.hjust = 0,
                                           keyheight = unit(16, "pt"),
                                           keywidth = unit(14, "pt"),
                                           byrow = TRUE,
                                           override.aes = list(colour = NA))) +
    scale_colour_identity() +
    labs(tag = tag, title = title) +
    coord_sf(xlim = lim_of(st)$x, ylim = lim_of(st)$y, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.tag          = element_text(face = "bold", size = 10),
      plot.tag.position = c(0.02, 0.98),
      legend.position   = "right",
      legend.title      = element_text(size = 8, face = "bold"),
      legend.text       = element_text(size = 8),
      plot.margin       = margin(2, 2, 2, 2)
    )
}

pa <- panel(STATES[1], "a", STATES[1])
pb <- panel(STATES[2], "b", STATES[2])


p <- (pa / pb) +
  plot_layout(nrow = 2, ncol = 1, heights = panel_heights, guides = "collect") &
  theme(legend.position = "right") & key_gap

# ---- save --------------------------------------------------------------------
# PDF is vector: use it in LaTeX. PNG for slides and previews, and TIFF with
# LZW compression for the journal upload, both at FIG_DPI and written with ragg
# when available, which gives cleaner anti-aliasing.
out <- file.path(FIG_DIR, "fig_susceptibility_map")

ggsave(paste0(out, ".pdf"), p, width = FIG_W, height = FIG_H,
       device = cairo_pdf)

has_ragg <- requireNamespace("ragg", quietly = TRUE)
ggsave(paste0(out, ".png"), p, width = FIG_W, height = FIG_H, dpi = FIG_DPI,
       bg = "white", device = if (has_ragg) ragg::agg_png else "png")
if (has_ragg) {
  ggsave(paste0(out, ".tif"), p, width = FIG_W, height = FIG_H, dpi = FIG_DPI,
         bg = "white", device = ragg::agg_tiff, compression = "lzw")
} else {
  ggsave(paste0(out, ".tif"), p, width = FIG_W, height = FIG_H, dpi = FIG_DPI,
         bg = "white", device = "tiff", compression = "lzw")
}
cat("\nsaved", basename(out), "as pdf, png and tif\n")

writeLines(paste0("fig_susceptibility_map:\n",
  "Reconstructed susceptible fraction by district at the end of the series. ",
  "(a) Uttar Pradesh; (b) West Bengal. Percentage of population, averaged over the ",
  "final ", SUSC_AVG_BIWEEKS, " biweeks, in four categories shared by both panels; both ",
  "states are drawn at the same scale. ",
  "Barddhaman appears undivided because the boundary data predate its 2017 ",
  "division; Alipurduar, created from Jalpaiguri in 2014, is shown with Jalpaiguri's ",
  "estimate."), file.path(FIG_DIR, "caption_susceptibility_map.txt"))

# ---- the districts in the top categories -------------------------------------
# Printed in full, since a label that collided was dropped from the map.
top_of <- function(band_col, val_col) {
  lv <- levels(mapdat[[band_col]])
  mapdat %>% st_drop_geometry() %>%
    filter(!is.na(.data[[band_col]]),
           as.character(.data[[band_col]]) == tail(lv, 1)) %>%
    arrange(desc(.data[[val_col]]))
}

for (s in STATES) {
  cat("\n===", s, "===\n")
  hs <- top_of("band_susc", "susc_pct") %>% filter(state == s)
  cat("  in the highest band (", nrow(hs), "): ",
      paste(hs$district, collapse = ", "), "\n", sep = "")
  if (HAVE_FC) {
    hr <- top_of("band_rate", "rate") %>% filter(state == s)
    cat("  highest forecast band (", nrow(hr), "): ",
        paste(hr$district, collapse = ", "), "\n", sep = "")
    cat("  in both: ",
        paste(intersect(hs$district, hr$district), collapse = ", "),
        "\n", sep = "")
  }
}

# Whether the two panels agree is worth knowing. If susceptibility and forecast
# rank districts differently, the forecast carries information the
# reconstruction does not.

if (HAVE_FC) {
  cmp <- mapdat %>% st_drop_geometry() %>%
    filter(is.finite(susc_pct), is.finite(rate))
  if (nrow(cmp) > 2)
    cat("\nrank correlation, susceptibility against forecast:",
        round(cor(cmp$susc_pct, cmp$rate, method = "spearman"), 3), "\n")
}

write_csv(mapdat %>% st_drop_geometry(),
          file.path(FC_DIR, "susceptibility_map_data.csv"))
