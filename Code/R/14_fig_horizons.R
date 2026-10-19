# 14_fig_horizons.R
# Observed against predicted across every horizon, for one district per state.
# Predictions cover the test period only. 
# Input  : sfnn_master.csv, predictions_k<K>.csv
# Output : figures/fig_horizons.{pdf,png,tif} 

source("/.../.../.../Code/R/00_config.R")

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(scales)
})


FIN <- file.path(OUT_DIR, "sfnn", "final")

KS      <- as.integer(strsplit(Sys.getenv("MS_KS", "1,4,6,12"), ",")[[1]])
VARIANT <- Sys.getenv("MS_VARIANT", "no_nfhs")
FIG_DPI <- as.integer(Sys.getenv("MS_FIG_DPI", "600"))

# from the config, so the same districts appear in every figure that names them
PICK <- FIG_DISTRICTS
if (nzchar(Sys.getenv("MS_UP"))) PICK[["Uttar Pradesh"]] <- Sys.getenv("MS_UP")
if (nzchar(Sys.getenv("MS_WB"))) PICK[["West Bengal"]]   <- Sys.getenv("MS_WB")

cat("horizons:", paste(KS, collapse = ", "), "| variant", VARIANT, "\n")

OBS_COL  <- "black"
PRED_COL <- "#D55E00"     # Okabe-Ito vermillion
VAL_FILL <- "#F4A582"     # validation period, as before
TST_FILL <- "#92C5DE"     # test period, as before

# The same cuts the loader uses, in decimal years so they sit on the time axis.
# A fiscal year starts at biweek 7, hence the offset.
VAL_FROM  <- BASE_YEAR + ((2016 - BASE_YEAR) * BIWEEKS_PER_YEAR + 6) / BIWEEKS_PER_YEAR
TEST_FROM <- BASE_YEAR + ((2017 - BASE_YEAR) * BIWEEKS_PER_YEAR + 6) / BIWEEKS_PER_YEAR

theme_clean <- function(base_size = 10) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      panel.grid        = element_blank(),
      axis.line         = element_blank(),
      axis.title        = element_text(size = rel(0.9), colour = "black"),
      axis.title.y      = element_text(size = rel(0.9), colour = "black", angle = 90,
                                       margin = margin(r = 5)),
      axis.text         = element_text(size = rel(0.82), colour = "black"),
      axis.ticks        = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(2.5, "pt"),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.text       = element_text(size = rel(0.9), colour = "black"),
      legend.key        = element_rect(fill = NA, colour = NA),
      legend.key.width  = unit(1.6, "cm"),
      legend.background = element_rect(fill = NA, colour = NA),
      plot.title        = element_text(size = rel(0.95), face = "bold", colour = "black",
                                       hjust = 0, margin = margin(b = 4)),
      plot.title.position = "plot",
      plot.margin       = margin(5, 7, 5, 5)
    )
}

has_ragg <- requireNamespace("ragg", quietly = TRUE)
save_all <- function(fig, stem, w, h) {
  ggsave(paste0(stem, ".pdf"), fig, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(stem, ".png"), fig, width = w, height = h, dpi = FIG_DPI, bg = "white",
         device = if (has_ragg) ragg::agg_png else "png")
  if (has_ragg) {
    ggsave(paste0(stem, ".tif"), fig, width = w, height = h, dpi = FIG_DPI, bg = "white",
           device = ragg::agg_tiff, compression = "lzw")
  } else {
    ggsave(paste0(stem, ".tif"), fig, width = w, height = h, dpi = FIG_DPI, bg = "white",
           device = "tiff", compression = "lzw")
  }
  cat("saved", basename(stem), "as pdf, png and tif\n")
}

# ---- the data ----------------------------------------------------------------
m <- read_csv(file.path(OUT_DIR, "sfnn_master.csv"), show_col_types = FALSE)

load_k <- function(k) {
  f <- file.path(FIN, paste0("predictions_k", k, ".csv"))
  if (!file.exists(f)) { cat("  no predictions for k =", k, "\n"); return(NULL) }
  read_csv(f, show_col_types = FALSE) %>%
    filter(variant == VARIANT) %>%
    select(uid, t_global, predicted) %>%
    mutate(k = k)
}

pred <- bind_rows(lapply(KS, load_k))
if (!nrow(pred)) stop("no predictions found for variant ", VARIANT)
KS <- sort(unique(pred$k))
cat("horizons available:", paste(KS, collapse = ", "), "\n")

# ---- one panel ---------------------------------------------------------------
panel_for <- function(d, k, show_y, show_period_labels = FALSE) {
  te <- d %>% filter(!is.na(predicted))
  r <- if (nrow(te) > 2 && sd(te$predicted) > 0 && sd(te$observed) > 0)
    cor(te$predicted, te$observed) else NA_real_
  rm_ <- if (nrow(te)) sqrt(mean((te$predicted - te$observed)^2)) else NA_real_
  ymax <- max(d$observed, d$predicted, na.rm = TRUE)

  ggplot(d, aes(time)) +
    annotate("rect", xmin = VAL_FROM, xmax = TEST_FROM, ymin = -Inf, ymax = Inf,
             fill = VAL_FILL, alpha = 0.16) +
    annotate("rect", xmin = TEST_FROM, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = TST_FILL, alpha = 0.20) +
    geom_vline(xintercept = c(VAL_FROM, TEST_FROM),
               colour = "#999999", linetype = 2, linewidth = 0.3) +
    geom_line(aes(y = observed, colour = "Observed"), linewidth = 0.6, lineend = "round") +
    geom_line(aes(y = predicted, colour = "Predicted"), linewidth = 0.65,
              linetype = "longdash", lineend = "round", na.rm = TRUE) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.2f\nRMSE = %.1f", r, rm_),
             hjust = 1.08, vjust = 1.25, size = 2.7, colour = "black", lineheight = 0.95) +
    scale_colour_manual(
      values = c(Observed = OBS_COL, Predicted = PRED_COL),
      guide = guide_legend(override.aes = list(
        linewidth = c(0.75, 0.75), linetype = c("solid", "longdash")))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                       breaks = seq(2008, 2020, 4)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)), labels = comma) +
    {if (show_period_labels)
      annotate("text",
               x = c(mean(c(min(d$time), VAL_FROM)), mean(c(VAL_FROM, TEST_FROM)),
                     mean(c(TEST_FROM, max(d$time)))),
               y = ymax * 1.03, label = c("train", "validation", "test"),
               size = 2.4, colour = "black", vjust = 0)
     else NULL} +
    labs(title = paste0("k = ", k), x = NULL,
         y = if (show_y) "Cases per biweek" else NULL) +
    theme_clean()
}

# ---- one row of panels for one state ------------------------------------------
NCOL <- if (length(KS) <= 4) length(KS) else ceiling(length(KS) / 2)
NROW <- ceiling(length(KS) / NCOL)

state_row <- function(st) {
  dist <- PICK[[st]]
  u <- make_uid(st, dist)
  if (!u %in% m$uid) { cat("\n", dist, "is not in the master for", st, ", skipping\n"); return(NULL) }

  obs <- m %>% filter(uid == u) %>% select(t_global, time, observed = cases)
  d <- bind_rows(lapply(KS, function(k) {
    obs %>% left_join(pred %>% filter(uid == u, k == !!k) %>% select(t_global, predicted),
                      by = "t_global") %>% mutate(k = k)
  }))

  cat("\n===", st, "|", dist, "===\n")
  print(as.data.frame(d %>% filter(!is.na(predicted)) %>% group_by(k) %>%
    summarise(biweeks = n(), observed_total = sum(observed),
              predicted_total = round(sum(predicted)),
              rmse = round(sqrt(mean((predicted - observed)^2)), 1),
              correlation = round(cor(predicted, observed), 3), .groups = "drop")))

  panels <- lapply(seq_along(KS), function(i)
    panel_for(d %>% filter(k == KS[i]), KS[i],
              show_y = ((i - 1) %% NCOL == 0), show_period_labels = (i == 1)))
  list(dist = dist, row = wrap_plots(panels, ncol = NCOL))
}

rows <- setNames(lapply(STATES, state_row), STATES)
rows <- rows[!sapply(rows, is.null)]

W <- 7.5                 # PLOS: width <= 7.5 in
H <- 4.1                 # two states stacked: 2H + 0.4 <= 8.75 in

# ---- one figure per state ------------------------------------------------------
for (st in names(rows)) {
  fig <- rows[[st]]$row + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save_all(fig, file.path(FIG_DIR, paste0("fig_horizons_", gsub("[^a-z0-9]+", "_", tolower(st)))),
           W, H)
}

# ---- both states in one figure, A above B, one legend ----------------------------
if (length(rows) == 2) {
  header <- function(txt) wrap_elements(full = grid::textGrob(
    txt, x = 0.005, hjust = 0, gp = grid::gpar(fontsize = 10.5, fontface = "bold")))
  st <- names(rows)
  fig <- header(paste0("A   ", rows[[1]]$dist, ", ", st[1])) /
    (rows[[1]]$row & theme(legend.position = "none")) /
    header(paste0("B   ", rows[[2]]$dist, ", ", st[2])) /
    (rows[[2]]$row + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
    plot_layout(heights = c(0.08, NROW, 0.08, NROW))
  save_all(fig, file.path(FIG_DIR, "fig_horizons"), W, 2 * H + 0.4)
}

cat("\nfigures in", FIG_DIR, "\n")


