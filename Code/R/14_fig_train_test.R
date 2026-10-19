# 14_fig_horizons.R
# Observed against predicted across every horizon, one figure per state.
# Predictions cover the test period only. The model was fitted on the earlier
# years, so drawing predictions there would show fit rather than forecast.
# Input  : sfnn_master.csv, predictions_k<K>.csv
# Output : figures/fig_horizons_<state>.pdf and .png

source("/.../.../.../Code/R/00_config.R")

library(tidyverse)
library(patchwork)
library(scales)

FIN <- file.path(OUT_DIR, "sfnn", "final")

KS      <- as.integer(strsplit(Sys.getenv("MS_KS", "1,4,6,12,20,26"), ",")[[1]])
VARIANT <- Sys.getenv("MS_VARIANT", "no_nfhs")

PICK <- c("Uttar Pradesh" = Sys.getenv("MS_UP", "Siddharth Nagar"),
          "West Bengal"   = Sys.getenv("MS_WB", "Murshidabad"))

cat("horizons:", paste(KS, collapse = ", "), "| variant", VARIANT, "\n")

OBS_COL  <- "#1A1A2E"
PRED_COL <- "#C0392B"

theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.background  = element_rect(fill = "#FAFAFA", colour = NA),
      panel.border      = element_rect(colour = "#CCCCCC", fill = NA,
                                       linewidth = 0.5),
      panel.grid.major  = element_line(colour = "#E5E5E5", linewidth = 0.35),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(size = rel(0.82), colour = "#444444"),
      axis.title.y      = element_text(size = rel(0.82), colour = "#444444",
                                       angle = 90, margin = margin(r = 5)),
      axis.text         = element_text(size = rel(0.74), colour = "#666666"),
      axis.ticks        = element_line(colour = "#BBBBBB", linewidth = 0.35),
      axis.ticks.length = unit(2.5, "pt"),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.text       = element_text(size = rel(0.85), colour = "#444444"),
      legend.key        = element_rect(fill = NA, colour = NA),
      legend.key.width  = unit(1.6, "cm"),
      legend.key.height = unit(0.45, "cm"),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.margin     = margin(t = 2, b = 0),
      plot.title        = element_text(size = rel(0.95), face = "bold",
                                       colour = "#222222", hjust = 0,
                                       margin = margin(b = 5)),
      plot.title.position = "plot",
      plot.margin       = margin(6, 8, 6, 6)
    )
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
panel_for <- function(d, k, show_y) {
  ggplot(d, aes(time)) +
    geom_line(aes(y = observed, colour = "Observed"),
              linewidth = 0.6, lineend = "round") +
    geom_line(aes(y = predicted, colour = "Predicted"),
              linewidth = 0.55, linetype = "longdash", lineend = "round",
              na.rm = TRUE) +
    scale_colour_manual(
      values = c(Observed = OBS_COL, Predicted = PRED_COL),
      guide = guide_legend(override.aes = list(
        linewidth = c(0.75, 0.6), linetype = c("solid", "longdash")))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                       breaks = seq(2008, 2020, 4)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)),
                       labels = comma) +
    labs(title = paste0("k = ", k), x = NULL,
         y = if (show_y) "Cases per biweek" else NULL) +
    theme_publication()
}

# ---- one figure per state ----------------------------------------------------
for (st in STATES) {
  dist <- PICK[[st]]
  u <- make_uid(st, dist)
  
  if (!u %in% m$uid) {
    cat("\n", dist, "is not in the master for", st, ", skipping\n")
    next
  }
  
  obs <- m %>% filter(uid == u) %>% select(t_global, time, observed = cases)
  
  # every horizon carries the full observed series with its own predictions,
  # so the panels share an x range
  d <- bind_rows(lapply(KS, function(k) {
    obs %>% left_join(pred %>% filter(uid == u, k == !!k) %>%
                        select(t_global, predicted),
                      by = "t_global") %>%
      mutate(k = k)
  }))
  
  cat("\n===", st, "|", dist, "===\n")
  print(as.data.frame(d %>% filter(!is.na(predicted)) %>%
                        group_by(k) %>%
                        summarise(biweeks = n(),
                                  observed_total = sum(observed),
                                  predicted_total = round(sum(predicted)),
                                  rmse = round(sqrt(mean((predicted - observed)^2)), 1),
                                  correlation = round(cor(predicted, observed), 3),
                                  .groups = "drop")))
  
  # the y label goes on the left column only
  ncol <- 3
  panels <- lapply(seq_along(KS), function(i)
    panel_for(d %>% filter(k == KS[i]), KS[i],
              show_y = ((i - 1) %% ncol == 0)))
  
  fig <- wrap_plots(panels, ncol = ncol) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  out <- file.path(FIG_DIR,
                   paste0("fig_horizons_",
                          gsub("[^a-z0-9]+", "_", tolower(st))))
  ggsave(paste0(out, ".pdf"), fig, width = 13, height = 7, device = cairo_pdf)
  ggsave(paste0(out, ".png"), fig, width = 13, height = 7, dpi = 500,
         bg = "white")
  cat("saved", basename(out), "\n")
}

cat("\nfigures in", FIG_DIR, "\n")
cat("\ncaption: Reported measles cases and model predictions across forecast",
    "horizons, for one district. Predictions cover the test period only, since",
    "the model was fitted on the earlier years. Panels have independent",
    "vertical scales.\n")