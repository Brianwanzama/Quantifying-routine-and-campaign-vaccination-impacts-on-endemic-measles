# 25_compare_sfnn_tsir.R
# SFNN against TSIR, district by district.
# Input  : final/predictions_k<K>.csv, tsir_forecast/tsir_predictions_w<W>.csv,
#          sfnn_master.csv
# Output : figures/compare/compare_rmse_*.pdf and .png

source("/.../.../.../Code/R/00_config.R")

library(tidyverse)
library(patchwork)
library(scales)

set.seed(42)

FIN     <- file.path(OUT_DIR, "sfnn", "final")
TSIR_F  <- file.path(OUT_DIR, "tsir_forecast")
OUT_FIG <- file.path(FIG_DIR, "compare")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# the horizons both methods ran; TSIR did 1, 4, 6, 12
KS      <- as.integer(strsplit(Sys.getenv("MS_KS", "1,4,6,12"), ",")[[1]])
VARIANT <- Sys.getenv("MS_VARIANT", "no_nfhs")
WINDOW  <- as.integer(Sys.getenv("MS_TSIR_WINDOW", "130"))

cat("horizons:", paste(KS, collapse = ", "),
    "| SFNN variant", VARIANT, "| TSIR window", WINDOW, "\n")

# ---- population, per district -------------------------------------------------
master <- read_csv(file.path(OUT_DIR, "sfnn_master.csv"), show_col_types = FALSE)
pop <- master %>%
  group_by(state, district, uid) %>%
  summarise(min_pop = min(population, na.rm = TRUE), .groups = "drop")

# ---- the SFNN, bound across horizons ------------------------------------------
read_sfnn <- function(k) {
  f <- file.path(FIN, paste0("predictions_k", k, ".csv"))
  if (!file.exists(f)) { cat("  no SFNN predictions for k =", k, "\n"); return(NULL) }
  read_csv(f, show_col_types = FALSE) %>%
    filter(variant == VARIANT) %>%
    transmute(uid, t_global, k = k, sfnn = predicted, observed)
}

sfnn <- bind_rows(lapply(KS, read_sfnn))
if (!nrow(sfnn)) stop("no SFNN predictions for variant ", VARIANT)
cat("SFNN rows:", nrow(sfnn), "| horizons:",
    paste(sort(unique(sfnn$k)), collapse = ", "), "\n")

# ---- TSIR ----------------------------------------------------------------------
tf <- file.path(TSIR_F, paste0("tsir_predictions_w", WINDOW, ".csv"))
if (!file.exists(tf)) stop("no TSIR predictions at ", tf)

tsir <- read_csv(tf, show_col_types = FALSE) %>%
  filter(k %in% KS, period == "test") %>%
  transmute(uid, t_global, k, tsir = pmax(predicted, 0))

cat("TSIR rows over the test period:", nrow(tsir), "\n")

# ---- join, and say what it cost ------------------------------------------------
dat <- sfnn %>%
  inner_join(tsir, by = c("uid", "t_global", "k")) %>%
  left_join(pop, by = "uid")

cat("\nrows after the join:", nrow(dat), "\n")
print(as.data.frame(
  bind_rows(
    sfnn %>% count(k, name = "sfnn_only"),
    tsir %>% count(k, name = "tsir_only"),
    dat  %>% count(k, name = "both")) %>%
    group_by(k) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")))

lost_d <- setdiff(unique(sfnn$uid), unique(dat$uid))
if (length(lost_d)) {
  cat("\ndistricts the SFNN has but the join drops:", length(lost_d), "\n")
  for (u in head(lost_d, 8)) cat("  ", u, "\n")
  cat("usually where TSIR failed to fit at most origins\n")
}

# ---- per district RMSE, on identical rows ---------------------------------------
rmse <- dat %>%
  filter(is.finite(sfnn), is.finite(tsir), is.finite(observed)) %>%
  group_by(state, district, uid, k) %>%
  summarise(n = n(),
            min_pop   = first(min_pop),
            tsir_rmse = sqrt(mean((observed - tsir)^2)),
            nn_rmse   = sqrt(mean((observed - sfnn)^2)),
            .groups = "drop") %>%
  mutate(k_label = factor(paste0("k = ", k), levels = paste0("k = ", KS)))

cat("\nmean RMSE by horizon:\n")
print(as.data.frame(rmse %>% group_by(k) %>%
  summarise(districts = n(),
            mean_tsir = round(mean(tsir_rmse), 2),
            mean_sfnn = round(mean(nn_rmse), 2),
            improvement_pct = round(100 * (mean(tsir_rmse) - mean(nn_rmse)) /
                                      mean(tsir_rmse), 1),
            sfnn_wins_pct = round(100 * mean(nn_rmse < tsir_rmse), 1),
            .groups = "drop")))

# The mean is carried by a few districts where TSIR ran away, so the median and
# a paired test say more about the typical district.
cat("\nmedian RMSE, and a paired test across districts:\n")
print(as.data.frame(rmse %>% group_by(k) %>%
  summarise(median_tsir = round(median(tsir_rmse), 2),
            median_sfnn = round(median(nn_rmse), 2),
            median_gain = round(median(tsir_rmse - nn_rmse), 2),
            p = signif(wilcox.test(tsir_rmse, nn_rmse, paired = TRUE)$p.value, 3),
            .groups = "drop")))

cat("\nby state:\n")
print(as.data.frame(rmse %>% group_by(state, k) %>%
  summarise(districts = n(),
            median_tsir = round(median(tsir_rmse), 2),
            median_sfnn = round(median(nn_rmse), 2),
            sfnn_wins_pct = round(100 * mean(nn_rmse < tsir_rmse), 1),
            .groups = "drop")))

# ---- axis limits ----------------------------------------------------------------
# The ninetieth percentile, with anything beyond squished to the edge. 
rmse_limit <- ceiling(quantile(c(rmse$tsir_rmse, rmse$nn_rmse), 0.90, na.rm = TRUE))
gain_limit <- ceiling(quantile(abs(rmse$tsir_rmse - rmse$nn_rmse), 0.90, na.rm = TRUE))
cat(sprintf("\naxis limits: rmse %.0f | gain %.0f\n", rmse_limit, gain_limit))

n_out <- sum(rmse$tsir_rmse > rmse_limit | rmse$nn_rmse > rmse_limit)
cat("points squished to the edge:", n_out, "of", nrow(rmse), "\n")

# ---- theme -----------------------------------------------------------------------
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(axis.line        = element_line(colour = "#1a1a1a", linewidth = 0.3),
          panel.border     = element_rect(colour = "#1a1a1a", fill = NA,
                                          linewidth = 0.4),
          panel.grid.major = element_line(colour = "#e8e8e8", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#f0f0f0", colour = "#1a1a1a",
                                          linewidth = 0.4),
          strip.text       = element_text(face = "bold", size = base_size * 0.85),
          legend.position  = "bottom",
          legend.key.width = unit(1.8, "cm"),
          legend.key.height = unit(0.3, "cm"),
          plot.title = element_blank(), plot.subtitle = element_blank(),
          plot.caption = element_blank(),
          plot.tag = element_text(size = base_size * 1.2, face = "bold"))
}

PAL <- c("#002A66", "#1a6aa8", "#f0c400", "#FBE045")

# ---- A: one method against the other ---------------------------------------------
# The diagonal is equality: below it the SFNN has the lower error.
pa <- ggplot(rmse, aes(tsir_rmse, nn_rmse, colour = log(min_pop))) +
  geom_abline(colour = "#888888", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.2, alpha = 0.85, shape = 16) +
  scale_color_gradientn(colours = PAL, name = "Log(Population)",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = 10, barheight = 0.5,
                           ticks.linewidth = 0.5)) +
  scale_x_continuous(limits = c(0, rmse_limit), oob = scales::squish) +
  scale_y_continuous(limits = c(0, rmse_limit), oob = scales::squish) +
  facet_wrap(~ k_label, ncol = length(KS)) +
  labs(x = expression(RMSE[TSIR]), y = expression(RMSE[SFNN])) +
  theme_pub()

ggsave(file.path(OUT_FIG, "compare_rmse_scatter.png"), pa,
       width = 9, height = 4, dpi = 400, bg = "white")
ggsave(file.path(OUT_FIG, "compare_rmse_scatter.pdf"), pa,
       width = 9, height = 4, device = cairo_pdf)

# ---- B: does the gain depend on district size? -----------------------------------
gain <- rmse %>% mutate(rmse_gain = tsir_rmse - nn_rmse)

pb <- ggplot(gain, aes(log(min_pop), rmse_gain, colour = log(min_pop))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#888888",
             linewidth = 0.4) +
  geom_point(size = 1.2, alpha = 0.85, shape = 16) +
  geom_smooth(method = "loess", colour = "#002A66", fill = "#002A66",
              alpha = 0.15, linewidth = 0.6, se = TRUE) +
  scale_color_gradientn(colours = PAL, name = "Log(Population)",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = 10, barheight = 0.5,
                           ticks.linewidth = 0.5)) +
  scale_y_continuous(limits = c(-gain_limit, gain_limit),
                     oob = scales::squish) +
  facet_wrap(~ k_label, ncol = length(KS)) +
  labs(x = "Log(Population)",
       y = expression(RMSE[TSIR] - RMSE[SFNN])) +
  theme_pub()

ggsave(file.path(OUT_FIG, "compare_rmse_gain.png"), pb,
       width = 9, height = 4, dpi = 400, bg = "white")
ggsave(file.path(OUT_FIG, "compare_rmse_gain.pdf"), pb,
       width = 9, height = 4, device = cairo_pdf)

# ---- combined, one legend ---------------------------------------------------------
both <- pa / guide_area() / pb +
  plot_layout(heights = c(3, 0.4, 3), guides = "collect") +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 11,
                                                        face = "bold")))

ggsave(file.path(OUT_FIG, "compare_rmse_combined.png"), both,
       width = 9, height = 8.5, dpi = 400, bg = "white")
ggsave(file.path(OUT_FIG, "compare_rmse_combined.pdf"), both,
       width = 9, height = 8.5, device = cairo_pdf)

# ---- does the gain track population? -----------------------------------------------
cat("\nrank correlation of the gain with population:\n")
print(as.data.frame(gain %>% group_by(k) %>%
  summarise(spearman = round(cor(rmse_gain, min_pop, method = "spearman",
                                 use = "complete.obs"), 3), .groups = "drop")))
cat("positive means the SFNN's advantage grows with district size\n")

write_csv(rmse, file.path(OUT_FIG, "compare_rmse_by_district.csv"))

cat("\nfigures in", OUT_FIG, "\n")
