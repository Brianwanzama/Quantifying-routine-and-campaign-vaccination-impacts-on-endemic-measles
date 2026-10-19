# 19_nfhs_comparison.R
# Does adding the static covariates improve forecasting?
# Per district paired comparison of the two models, on RMSE in cases over the
# test period, at each horizon.
# The Wilcoxon signed rank test is used rather than a t test because the
# per district RMSE differences are heavily skewed: a few large districts
# dominate the mean while most sit near zero.
#
# Input  : district_errors_k<K>.csv
# Output : tables/nfhs_comparison.csv, figures/fig_nfhs_*.pdf and .png

source("/.../.../.../Code/R/00_config.R")

library(tidyverse)
library(scales)

FIN <- file.path(OUT_DIR, "sfnn", "final")
KS  <- as.integer(strsplit(Sys.getenv("MS_KS", "1,4,6,12"), ",")[[1]])

read_k <- function(k) {
  f <- file.path(FIN, paste0("district_errors_k", k, ".csv"))
  if (!file.exists(f)) { cat("  no errors for k =", k, "\n"); return(NULL) }
  read_csv(f, show_col_types = FALSE) %>% mutate(k = k)
}

e <- bind_rows(lapply(KS, read_k))
if (!nrow(e)) stop("no district errors found")
KS <- sort(unique(e$k))

# ---- the pairs ---------------------------------------------------------------
# gain is positive where the covariates lower the error, so a positive median
# means they help.
pairs <- e %>%
  select(state, district, k, variant, rmse) %>%
  pivot_wider(names_from = variant, values_from = rmse) %>%
  filter(is.finite(with_nfhs), is.finite(no_nfhs)) %>%
  mutate(gain = no_nfhs - with_nfhs,
         helps = gain > 0)

cat("pairs per horizon:", nrow(pairs) / length(KS), "\n")

# ---- the test ----------------------------------------------------------------
res <- pairs %>%
  group_by(k) %>%
  summarise(
    n            = n(),
    median_with  = median(with_nfhs),
    median_without = median(no_nfhs),
    median_gain  = median(gain),
    win_rate     = 100 * mean(helps),
    p_value      = wilcox.test(with_nfhs, no_nfhs, paired = TRUE)$p.value,
    .groups = "drop") %>%
  mutate(verdict = case_when(
    p_value >= 0.05                   ~ "no difference",
    median_gain > 0                   ~ "covariates help",
    TRUE                              ~ "covariates hurt"))

cat("\n---- paired by district, at each horizon ----\n")
print(as.data.frame(res %>%
                      transmute(k, n,
                                median_rmse_with = round(median_with, 2),
                                median_rmse_without = round(median_without, 2),
                                median_gain = round(median_gain, 3),
                                win_rate = round(win_rate, 1),
                                p = signif(p_value, 3),
                                verdict)))

cat("\nHow to read this. The gain is a district's forecasting error without the\n")
cat("covariates minus its error with them, in cases of RMSE. So a positive\n")
cat("gain means including the covariates lowered the error and helped; a\n")
cat("negative gain means they raised it. The win rate is the share of\n")
cat("districts where they helped, and the p value tests whether the paired\n")
cat("difference across all districts is distinguishable from none.\n")

# by state, since the two may differ
cat("\n---- by state ----\n")
print(as.data.frame(pairs %>% group_by(state, k) %>%
                      summarise(n = n(), median_gain = round(median(gain), 3),
                                win_rate = round(100 * mean(helps), 1),
                                p = signif(wilcox.test(with_nfhs, no_nfhs, paired = TRUE)$p.value, 3),
                                .groups = "drop")))

write_csv(res, file.path(TAB_DIR, "nfhs_comparison.csv"))
write_csv(pairs, file.path(TAB_DIR, "nfhs_comparison_pairs.csv"))

# ---- does the effect depend on district size? --------------------------------
# The hypothesis worth testing: covariates might help where the case history is
# sparse and the network has little to learn from, and add nothing where the
# series is informative on its own.
m <- read_csv(file.path(OUT_DIR, "sfnn_master.csv"), show_col_types = FALSE) %>%
  group_by(state, district) %>%
  summarise(pop = mean(population),
            mean_cases = mean(cases), .groups = "drop")

pg <- pairs %>% left_join(m, by = c("state", "district"))

cat("\n---- does the gain depend on population? ----\n")
print(as.data.frame(pg %>% group_by(k) %>%
                      summarise(spearman_pop = round(cor(gain, pop, method = "spearman"), 3),
                                spearman_cases = round(cor(gain, mean_cases, method = "spearman"), 3),
                                .groups = "drop")))
cat("a positive figure means the covariates help more in larger districts\n")

theme_pub <- theme_classic(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold", size = 8.5),
        axis.text  = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 9, colour = "grey15"),
        legend.position = "bottom", legend.title = element_blank())

# ---- the gain against population ---------------------------------------------
p1 <- ggplot(pg, aes(pop, gain, colour = state)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50",
             linewidth = 0.35) +
  geom_point(alpha = 0.55, size = 1) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.5, alpha = 0.15) +
  # horizons in numeric order (1, 4, 6, 12), not alphabetical (1, 12, 4, 6)
  facet_wrap(~ factor(paste0("k = ", k), levels = paste0("k = ", sort(unique(k)))),
             scales = "free_y", nrow = 2) +
  scale_x_log10(labels = label_number(scale = 1e-6, suffix = "M")) +
  scale_colour_manual(values = setNames(c("#2166AC", "#B2182B"), STATES)) +
  labs(x = "Population, log scale",
       y = "Gain from the covariates (cases of RMSE)\npositive = they help") +
  theme_pub

ggsave(file.path(FIG_DIR, "fig_nfhs_gain_vs_population.pdf"), p1,
       width = 7.5, height = 5)
ggsave(file.path(FIG_DIR, "fig_nfhs_gain_vs_population.png"), p1,
       width = 7.5, height = 5, dpi = 400, bg = "white")
cat("\nsaved fig_nfhs_gain_vs_population\n")

# ---- the distribution of the gain --------------------------------------------
p2 <- ggplot(pairs, aes(gain, factor(k, levels = rev(KS)), fill = state)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50",
             linewidth = 0.35) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.3, alpha = 0.75,
               position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = setNames(c("#2166AC", "#B2182B"), STATES)) +
  labs(x = "Gain from the covariates (cases of RMSE), positive = they help",
       y = "Horizon, biweeks") +
  theme_pub

ggsave(file.path(FIG_DIR, "fig_nfhs_gain_distribution.pdf"), p2,
       width = 7, height = 4.5)
ggsave(file.path(FIG_DIR, "fig_nfhs_gain_distribution.png"), p2,
       width = 7, height = 4.5, dpi = 400, bg = "white")
cat("saved fig_nfhs_gain_distribution\n")

