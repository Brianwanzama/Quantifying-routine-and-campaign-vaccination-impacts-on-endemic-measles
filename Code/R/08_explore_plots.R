# 08_explore_plots.R
# What the two states look like before the SFNN split is chosen.
# Output : figures/fig_explore_*.png and .pdf, and a by year table

source("/.../.../.../Code/R/00_config.R")

library(tidyverse)
library(scales)
library(patchwork)

m <- read_csv(file.path(OUT_DIR, "sfnn_master.csv"), show_col_types = FALSE)
cat("master:", nrow(m), "rows |", n_distinct(m$uid), "districts\n")

theme_pub <- theme_classic(base_size = 10) +
  theme(
    axis.title        = element_text(size = 9.5, colour = "grey15"),
    axis.text         = element_text(size = 8.5, colour = "black"),
    axis.line         = element_line(colour = "grey20", linewidth = 0.5),
    axis.ticks        = element_line(colour = "grey20", linewidth = 0.4),
    strip.background  = element_rect(fill = "grey94", colour = NA),
    strip.text        = element_text(face = "bold", size = 8.5),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8.5),
    plot.margin       = margin(4, 8, 4, 6)
  )
theme_set(theme_pub)

pal_state <- setNames(c("#2166AC", "#B2182B"), STATES)

sv <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0("fig_explore_", name, ".pdf")), p,
         width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0("fig_explore_", name, ".png")), p,
         width = w, height = h, dpi = 400)
  cat("  saved fig_explore_", name, "\n", sep = "")
}

# ---- the table the split rests on --------------------------------------------
yr <- m %>%
  group_by(state, year) %>%
  summarise(districts = n_distinct(uid),
            biweeks   = n_distinct(t_global),
            total     = sum(cases, na.rm = TRUE),
            mean      = mean(cases, na.rm = TRUE),
            sd        = sd(cases, na.rm = TRUE),
            zero      = mean(cases == 0, na.rm = TRUE),
            max       = max(cases, na.rm = TRUE),
            .groups = "drop")

cat("\n---- cases by year ----\n")
print(as.data.frame(yr %>%
                      transmute(state, year, biweeks, total,
                                mean = round(mean, 2), sd = round(sd, 2),
                                zero = round(zero, 3), max)))

write_csv(yr, file.path(TAB_DIR, "explore_cases_by_year.csv"))

# A year is flagged quiet if its mean is far below the state's own median year.
# That is what forced the year list split in the earlier work, so the same test
# is applied here rather than assumed not to matter.
quiet <- yr %>%
  group_by(state) %>%
  mutate(ref = median(mean), ratio = mean / ref) %>%
  ungroup() %>%
  filter(ratio < 0.25)

cat("\n---- years far below their state's median ----\n")
if (nrow(quiet)) {
  print(as.data.frame(quiet %>%
                        transmute(state, year, mean = round(mean, 2),
                                  state_median = round(ref, 2), ratio = round(ratio, 3))))
  cat("these would need excluding, or the split listing years rather than\n")
  cat("cutting at a date\n")
} else {
  cat("none. A date cut will serve, and the same one can serve both states\n")
}

# ---- 1. the state series -----------------------------------------------------
st <- m %>% group_by(state, time) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(st, aes(time, cases, colour = state)) +
  geom_line(linewidth = 0.4) +
  # free rather than free_y, so each panel carries its own x axis as well as
  # its own y and can be read without reference to the other
  facet_wrap(~ state, ncol = 1, scales = "free") +
  scale_colour_manual(values = pal_state, guide = "none") +
  scale_x_continuous(breaks = seq(2008, 2020, 2)) +
  scale_y_continuous(labels = comma) +
  labs(x = "year", y = "reported cases per biweek")
sv(p1, "state_series", 8, 5.5)

# ---- 2. the seasonal shape ---------------------------------------------------
# Whether the two states peak at the same time of year matters for the seasonal
# features, and for whether one model can serve both.
seas <- m %>%
  group_by(state, biweek) %>%
  summarise(mean = mean(cases, na.rm = TRUE), .groups = "drop")

p2 <- ggplot(seas, aes(biweek, mean, colour = state)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = pal_state) +
  scale_x_continuous(breaks = seq(1, 26, 5)) +
  labs(x = "biweek of the year", y = "mean cases per district")
sv(p2, "seasonality", 6, 3.5)

# ---- 3. how sparse, district by district -------------------------------------
# The share of biweeks with no cases is what makes the forecasting hard, and it
# is worth seeing whether the two states differ in it.
sp <- m %>% group_by(state, uid, district) %>%
  summarise(zero = mean(cases == 0, na.rm = TRUE),
            total = sum(cases, na.rm = TRUE),
            pop = mean(population), .groups = "drop")

cat("\n---- share of biweeks with no cases ----\n")
print(as.data.frame(sp %>% group_by(state) %>%
                      summarise(districts = n(),
                                min = round(min(zero), 3), median = round(median(zero), 3),
                                max = round(max(zero), 3), .groups = "drop")))

p3 <- ggplot(sp, aes(zero, fct_rev(state), fill = state)) +
  geom_boxplot(width = 0.5, linewidth = 0.35, outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = pal_state, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "biweeks with no reported cases", y = NULL)
sv(p3, "sparsity", 6, 2.6)

# ---- 4. the reconstructed susceptible fraction -------------------------------
su <- m %>% group_by(state, time) %>%
  summarise(frac = 100 * median(susc / population), .groups = "drop")

p4 <- ggplot(su, aes(time, frac, colour = state)) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = pal_state) +
  scale_x_continuous(breaks = seq(2008, 2020, 2)) +
  labs(x = "year", y = "susceptible, percent of population")
sv(p4, "susceptible", 7, 3.2)

cat("\nsusceptible fraction, by state:\n")
print(as.data.frame(m %>% group_by(state) %>%
                      summarise(min = round(100 * min(susc / population), 2),
                                median = round(100 * median(susc / population), 2),
                                max = round(100 * max(susc / population), 2), .groups = "drop")))

# ---- 5. coverage over the period ---------------------------------------------
cv <- m %>%
  select(state, time, MCV1 = v1, MCV2 = v2) %>%
  pivot_longer(c(MCV1, MCV2), names_to = "dose", values_to = "p") %>%
  group_by(state, time, dose) %>%
  summarise(mid = 100 * median(p, na.rm = TRUE),
            lo  = 100 * quantile(p, 0.25, na.rm = TRUE),
            hi  = 100 * quantile(p, 0.75, na.rm = TRUE), .groups = "drop")

p5 <- ggplot(cv, aes(time, mid, colour = dose, fill = dose)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.5) +
  # the x axis is repeated per panel; the y is held at 0 to 100 deliberately,
  # since coverage is a percentage and the two states should be read against
  # the same range
  facet_wrap(~ state, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = c(MCV1 = "#2166AC", MCV2 = "#7FCDBB")) +
  scale_fill_manual(values = c(MCV1 = "#2166AC", MCV2 = "#7FCDBB")) +
  scale_x_continuous(breaks = seq(2008, 2020, 4)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "year", y = "routine coverage, percent")
sv(p5, "coverage", 8, 3.4)

# ---- 6. cases and coverage together ------------------------------------------
# Both panels keep their own axes, so either can be read alone. The cases block
# takes more height, since its variation is the thing to see.
p6 <- p1 / p5 + plot_layout(heights = c(1.6, 1))
sv(p6, "cases_and_coverage", 9, 8.5)

cat("\nfigures in", FIG_DIR, "\n")
