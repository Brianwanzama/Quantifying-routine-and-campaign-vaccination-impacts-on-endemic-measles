# 16_fig_explain.R
# What the network attends to, from the grouped Shapley attributions.
# Input  : explain/shap_attributions_k<K>.csv, shap_importance_k<K>.csv
# Output : figures/explain/*.pdf and .png

source("/.../.../.../Code/R/00_config.R")
library(tidyverse)
library(scales)

EXPLAIN <- file.path(OUT_DIR, "sfnn", "explain")
OUT_FIG <- file.path(FIG_DIR, "explain")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

KS      <- as.integer(strsplit(Sys.getenv("MS_KS", "1,4,6,12"), ",")[[1]])
VARIANT <- Sys.getenv("MS_VARIANT", "with_nfhs")

# the stack shows one horizon; the bars show all of them
K_STACK <- as.integer(Sys.getenv("MS_K_STACK", as.character(max(KS))))

# from the config, so the stack and the horizon figure show the same districts
PICK <- FIG_DISTRICTS
if (nzchar(Sys.getenv("MS_UP"))) PICK[["Uttar Pradesh"]] <- Sys.getenv("MS_UP")
if (nzchar(Sys.getenv("MS_WB"))) PICK[["West Bengal"]]   <- Sys.getenv("MS_WB")

cat("variant", VARIANT, "| horizons", paste(KS, collapse = ", "),
    "| the stack shows k =", K_STACK, "\n")

# ---- the groups, and how they are named on a figure --------------------------
GROUP_LABELS <- c(
  own            = "Local incidence",
  susc           = "Susceptibles",
  nearest        = "Nearby incidence",
  core           = "Large districts",
  births         = "Births",
  pop            = "Population",
  vacc           = "MCV1 and MCV2",
  campaign       = "MR campaign",
  nfhs_supply    = "NFHS: susceptible supply",
  nfhs_contact   = "NFHS: contact structure",
  nfhs_access    = "NFHS: service access",
  nfhs_nutrition = "NFHS: nutrition and outreach")

pal <- c("Local incidence"              = "#1B7837",
         "Susceptibles"                 = "#762A83",
         "Nearby incidence"             = "#0072B2",
         "Large districts"              = "#56B4E9",
         "Births"                       = "#8C510A",
         "Population"                   = "#C994C7",
         "MCV1 and MCV2"                = "#B15928",
         "MR campaign"                  = "#E31A1C",
         "NFHS: susceptible supply"     = "#E69F00",
         "NFHS: contact structure"      = "#FDBF6F",
         "NFHS: service access"         = "#FFD92F",
         "NFHS: nutrition and outreach" = "#BFBFBF")

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(panel.border     = element_rect(fill = NA, colour = "black",
                                          linewidth = 0.5),
          strip.background = element_rect(fill = "grey92", colour = "black",
                                          linewidth = 0.4),
          strip.text       = element_text(face = "bold", size = base_size - 2),
          axis.text        = element_text(colour = "black",
                                          size = base_size - 2.5),
          axis.title       = element_text(size = base_size - 1.5,
                                          colour = "grey15"),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          legend.text      = element_text(size = base_size - 2.5),
          legend.key.size  = unit(9, "pt"),
          plot.title       = element_blank(),
          plot.caption     = element_blank())
}

# ---- read --------------------------------------------------------------------
read_attr <- function(k) {
  f <- file.path(EXPLAIN, paste0("shap_attributions_k", k, ".csv"))
  if (!file.exists(f)) { cat("  no attributions for k =", k, "\n"); return(NULL) }
  read_csv(f, show_col_types = FALSE) %>% filter(variant == VARIANT)
}

att <- bind_rows(lapply(KS, read_attr))
if (!nrow(att)) stop("no attributions for variant ", VARIANT)
KS <- sort(unique(att$k))
cat("horizons available:", paste(KS, collapse = ", "), "\n")

# the explain script prefixes the identifiers, since one group is also called
# state and a bare column would have overwritten its attribution
have <- intersect(names(GROUP_LABELS), names(att))
cat("groups present:", paste(have, collapse = ", "), "\n")

sel <- tibble(state = STATES, district = unname(PICK[STATES])) %>%
  mutate(uid = make_uid(state, district))

# ---- 1. the stack, composition through the series ----------------------------
if (!K_STACK %in% KS)
  cat("k =", K_STACK, "is not among the horizons read; the stack will be empty\n")

long <- att %>%
  filter(id_uid %in% sel$uid, k == K_STACK) %>%
  select(uid = id_uid, state = id_state, district = id_district,
         time = id_time, k, all_of(have)) %>%
  pivot_longer(all_of(have), names_to = "group", values_to = "shap")

if (!nrow(long)) {
  cat("\nneither chosen district appears in the attributions.\n")
  cat("available:", paste(head(unique(att$id_district), 10), collapse = ", "), "\n")
} else {
  share <- long %>%
    group_by(k, uid, time) %>%
    mutate(share = abs(shap) / sum(abs(shap), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(group = factor(unname(GROUP_LABELS[group]),
                          levels = unname(GROUP_LABELS)),
           panel_d = paste0(district, " (", state, ")"))
  
  p1 <- ggplot(share, aes(time, share, fill = group)) +
    geom_area(position = "stack", alpha = 0.9, colour = NA) +
    facet_wrap(~ panel_d, ncol = 2) +
    scale_fill_manual(values = pal, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_x_continuous(breaks = pretty_breaks(5)) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = "Year", y = "Share of absolute attribution") +
    theme_pub()
  
  f1 <- file.path(OUT_FIG, paste0("shap_stack_k", K_STACK, "_", VARIANT))
  ggsave(paste0(f1, ".png"), p1, width = 10, height = 4.5, dpi = 400,
         bg = "white")
  ggsave(paste0(f1, ".pdf"), p1, width = 10, height = 4.5, device = cairo_pdf)
  cat("saved", basename(f1), "\n")
  
  cat("\nmean share by group at k =", K_STACK, ", the two districts:\n")
  print(as.data.frame(share %>%
                        group_by(district, group) %>%
                        summarise(m = round(100 * mean(share, na.rm = TRUE), 1),
                                  .groups = "drop") %>%
                        pivot_wider(names_from = district, values_from = m)))
}

# ---- 2. the bars, mean importance by group and horizon -----------------------
read_imp <- function(k) {
  f <- file.path(EXPLAIN, paste0("shap_importance_k", k, ".csv"))
  if (!file.exists(f)) return(NULL)
  read_csv(f, show_col_types = FALSE)
}

imp <- bind_rows(lapply(KS, read_imp))
if (!is.null(imp) && nrow(imp)) {
  d <- imp %>%
    filter(group %in% names(GROUP_LABELS)) %>%
    mutate(label = factor(unname(GROUP_LABELS[group]),
                          levels = rev(unname(GROUP_LABELS))),
           panel_k = factor(paste0("k = ", k), levels = paste0("k = ", KS)))
  
  p2 <- ggplot(d, aes(share, label, fill = variant)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    facet_wrap(~ panel_k, nrow = 1) +
    scale_fill_manual(values = c(with_nfhs = "#2166AC", no_nfhs = "#B2182B"),
                      labels = c(with_nfhs = "with covariates",
                                 no_nfhs = "without")) +
    scale_x_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Share of absolute attribution", y = NULL) +
    theme_pub()
  
  f2 <- file.path(OUT_FIG, "shap_importance_by_horizon")
  ggsave(paste0(f2, ".png"), p2, width = 14, height = 5, dpi = 400, bg = "white")
  ggsave(paste0(f2, ".pdf"), p2, width = 14, height = 5, device = cairo_pdf)
  cat("saved", basename(f2), "\n")
  
  
  w <- d %>% filter(variant == "with_nfhs")
  if (nrow(w)) {
    cat("\nthe four covariate pathways together, by horizon:\n")
    print(as.data.frame(w %>%
                          mutate(is_nfhs = str_starts(group, "nfhs_")) %>%
                          group_by(panel_k) %>%
                          summarise(nfhs_share = round(sum(share[is_nfhs]), 1),
                                    largest_other = group[which.max(ifelse(is_nfhs, -Inf, share))],
                                    largest_other_share = round(max(share[!is_nfhs]), 1),
                                    .groups = "drop")))
  }
  
  write_csv(d, file.path(TAB_DIR, "shap_importance_all_horizons.csv"))
}

cat("\nfigures in", OUT_FIG, "\n")

