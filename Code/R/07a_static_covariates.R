# extract_covariates.R
#   IR file  : mother and household variables
#   KR file  : child variables, restricted to living children in the
#              standard age window for each indicator
#   All means are weighted by the individual sample weight.

library(haven)
library(dplyr)
library(readr)
library(labelled)

# ---- paths and settings ----------------------------------------------------
IR_PATH  <- "/.../.../NFHS5/IAIR7EFL.DTA"
KR_PATH  <- "/.../.../NFHS5/IAKR7EFL.DTA"
OUT_DIR  <- "/.../.../Data/NFHS5_State_Covariates"
AS_PERCENT <- TRUE          # write proportions as percentages

study_states <- c("uttar pradesh", "west bengal")

ir <- read_dta(IR_PATH)
kr <- read_dta(KR_PATH)
cat("IR rows:", nrow(ir), " KR rows:", nrow(kr), "\n")

# ============================================================
# IR covariates: mother and household
# ============================================================
ir_cov <- ir %>%
  mutate(
    state    = as.character(to_factor(v024)),
    district = as.character(to_factor(sdist)),
    wt       = v005 / 1000000,
    
    # maternal schooling: ten or more completed years
    women10plus = if_else(v133 >= 10 & v133 < 95, 1, 0),
    
    # lowest wealth quintile
    poorest = if_else(v190 == 1, 1, 0),
    
    # urban residence
    urban = if_else(v025 == 1, 1, 0),
    
    # scheduled caste or scheduled tribe
    sc_st = case_when(
      s116 %in% c(1, 2) ~ 1,
      s116 %in% c(3, 4) ~ 0,
      TRUE              ~ NA_real_
    )
  ) %>%
  filter(state %in% study_states) %>%
  group_by(state, district) %>%
  summarise(
    women10plus  = weighted.mean(women10plus, wt, na.rm = TRUE),
    poorest      = weighted.mean(poorest,     wt, na.rm = TRUE),
    hhsize       = weighted.mean(v136,        wt, na.rm = TRUE),
    urbanization = weighted.mean(urban,       wt, na.rm = TRUE),
    sc_st        = weighted.mean(sc_st,       wt, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# IR: institutional birth, from the most recent birth (m15_1)
#   Codes 21 to 33 are health facilities. 10 to 13 are home.
#   Restricted to women with a birth in the five years before the survey.
# ============================================================
ir_delivery <- ir %>%
  mutate(
    state    = as.character(to_factor(v024)),
    district = as.character(to_factor(sdist)),
    wt       = v005 / 1000000,
    inst_birth = case_when(
      m15_1 >= 20 & m15_1 <= 33 ~ 1,
      m15_1 >= 10 & m15_1 <= 13 ~ 0,
      m15_1 == 96               ~ 0,
      TRUE                      ~ NA_real_
    ),
    # maternal immunisation: two or more tetanus injections during the pregnancy
    # m1 codes 0 to 7 are counts, 8 is don't know
    mat_tt = case_when(
      m1_1 >= 2 & m1_1 <= 7 ~ 1,
      m1_1 >= 0 & m1_1 <  2 ~ 0,
      TRUE                  ~ NA_real_
    )
  ) %>%
  filter(state %in% study_states, !is.na(inst_birth)) %>%
  group_by(state, district) %>%
  summarise(
    institutional_birth   = weighted.mean(inst_birth, wt, na.rm = TRUE),
    maternal_immunisation = weighted.mean(mat_tt,     wt, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# KR covariates: children
#   b5 == 1 keeps living children only.
#   Vaccination indicators use children aged 12 to 23 months, which is the
#   standard denominator. Younger children have not reached the age for the
#   dose and would otherwise be counted as unvaccinated.
#   Vaccination codes: 0 no, 1 date on card, 2 mother reported,
#   3 marked on card, 4 identified, 8 don't know.
# ============================================================
kr_base <- kr %>%
  mutate(
    state    = as.character(to_factor(v024)),
    district = as.character(to_factor(sdist)),
    wt       = v005 / 1000000
  ) %>%
  filter(state %in% study_states, b5 == 1)

vax_yes <- function(x) case_when(
  x %in% c(1, 2, 3, 4) ~ 1,
  x == 0               ~ 0,
  TRUE                 ~ NA_real_
)

# vaccination block, children 12 to 23 months
kr_vax <- kr_base %>%
  filter(hw1 >= 12, hw1 <= 23) %>%
  mutate(
    mcv1        = vax_yes(h9),
    penta1      = vax_yes(h51),
    penta3      = vax_yes(h53),
    bcg         = vax_yes(h2),
    measles_zero = if_else(mcv1 == 0, 1, if_else(mcv1 == 1, 0, NA_real_)),
    penta_zero   = if_else(penta1 == 0, 1, if_else(penta1 == 1, 0, NA_real_)),
    # full immunisation: BCG, three pentavalent, three polio, measles
    full_imm = if_else(
      vax_yes(h2) == 1 &
        vax_yes(h51) == 1 & vax_yes(h52) == 1 & vax_yes(h53) == 1 &
        vax_yes(h4)  == 1 & vax_yes(h6)  == 1 & vax_yes(h8)  == 1 &
        vax_yes(h9)  == 1, 1, 0)
  ) %>%
  group_by(state, district) %>%
  summarise(
    measles_zero = weighted.mean(measles_zero, wt, na.rm = TRUE),
    penta_zero   = weighted.mean(penta_zero,   wt, na.rm = TRUE),
    full_imm     = weighted.mean(full_imm,     wt, na.rm = TRUE),
    n_children_vax = n(),
    .groups = "drop"
  )

# nutrition, children under five, WHO z-scores in hw70 (height for age)
# values 9996 to 9998 are flags and are excluded
kr_nutr <- kr_base %>%
  mutate(
    haz = if_else(hw70 < 9000, hw70 / 100, NA_real_),
    stunted = if_else(haz < -2, 1, if_else(!is.na(haz), 0, NA_real_))
  ) %>%
  filter(!is.na(stunted)) %>%
  group_by(state, district) %>%
  summarise(
    stunting = weighted.mean(stunted, wt, na.rm = TRUE),
    .groups = "drop"
  )

# breastfeeding: put to the breast within one hour of birth
# m34 codes: 0 immediately, 100 within first hour, 101 to 199 hours,
# 201 to 299 days. Early initiation is 0 or 100.
kr_bf <- kr_base %>%
  filter(hw1 <= 23) %>%
  mutate(
    early_bf = case_when(
      m34 %in% c(0, 100)      ~ 1,
      m34 > 100 & m34 < 300   ~ 0,
      TRUE                    ~ NA_real_
    )
  ) %>%
  filter(!is.na(early_bf)) %>%
  group_by(state, district) %>%
  summarise(
    early_breastfeeding = weighted.mean(early_bf, wt, na.rm = TRUE),
    .groups = "drop"
  )

# vitamin A in the last six months, children 6 to 59 months
kr_vita <- kr_base %>%
  filter(hw1 >= 6, hw1 <= 59) %>%
  mutate(
    vita = case_when(
      h34 == 1 ~ 1,
      h34 == 0 ~ 0,
      TRUE     ~ NA_real_
    )
  ) %>%
  filter(!is.na(vita)) %>%
  group_by(state, district) %>%
  summarise(
    vitamin_a = weighted.mean(vita, wt, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# Merge
# ============================================================
district_covariates <- ir_cov %>%
  left_join(ir_delivery, by = c("state", "district")) %>%
  left_join(kr_vax,      by = c("state", "district")) %>%
  left_join(kr_nutr,     by = c("state", "district")) %>%
  left_join(kr_vita,     by = c("state", "district")) %>%
  left_join(kr_bf,       by = c("state", "district"))

if (AS_PERCENT) {
  pct_cols <- c("women10plus", "poorest", "urbanization", "sc_st",
                "institutional_birth", "maternal_immunisation",
                "measles_zero", "penta_zero", "full_imm",
                "stunting", "vitamin_a", "early_breastfeeding")
  district_covariates <- district_covariates %>%
    mutate(across(all_of(pct_cols), ~ .x * 100))
}

# ---- report ----------------------------------------------------------------
cat("\ndistricts per state:\n")
print(table(district_covariates$state))
cat("\nmissing values per covariate:\n")
print(colSums(is.na(district_covariates)))
cat("\nsummary:\n")
print(summary(district_covariates[, -(1:2)]))

cat("\nchildren aged 12 to 23 months per district:\n")
print(summary(district_covariates$n_children_vax))
low <- district_covariates %>% filter(n_children_vax < 25) %>%
  select(state, district, n_children_vax, measles_zero, penta_zero, full_imm)
if (nrow(low)) {
  cat("\ndistricts with fewer than 25 children, vaccination covariates imprecise:\n")
  print(as.data.frame(low))
} else {
  cat("\nno district has fewer than 25 children in the vaccination window\n")
}

# ============================================================
# Export
# ============================================================
dir.create(OUT_DIR, showWarnings = FALSE)
for (st in study_states) {
  out <- district_covariates %>% filter(state == st) %>% arrange(district)
  write_csv(out, file.path(OUT_DIR,
                           paste0(gsub(" ", "_", st), "_district_covariates.csv")))
}
write_csv(district_covariates, file.path(OUT_DIR, "all_states_district_covariates.csv"))
cat("\nwrote", length(study_states), "state files plus a combined file to", OUT_DIR, "\n")

# ============================================================
# Collinearity check on the covariate set
# ============================================================
cov_names <- c("women10plus", "poorest", "hhsize", "urbanization", "sc_st",
               "institutional_birth", "maternal_immunisation",
               "measles_zero", "penta_zero", "full_imm",
               "stunting", "vitamin_a", "early_breastfeeding")
num <- district_covariates[, cov_names]

cat("\n=== Pearson correlation ===\n")
cc <- cor(num, use = "pairwise.complete.obs")
print(round(cc, 2))

cat("\n=== pairs with |r| >= 0.7 ===\n")
nm <- colnames(cc); found <- FALSE
for (i in 1:(ncol(cc) - 1)) for (j in (i + 1):ncol(cc)) {
  if (is.finite(cc[i, j]) && abs(cc[i, j]) >= 0.7) {
    cat(sprintf("  %-22s <-> %-22s r = %+.2f\n", nm[i], nm[j], cc[i, j]))
    found <- TRUE
  }
}
if (!found) cat("  none\n")

cat("\n=== VIF ===\n")
comp <- num[complete.cases(num), ]
vif1 <- function(v) {
  f <- as.formula(paste0("`", v, "` ~ ",
                         paste(sprintf("`%s`", setdiff(names(comp), v)), collapse = " + ")))
  r2 <- summary(lm(f, data = comp))$r.squared
  if (r2 >= 1) Inf else 1 / (1 - r2)
}
print(round(sort(sapply(names(comp), vif1), decreasing = TRUE), 2))
cat("\nVIF above 5 indicates redundancy with the other covariates.\n")

