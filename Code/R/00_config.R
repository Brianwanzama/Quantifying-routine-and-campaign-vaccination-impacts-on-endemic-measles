# 00_config.R

ROOT     <- "/.../.../.../"
DATA_DIR <- file.path(ROOT, "Data")
OUT_DIR  <- file.path(ROOT, "Derived")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TAB_DIR  <- file.path(OUT_DIR, "tables")

for (d in c(OUT_DIR, FIG_DIR, TAB_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- the states, and where each file sits ------------------------------------
# The two states name their files differently, so each is listed rather than
# built from a pattern. West Bengal has no campaign, hence NA.
STATES <- c("Uttar Pradesh", "West Bengal")

FILES <- list(
  `Uttar Pradesh` = list(
    coords   = "Uttar Pradesh/coordinates.csv",
    cases    = "Uttar Pradesh/Uttar Pradesh_cases.csv",
    births   = "Uttar Pradesh/Uttar Pradesh_live_births.csv",
    popn     = "Uttar Pradesh/Uttar Pradesh_popn.csv",
    mcv1     = "Uttar Pradesh/Uttar Pradesh_mcv1.csv",
    mcv2     = "Uttar Pradesh/Uttar Pradesh_mcv2.csv",
    campaign = "Uttar Pradesh/Uttar Pradesh_mr_campaign.csv"),
  `West Bengal` = list(
    coords   = "West Bengal/coordinates.csv",
    cases    = "West Bengal/Cases_west_bengal.csv",
    births   = "West Bengal/Live_births_west_bengal.csv",
    popn     = "West Bengal/popn_west_bengal.csv",
    mcv1     = "West Bengal/West Bengal_mcv1.csv",
    mcv2     = "West Bengal/West Bengal_mcv2.csv",
    campaign = NA_character_)
)

file_of <- function(state, what) {
  f <- FILES[[state]][[what]]
  if (is.na(f)) return(NA_character_)
  file.path(DATA_DIR, f)
}

# ---- the biweek clock --------------------------------------------------------
# 14 day blocks from 1 January, the count resetting each year, so biweek 7
# contains 1 April. t_global = (year - BASE_YEAR) * 26 + biweek.
BASE_YEAR        <- 2008L
BIWEEK_DAYS      <- 14L
BIWEEKS_PER_YEAR <- 26L


# April 2008 to March 2020. The first six biweeks of 2008 fall before the data
# and are dropped; biweeks 1 to 6 of later years are the January to March tail
# of an ongoing year and are kept.
STUDY_START_DATE  <- as.Date("2008-04-01")
STUDY_END_DATE    <- as.Date("2020-03-31")
FIRST_BIWEEK_KEPT <- 7L

# ---- rounding ----------------------------------------------------------------
# Monthly counts are apportioned across biweeks by day overlap, which gives
# fractions. Stochastic rounding preserves the total in expectation.
RNG_SEED <- 2026L

# ---- stocks ------------------------------------------------------------------
# An annual figure describes the year as a whole, so population and coverage
# are anchored at mid year and read onto biweek midpoints.
POP_ANCHOR_MONTH <- 7L
POP_ANCHOR_DAY   <- 2L
COVERAGE_IS_PERCENT <- TRUE

# ---- routine vaccination -----------------------------------------------------
E1  <- 0.84   
E12 <- 0.941   

# Coverage shifted so each dose meets the birth cohort reaching that age.
SHIFT_V1 <- 20L
SHIFT_V2 <- 32L

# ---- campaign ----------------------------------------------------------------
E_CAMPAIGN   <- 0.85
CAMPAIGN_LAG <- 1L     # a dose protects about 14 days later
KAPPA_CAP    <- 1.0    # cannot reach more than the whole age band

# ---- the mean susceptible level ----------------------------------------------

# The fallback is above the endemic equilibrium fraction 1 / R0, which for
# measles gives 0.056 to 0.083. Equilibrium is where transmission is exactly
# self sustaining, so it is a floor; where births replenish the pool faster
# than transmission draws it down the standing level sits above it.
SBAR_GRID_MAX <- 0.6
SBAR_FALLBACK <- 0.1
SBAR_EDGE_TOL <- 2L      # positions from the grid end that count as the edge

# ---- reconstruction ----------------------------------------------------------
MISSING_CASES_ARE_ZERO <- TRUE   # under reporting is absorbed by rho
ALPHA_FIXED <- NA                # NA estimates it

# ---- counterfactual runs -----------------------------------------------------

CF_PROCESS_NOISE  <- FALSE           # FALSE: replicates are deterministic and vary only
                                     # in refitted parameters, matching the figures
CF_N_BOOT         <- 1000L          # bootstrap replicates
CF_BLOCK_LENGTH   <- 26L            # residual block length, one year of biweeks
CF_START_YEARS_UP <- c(2012L, 2013L) # Uttar Pradesh from April of these years,
                                     # main reconstruction
CF_REFIT_YEARS_UP <- c(2013L)        # Uttar Pradesh reconstruction refitted from
                                     # April of this year by 06b_tsir_refit_up.R
CF_RERUN          <- TRUE             # 20_run_all.R: redo runs whose results exist
# Runs that feed Table 2, and the period each covers. West Bengal comes from the
# main run; Uttar Pradesh from the reconstruction refitted from 2013, because its
# reporting before late 2012 was too sparse (Fig 1).
CF_TABLE_RUN <- c("West Bengal"   = "scenarios",
                  "Uttar Pradesh" = "scenarios_refit_up_from2013")
CF_TABLE_PERIOD <- c("West Bengal"   = "April 2008 to March 2020",
                     "Uttar Pradesh" = "April 2013 to March 2020")

CF_STACK_RUN <- c("Uttar Pradesh" = "scenarios", "West Bengal" = "scenarios")

# What the scenario figure draws. "point": the deterministic simulation from the
# fitted model, without random variation
CF_STACK_SOURCE <- "median"

# ---- figures -----------------------------------------------------------------

# susceptible fraction, percent of population
SUSC_BREAKS <- c(-Inf, 3, 6, 10, Inf)
# The bands are named by their range rather than by a word, since the legend
# title already says what the quantity is and "high" says less than "6-10".
SUSC_LABELS <- c("<3", "3\u20136", "6\u201310", "\u226510")
SUSC_PALETTE <- c("#FEF0D9", "#FDCC8A", "#E34A33", "#900000")

# forecast incidence, cases per 100,000 over the forecast year
RATE_BREAKS <- c(-Inf, 5, 20, 50, Inf)
RATE_LABELS <- c("<5", "5\u201320", "20\u201350", "\u226550")
RATE_PALETTE <- c("#FEF0D9", "#FDCC8A", "#E34A33", "#900000")
CF_MAP_BREAKS <- c(-Inf, 10, 20, 30, Inf)

# How many of the highest bands get their districts named, counted from the
# top of the scale. 
MAP_LABEL_TOP_N <- 4L


MAP_LABEL_SIZE   <- 1.4
MAP_LABEL_COLOUR <- "black"

# how many biweeks at the end of the series the susceptible fraction averages
# over. One month, so a single noisy final value does not set the category.
SUSC_AVG_BIWEEKS <- 2L

# the districts drawn in the horizon figures, one per state
FIG_DISTRICTS <- c("Uttar Pradesh" = "Siddharth Nagar",
                   "West Bengal"   = "Murshidabad")
CF_FIG_DISTRICTS <- c("Uttar Pradesh" = "Siddharth Nagar",
                      "West Bengal"   = "Dakshin Dinajpur")

# ---- keys --------------------------------------------------------------------
# Two states, so district alone could collide. Every join keys on both.
UID_SEP  <- "__"
make_uid <- function(state, district) paste(state, district, sep = UID_SEP)

set.seed(RNG_SEED)
message("config loaded | ", length(STATES), " states | ",
        format(STUDY_START_DATE), " to ", format(STUDY_END_DATE))
