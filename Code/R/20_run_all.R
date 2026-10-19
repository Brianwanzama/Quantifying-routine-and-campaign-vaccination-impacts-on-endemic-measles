# 20_run_all.R
# Runs the counterfactual analysis (script 20) for every run the manuscript
# reports, each into its own folder:
# Settings, all in 00_config.R:
#   CF_RUNS      which of the three to do, e.g. c("up_from", "up_refit") to redo
#                only the Uttar Pradesh runs; all three if not set
#   CF_RERUN     TRUE redoes runs whose results already exist
#   CF_N_BOOT    bootstrap replicates


source("/.../.../.../Code/R/00_config.R")

CODE_DIR <- file.path(ROOT, "Code", "R")
F20      <- file.path(CODE_DIR, "20_scenario_analysis.R")
if (!file.exists(F20)) stop("script 20 not found at ", F20, ": check ROOT in 00_config.R")

RUNS <- if (exists("CF_RUNS")) CF_RUNS else c("main", "up_from", "up_refit")
ENV  <- c("MS_CF_START_FY_UP", "MS_CF_STATES", "MS_RECON_START_FY_UP", "MS_RECON_TAG")
clear <- function() Sys.unsetenv(ENV)
done <- function(folder) file.exists(file.path(OUT_DIR, folder, "bootstrap_state.csv"))

run_log <- list()
attempt <- function(label, folder, env = character()) {
  if (!CF_RERUN && done(folder)) {
    run_log[[label]] <<- "skipped, results exist (CF_RERUN is FALSE)"
    cat("\n", label, ": results exist, skipped\n", sep = "")
    return(invisible())
  }
  cat("\n######## ", label, " ########\n", sep = "")
  clear(); Sys.setenv(MS_N_BOOT = CF_N_BOOT)
  if (length(env)) do.call(Sys.setenv, as.list(env))
  t0  <- Sys.time()
  res <- tryCatch({ source(F20, local = new.env()); "completed" },
                  error = function(e) paste("FAILED:", conditionMessage(e)))
  clear()
  if (res == "completed" && !done(folder)) res <- "FAILED: no results written"
  run_log[[label]] <<- sprintf("%s (%.1f min)", res, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

if ("main" %in% RUNS)
  attempt("main: both states from 2008", "scenarios")

if ("up_from" %in% RUNS)
  for (fy in CF_START_YEARS_UP)
    attempt(paste("Uttar Pradesh from", fy), paste0("scenarios_up_from", fy),
            c(MS_CF_START_FY_UP = fy, MS_CF_STATES = "Uttar Pradesh"))

if ("up_refit" %in% RUNS)
  for (fy in CF_REFIT_YEARS_UP) {
    if (!file.exists(file.path(OUT_DIR, paste0("tsir_panel_up_from", fy, ".csv")))) {
      run_log[[paste("Uttar Pradesh refitted from", fy)]] <-
        "FAILED: no refitted reconstruction; run 06b_tsir_refit_up.R first"
      next
    }
    attempt(paste("Uttar Pradesh refitted from", fy), paste0("scenarios_refit_up_from", fy),
            c(MS_RECON_TAG = paste0("_up_from", fy), MS_CF_STATES = "Uttar Pradesh"))
  }

clear()
cat("\n================ counterfactual runs ================\n")
for (k in names(run_log)) cat(sprintf("  %-40s %s\n", k, run_log[[k]]))
if (any(grepl("^FAILED", unlist(run_log))))
  cat("\nat least one run failed; its results are not current\n") else
  cat("\nall runs current; now run 21, 29 and 28\n")
