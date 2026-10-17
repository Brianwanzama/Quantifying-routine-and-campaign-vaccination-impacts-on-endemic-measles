# Quantifying-routine-and-campaign-vaccination-impacts-on-endemic-measles
**Extended TSIR modelling and spatial neural network forecasting**

Code for the manuscript of the same title, by Brian Wanzama and Siuli
Mukhopadhyay (Department of Mathematics, IIT Bombay).

The analysis covers 95 districts (75 in Uttar Pradesh, 20 in West Bengal) from
April 2008 to March 2020 on a biweekly time step. It has three parts.

1. **Susceptible reconstruction.** A time-series SIR (TSIR) model is extended to
   two routine doses of measles-containing vaccine and the 2018–19
   measles–rubella campaign in Uttar Pradesh.
2. **Counterfactual scenarios.** The fitted model estimates the infections the
   first dose, the second dose and the campaign each averted, with block
   bootstrap intervals.
3. **Forecasting.** A spatial feature neural network (SFNN), following Madden et
   al. (2024), forecasts reported cases 1, 4, 6 and 12 biweeks ahead. It is
   compared with a rolling TSIR forecast, tested with and without 13 district
   indicators from NFHS-5, and explained with grouped Shapley values.

The derivation of the reconstruction and the campaign factor is in
`S1_Appendix.pdf`.

---

## Folder layout

```
Manuscript/
├── Code/
│   ├── R/          00_config.R … 31_plos_figure_files.R
│   └── python/     09_ … 13_, 17_ (.py and .sh)
├── Data/           raw inputs, one folder per state (see below)
└── Derived/        everything the code writes
    ├── figures/    (with compare/ and explain/ subfolders)
    ├── tables/     (with supplement/)
    ├── sfnn/       prefit/, tuning/, final/, explain/
    ├── scenarios*/ one folder per counterfactual run
    ├── tsir_forecast/
    ├── supporting_information/   S1_Appendix, S1–S2 Fig, S1–S6 Table (PDF)
    └── plos_figures/             Fig1.tif … Fig9.tif
```

Every path and parameter is set in `00_config.R`. Set `ROOT` there to the
`Manuscript` folder on your machine. Scripts that are run on their own source the
config by its full path, so change that line too if the folder moves.

## Input data

| File (per state, under `Data/<state>/`) | Content | Source |
|---|---|---|
| `*_cases.csv` | Monthly reported measles cases by district | Health Management Information System (HMIS), MoHFW |
| `*_live_births.csv` | Monthly live births by district | HMIS |
| `*_popn.csv` | Annual district population, 2012–2031, extrapolated back to 2008 | IIPS district projections (Dhar 2022) |
| `*_mcv1.csv`, `*_mcv2.csv` | Annual first- and second-dose coverage | IHME subnational vaccine coverage estimates |
| `Uttar Pradesh_mr_campaign.csv` | Campaign start and end dates, children vaccinated, population aged 9 months–15 years, by district | District campaign reports *(name the source)* |
| `coordinates.csv` | District centres; also the reference list of district names | GADM 3.6 |
| `Derived/all_states_district_covariates.csv` | 13 district indicators | NFHS-5 (2019–21) unit-level data, DHS Program |

The NFHS-5 microdata are available from the DHS Program after registration and
are not redistributed here. West Bengal held no campaign in the study period,
so it has no campaign file.

## Software

- **R 4.4.1**, with `tsiR` 0.4.3, `sf` 1.0.24, `tidyverse`, `lubridate`,
  `scales`, `patchwork`, `ragg`, `gridExtra`, `gtable`, `magick` and `png`.
- **Python 3**, with PyTorch, Ray Tune (ASHA scheduler), Captum, NumPy and pandas.
  *Record the Python, PyTorch, Ray and Captum versions here.*
- The network was trained on a Linux server (Ubuntu 24.04, Intel Xeon Gold
  6226R, two NVIDIA RTX A5000 GPUs). The R scripts run on Windows or Linux.

## Running the pipeline

Run the scripts in this order. Each prints what it read and wrote.

### 1. Data preparation (R)

| Script | Does | Writes |
|---|---|---|
| `01a_audit_names.R` | Checks district names in every file against the reference list. Run once; writes nothing | — |
| `01_district_names.R` | Name harmonisation, sourced by the other scripts | — |
| `02_flows_to_biweekly.R` | Monthly cases and births to biweeks, by day overlap | `cases_biweekly.csv`, `births_biweekly.csv` |
| `03_stocks_to_biweekly.R` | Annual population and coverage interpolated to biweeks | `popn_biweekly.csv`, `coverage_biweekly.csv` |
| `04_campaign_factors.R` | Campaign survival factor Φ for every district and biweek | `campaign_factors.csv` |

### 2. Reconstruction and counterfactuals (R)

| Script | Does | Writes |
|---|---|---|
| `05_tsir_engine.R` | The `tsiR` fitting function, extended with the campaign term. Sourced, not run | — |
| `06_tsir_reconstruction.R` | Reconstruction for every district, April 2008–March 2020 | `tsir_panel.csv`, `tsir_fit_log.csv`, `tsir_beta.csv` |
| `06b_tsir_refit_up.R` | The same for Uttar Pradesh, refitted from April 2013 (used for its counterfactual) | `*_up_from2013.csv` |
| `20_run_all.R` | Runs `20_scenario_analysis.R` for every run in `CF_RUNS`: the main run, Uttar Pradesh from later starts, and the 2013 refit. Bootstrap of 1,000 replicates | `scenarios*/` |
| `21_fig_scenarios.R` | Table 2, the calibration table, Figs 3–5 and S1 Fig | `tables/`, `figures/fig_cf_*` |

### 3. Forecasting (Python, on the GPU machine)

| Script | Does |
|---|---|
| `07_build_master.R` (R) | One long file of every district and biweek for the network: `sfnn_master.csv` |
| `09_sfnn_loader.py` | Feature matrices for each horizon and lag window |
| `10_sfnn_functions.py` | Network, data wrapper and training loops; imported by the others |
| `11_run_tuning.sh` → `11_sfnn_raytune.py` | Hyperparameter search per horizon, with and without the covariates |
| `12_run_training.sh` → `12_sfnn_train.py` | Final models, three seeds (0, 1, 2), test-period predictions |
| `13_run_explain.sh` → `13_sfnn_explain.py` | Grouped Shapley values on 3,000 test rows |

### 4. TSIR forecast and comparisons (R)

| Script | Does |
|---|---|
| `22_tsir_run.R`, `23_tsir_process.R` | Rolling-origin TSIR forecast with a 130-biweek window, merged to the SFNN rows |
| `24_tsir_window_sensitivity.R` | The same at other windows |
| `25_compare_sfnn_tsir.R` | SFNN against TSIR by district (Fig 7) |
| `15_district_performance.R` | Correlation and RMSE by district and horizon (Table 3) |
| `19_nfhs_comparison.R` | With and without the covariates, Wilcoxon test (Table 4, S2 Fig) |
| `27_persistence_benchmark.R` | Both methods against carrying the last value forward |

### 5. Figures, tables and numbers (R)

| Script | Writes |
|---|---|
| `08_explore_plots.R` | Fig 1 |
| `18_map_forecast.R` | Fig 2 (susceptible fraction map) |
| `14_fig_horizons.R` | Fig 6 |
| `16_fig_explain.R` | Figs 8 and 9 |
| `28_manuscript_numbers.R` | Every number quoted in the text: `tables/manuscript_numbers.txt` |
| `29_supplementary_tables.R` | S1–S6 Table as PDFs |
| `30_supporting_files.R` | All supporting files, named for upload |
| `31_plos_figure_files.R` | Fig1.tif–Fig9.tif (LZW TIFF, at most 7.5 × 8.75 in), with a check table |

## Where each manuscript item comes from

| Item | Script |
|---|---|
| Table 1 (covariates by pathway) | written in the manuscript; indicators from `07_build_master.R` |
| Table 2 (infections averted) | `21_fig_scenarios.R` |
| Table 3 (forecast performance) | `15_district_performance.R`, `28_manuscript_numbers.R` |
| Table 4 (with and without covariates) | `19_nfhs_comparison.R` |
| Fig 1 | `08_explore_plots.R` |
| Figs 2 and 4 (maps) | `18_map_forecast.R`, `21_fig_scenarios.R` |
| Figs 3 and 5 | `21_fig_scenarios.R` |
| Fig 6 | `14_fig_horizons.R` |
| Fig 7 | `25_compare_sfnn_tsir.R` |
| Figs 8 and 9 | `16_fig_explain.R` |
| S1 Fig | `21_fig_scenarios.R` |
| S2 Fig | `19_nfhs_comparison.R` |
| S1–S6 Table | `29_supplementary_tables.R` |
| All numbers in the text | `28_manuscript_numbers.R` |

## Key settings (`00_config.R`)

| Setting | Value | Meaning |
|---|---|---|
| `E1`, `E12` | 0.84, 0.941 | Field effectiveness of one and two doses |
| `SHIFT_V1`, `SHIFT_V2` | 20, 32 biweeks | Coverage shifted to the age at each dose |
| `E_CAMPAIGN` | 0.85 | Effectiveness of a campaign dose |
| `CAMPAIGN_LAG` | 1 biweek | Time from dose to protection |
| `SBAR_GRID_MAX`, `SBAR_FALLBACK` | 0.6, 0.1 | Range searched for the mean susceptible level, and the value used where it is not identified |
| `CF_N_BOOT`, `CF_BLOCK_LENGTH` | 1,000, 26 | Bootstrap replicates and block length |
| `CF_PROCESS_NOISE` | `FALSE` | Replicates are deterministic and differ only in refitted parameters |
| `CF_TABLE_RUN` | WB: main run; UP: 2013 refit | Runs that feed Table 2 |
| `FIG_DISTRICTS` | Siddharth Nagar, Murshidabad | Districts drawn in Figs 5, 6 and 9 |
| `RNG_SEED` | 2026 | Seed for rounding and the bootstrap |

## Notes

- `06_tsir_reconstruction.R` is the reference reconstruction and should not be
  changed. `06b` is a copy that only restricts Uttar Pradesh to April 2013 onward.
- `I_t = ρ C_t`: ρ is the reporting rate (ρ > 1 under-reporting). The column
  `rho` in the output files holds the fraction reported, 1/ρ.
- `17_sfnn_forecast.py` forecasts beyond March 2020. It is not used in the
  manuscript.
- `20b_scenario_bootstrap.R`, `14_fig_train_test.R`, `26_tables.R` and
  `29_supplementary_tables_latex.R` are earlier or alternative versions and are
  not needed to reproduce the manuscript.

## Citation

Wanzama B, Mukhopadhyay S. Quantifying routine and campaign vaccination impacts
on endemic measles: Extended TSIR modelling and spatial neural network
forecasting. *Manuscript in preparation.*

Contact: Siuli Mukhopadhyay, siuli@math.iitb.ac.in
