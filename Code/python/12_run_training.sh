#!/usr/bin/env bash
# =============================================================================
# 12_run_training.sh
#
# Trains the final model at each horizon, both variants, from the
# configurations the search settled on.
#
# Both variants are trained in one call per horizon, so they share the data and
# the seeds and differ only in whether the covariates were present.
#
# Several seeds, because a single fit gives one number and a small difference
# between variants can be initialisation noise rather than anything real.
#
# Usage
#   bash 12_run_training.sh
# =============================================================================
set -euo pipefail

ROOT="/.../.../.../"
CODE="$ROOT/Code/python"
PREFIT="$ROOT/Derived/sfnn/prefit"
TUNING="$ROOT/Derived/sfnn/tuning"
RESULTS="$ROOT/Derived/sfnn/final"
LOGS="$ROOT/Derived/sfnn/logs"

KS=(1 4 6 12 20 26)
SEEDS=3
EPOCHS=300
PATIENCE=15

mkdir -p "$RESULTS" "$LOGS"
cd "$CODE"

STAMP=$(date +%Y%m%d_%H%M%S)
RUNLOG="$LOGS/training_${STAMP}.log"

echo "checking the tuning results"
missing=0
for k in "${KS[@]}"; do
  for v in with_nfhs no_nfhs; do
    f="$TUNING/$v/tuning_k${k}.csv"
    if [[ -f "$f" ]]; then
      echo "  k=$k $v      found"
    else
      echo "  k=$k $v      MISSING"
      missing=1
    fi
  done
done
if (( missing )); then
  echo ""
  echo "run the search first for the studies above"
  exit 1
fi

{
  echo "started $(date)"
  echo "horizons ${KS[*]} | seeds $SEEDS | epoch cap $EPOCHS"
  for k in "${KS[@]}"; do
    echo ""
    echo "##############################################################"
    echo "# horizon k=$k"
    echo "##############################################################"
    python 12_sfnn_train.py \
      --k "$k" \
      --seeds "$SEEDS" --epochs "$EPOCHS" --patience "$PATIENCE" \
      --prefit-dir "$PREFIT" --tuning-dir "$TUNING" --results-dir "$RESULTS"
  done
  echo ""
  echo "finished $(date)"
} 2>&1 | tee "$RUNLOG"

echo ""
echo "=============================================================="
echo "across every horizon"
echo "=============================================================="
python - <<'PY' 2>&1 | tee -a "$RUNLOG"
import glob, os
import pandas as pd

RESULTS = "/home/brain/Msc_project/Manuscript/Derived/sfnn/final"
files = sorted(glob.glob(os.path.join(RESULTS, "summary_k*.csv")))
if not files:
    print("no summaries found")
    raise SystemExit

s = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

cols = [c for c in ["k", "variant", "t_lag", "n_features", "val_mse_mean",
                    "test_mse_mean", "test_mse_sd", "rmse_counts_mean",
                    "correlation_mean", "sd_ratio_mean", "skill_mean",
                    "skill_vs_persistence_mean", "best_epoch_mean"]
        if c in s.columns]
print(s[cols].sort_values(["k", "variant"]).round(4).to_string(index=False))

print("\nthe test error is a held out estimate: the architecture came from")
print("validation and training stopped there too.")

print("\nwhat the static covariates changed:")
for k in sorted(s["k"].unique()):
    sub = s[s.k == k]
    if sub["variant"].nunique() < 2:
        continue
    w = sub[sub.variant == "with_nfhs"].iloc[0]
    n = sub[sub.variant == "no_nfhs"].iloc[0]
    d = 100 * (n.test_mse_mean - w.test_mse_mean) / n.test_mse_mean
    tol = 2 * max(w.test_mse_sd, n.test_mse_sd)
    verdict = ("within the seed variation"
               if abs(w.test_mse_mean - n.test_mse_mean) < tol
               else ("lower error with them" if d > 0 else "lower error without"))
    print(f"  k={k:<3d} {d:+6.2f} percent   {verdict}")

flat = s[s.sd_ratio_mean < 0.3]
if len(flat):
    print("\npredictions varying far less than the data:")
    for _, r in flat.iterrows():
        print(f"  k={int(r.k)} {r.variant}: spread {r.sd_ratio_mean:.3f}, "
              f"correlation {r.correlation_mean:+.3f}")
    print("  the error there is not informative about forecasting skill")

if "skill_vs_persistence_mean" in s.columns:
    beaten = s[s.skill_vs_persistence_mean < 0]
    if len(beaten):
        print("\nbeaten by carrying the last observed value forward:")
        for _, r in beaten.iterrows():
            print(f"  k={int(r.k)} {r.variant}: "
                  f"{r.skill_vs_persistence_mean:+.3f}")

good = s[(s.skill_mean > 0.1) & (s.sd_ratio_mean > 0.4)]
if "skill_vs_persistence_mean" in s.columns:
    good = good[good.skill_vs_persistence_mean > 0]
if len(good):
    print("\nhorizons where the model forecasts and beats persistence:")
    for _, r in good.iterrows():
        print(f"  k={int(r.k)} {r.variant}: skill vs mean {r.skill_mean:+.3f}, "
              f"vs persistence {r.skill_vs_persistence_mean:+.3f}, "
              f"rmse {r.rmse_counts_mean:.1f} cases")
PY

echo ""
echo "results in $RESULTS"
echo "log $RUNLOG"
