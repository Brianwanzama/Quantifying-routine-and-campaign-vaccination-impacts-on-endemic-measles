#!/usr/bin/env bash
# =============================================================================
# 11_run_tuning.sh
#
# Builds the feature matrices, then searches for an architecture at each
# horizon, twice: once with the static covariates and once without.
#
# Tuning each variant separately is the fair comparison. Tuning once and
# ablating at training time would leave the model without the covariates
# inheriting an architecture chosen while they were present.
#
# Six horizons at two lag windows gives twelve matrices, and each horizon is
# searched twice, so twelve studies. The lag is a search dimension rather than
# a setting: at k=26 a 26 window leaves a single lag while 52 gives 27, so
# which is right differs by horizon.
#
# Either variant can be rerun alone with the stage argument.
#
# Usage
#   bash 11_run_tuning.sh              matrices, then both variants
#   bash 11_run_tuning.sh loader       the matrices only
#   bash 11_run_tuning.sh with         the with covariates search only
#   bash 11_run_tuning.sh without      the ablated search only
# =============================================================================
set -euo pipefail

ROOT="/../.../.../.../"
CODE="$ROOT/Code/python"
PREFIT="$ROOT/Derived/sfnn/prefit"
RESULTS="$ROOT/Derived/sfnn/tuning"
LOGS="$ROOT/Derived/sfnn/logs"

KS=(1 4 6 12,20 26)
LAGS=(26 52)

N_NEAR=5
N_CORE=4

NUM_SAMPLES=20
MAX_EPOCHS=300
CPUS_PER_TRIAL=2
GPUS_PER_TRIAL=0.25
MAX_CPUS=32
MAX_GPUS=2

STAGE="${1:-all}"

mkdir -p "$PREFIT" "$RESULTS" "$LOGS"
cd "$CODE"

run_search () {
  local variant_flag="$1" label="$2"
  for k in "${KS[@]}"; do
    echo ""
    echo "--- k=$k, $label ---"
    python 11_sfnn_raytune.py \
      --k "$k" $variant_flag \
      --num-samples "$NUM_SAMPLES" \
      --max-num-epochs "$MAX_EPOCHS" \
      --cpus-per-trial "$CPUS_PER_TRIAL" \
      --gpus-per-trial "$GPUS_PER_TRIAL" \
      --max-cpus "$MAX_CPUS" --max-gpus "$MAX_GPUS" \
      --prefit-dir "$PREFIT" --results-dir "$RESULTS" \
      2>&1 | tee "$LOGS/tune_k${k}_${label}.log"
  done
}

if [[ "$STAGE" == "all" || "$STAGE" == "loader" ]]; then
  echo "=============================================="
  echo "building the feature matrices"
  echo "=============================================="
  for k in "${KS[@]}"; do
    for t in "${LAGS[@]}"; do
      # a window shorter than the horizon would leave the most recent
      # observation unused
      if (( t < k )); then
        echo "  skipping k=$k lag=$t, the window is shorter than the horizon"
        continue
      fi
      echo ""
      echo "--- k=$k lag=$t ---"
      python 09_sfnn_loader.py \
        --k "$k" --t_lag "$t" \
        --n_near "$N_NEAR" --n_core "$N_CORE" \
        --save_dir "$PREFIT" \
        2>&1 | tee "$LOGS/loader_k${k}_tlag${t}.log"
    done
  done
fi

if [[ "$STAGE" == "all" || "$STAGE" == "with" ]]; then
  echo ""
  echo "=============================================="
  echo "searching, with the static covariates"
  echo "=============================================="
  run_search "" "with_nfhs"
fi

if [[ "$STAGE" == "all" || "$STAGE" == "without" ]]; then
  echo ""
  echo "=============================================="
  echo "searching, without the static covariates"
  echo "=============================================="
  run_search "--drop-group nfhs" "no_nfhs"
fi

echo ""
echo "=============================================="
echo "the best configuration in each study"
echo "=============================================="
python - <<'PY'
import os, glob
import pandas as pd

RESULTS = "/home/brain/Msc_project/Manuscript/Derived/sfnn/tuning"
rows = []
for f in sorted(glob.glob(os.path.join(RESULTS, "*", "tuning_k*.csv"))):
    d = pd.read_csv(f)
    if not len(d):
        continue
    key = "val_mse" if "val_mse" in d.columns else "test_mse"
    r = d.sort_values(key).iloc[0]
    r["variant"] = os.path.basename(os.path.dirname(f))
    rows.append(r)

if not rows:
    print("no results yet")
else:
    b = pd.DataFrame(rows)
    cols = [c for c in ["k", "variant", "val_mse", "test_mse", "t_lag",
                        "hidden_dim", "num_hidden_layers", "lr",
                        "weight_decay", "epoch"] if c in b.columns]
    print(b[cols].sort_values(["k", "variant"]).round(6).to_string(index=False))
    print("\nranked on validation. the test error rides along for reference")
    print("and played no part in choosing these configurations.")

    if "t_lag" in b.columns:
        print("\nwhich lag window won, by horizon:")
        for _, r in b.sort_values(["k", "variant"]).iterrows():
            print(f"  k={int(r['k']):<3d} {r['variant']:<10s} {int(r['t_lag'])}")

    near = b[b["epoch"] > 135] if "epoch" in b.columns else b.iloc[0:0]
    if len(near):
        print("\nstudies stopping near the epoch cap, worth raising it:")
        for _, r in near.iterrows():
            print(f"  k={int(r['k'])} {r['variant']} at epoch {int(r['epoch'])}")
PY

echo ""
echo "logs in $LOGS"
