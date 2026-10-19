#!/usr/bin/env bash
# =============================================================================
# 13_run_explain.sh
#
# Grouped Shapley attribution at each horizon, both variants.
#
# Needs the models the training script saved, so run 12 first.
#
# Attribution is over a subsample of test rows. Shapley sampling costs roughly
# groups times samples forward passes per row, so the full test set of 7,505
# would be slow and a few thousand rows give a stable mean.
#
# Usage
#   bash 13_run_explain.sh
# =============================================================================
set -euo pipefail

ROOT="/.../.../.../"
CODE="$ROOT/Code/python"
PREFIT="$ROOT/Derived/sfnn/prefit"
FINAL="$ROOT/Derived/sfnn/final"
EXPLAIN="$ROOT/Derived/sfnn/explain"
LOGS="$ROOT/Derived/sfnn/logs"

# Every horizon that was trained. The working ones, k=4 to k=12, show what a
# forecasting model attends to; k=1, k=20 and k=26 lose to persistence, and
# what a failing model leans on is worth seeing rather than assuming.
KS=(1 4 6 12 20 26)

N_SAMPLES=25
N_ROWS=3000        # 0 for every test row
BATCH_ROWS=500

mkdir -p "$EXPLAIN" "$LOGS"
cd "$CODE"

echo "checking the models"
missing=0
for k in "${KS[@]}"; do
  for v in with_nfhs no_nfhs; do
    f="$FINAL/models/k${k}_${v}.pt"
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
  echo "run 12_run_training.sh first; it saves the models"
  exit 1
fi

STAMP=$(date +%Y%m%d_%H%M%S)
{
  for k in "${KS[@]}"; do
    echo ""
    echo "##############################################################"
    echo "# horizon k=$k"
    echo "##############################################################"
    python 13_sfnn_explain.py \
      --k "$k" \
      --n-samples "$N_SAMPLES" \
      --n-rows "$N_ROWS" \
      --batch-rows "$BATCH_ROWS" \
      --prefit-dir "$PREFIT" \
      --final-dir "$FINAL" \
      --out-dir "$EXPLAIN"
  done
} 2>&1 | tee "$LOGS/explain_${STAMP}.log"

echo ""
echo "=============================================================="
echo "importance across horizons"
echo "=============================================================="
python - <<'PY' 2>&1 | tee -a "$LOGS/explain_${STAMP}.log"
import glob, os
import pandas as pd

EXPLAIN = "/home/brain/Msc_project/Manuscript/Derived/sfnn/explain"
files = sorted(glob.glob(os.path.join(EXPLAIN, "shap_importance_k*.csv")))
if not files:
    print("nothing found")
    raise SystemExit

d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

for v in sorted(d["variant"].unique()):
    sub = d[d["variant"] == v]
    w = sub.pivot_table(index="group", columns="k", values="share")
    w = w.loc[w.mean(axis=1).sort_values(ascending=False).index]
    print(f"\nshare of attribution, {v}, by horizon:")
    print(w.round(1).to_string())

# whether the covariate pathways carry anything, at every horizon
w = d[d["variant"] == "with_nfhs"]
if len(w):
    print("\nthe four covariate pathways together, by horizon:")
    for k in sorted(w["k"].unique()):
        sub = w[w.k == k]
        nf = sub[sub["group"].str.startswith("nfhs_")]["share"].sum()
        print(f"  k={k:<3d} {nf:.1f} percent")
    print("\nread against the ablation: if these are small and the ablation")
    print("showed no benefit, the two agree.")
PY

echo ""
echo "results in $EXPLAIN"
