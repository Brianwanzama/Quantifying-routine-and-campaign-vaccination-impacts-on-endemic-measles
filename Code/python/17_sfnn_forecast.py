#!/usr/bin/env python3
# =============================================================================
# 17_sfnn_forecast.py
#
# Trains on every year and forecasts one year past the end of the data.
#
# Different from the evaluation. The target biweek has not happened, so nothing
# exists to score against. What comes out is a forecast, not a measurement of
# accuracy.
#
# At k=26 the target is 26 biweeks past the forecast origin, which is one year.
# Forecasting from the last observed biweek, March 2020, puts the target in
# early 2021; forecasting from each biweek of the final year in turn gives a
# prediction for every biweek of 2020-21.
#
# Training on everything removes the validation split, so the epoch count comes
# from the tuned runs where it was chosen on a year resembling the forecast
# period. Holding a year back to preserve the stopping rule would mean
# forecasting with a model that never saw the most recent data.
#
# The model is refitted several times and the forecasts averaged. The spread is
# over initialisations, not a forecast interval, and should be described that
# way.
#
# One caution belongs with the output. The forecast period is April 2020
# onward, when surveillance was disrupted. The model cannot know that, so what
# it produces is what would have been expected under the reporting of previous
# years rather than what was observed.
#
# Output : forecast_k<K>.csv, forecast_summary.csv
# =============================================================================
import argparse
import os

import numpy as np
import pandas as pd
import torch
from torch import nn

import importlib
msf = importlib.import_module("10_sfnn_functions")

ROOT = "/.../.../.../Derived"
PREFIT_DEFAULT   = f"{ROOT}/sfnn/prefit"
TUNING_DEFAULT   = f"{ROOT}/sfnn/tuning"
FINAL_DEFAULT    = f"{ROOT}/sfnn/final"
FORECAST_DEFAULT = f"{ROOT}/sfnn/forecast"


def best_config(tuning_dir, variant, k):
    path = os.path.join(tuning_dir, variant, f"tuning_k{k}.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"no tuning results at {path}")
    d = pd.read_csv(path)
    if not len(d):
        raise ValueError(f"{path} is empty")
    key = "val_mse" if "val_mse" in d.columns else "test_mse"
    return d.sort_values(key).iloc[0].to_dict()


def tuned_epochs(final_dir, variant, k, fallback=60):
    """How long the tuned runs trained before early stopping ended them.

    Training on everything leaves nothing to stop on, so the schedule comes
    from the evaluation runs where it was chosen properly.
    """
    path = os.path.join(final_dir, f"summary_k{k}.csv")
    if not os.path.exists(path):
        print(f"  no summary at {path}, using {fallback} epochs")
        return fallback
    d = pd.read_csv(path)
    row = d[d["variant"] == variant]
    if not len(row) or "best_epoch_mean" not in d.columns:
        print(f"  no stopping epoch for {variant}, using {fallback}")
        return fallback
    return max(int(round(float(row["best_epoch_mean"].iloc[0]))), 5)


def main():
    ap = argparse.ArgumentParser(description="forecast beyond the data")
    ap.add_argument("--k", type=int, default=26, help="horizon in biweeks")
    ap.add_argument("--variant", default="no_nfhs",
                    choices=["with_nfhs", "no_nfhs"])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--n-origins", type=int, default=26,
                    help="how many of the final biweeks to forecast from; 26 "
                         "gives a prediction for every biweek of the year "
                         "ahead, 1 forecasts from the last observation only")
    ap.add_argument("--prefit-dir", default=PREFIT_DEFAULT)
    ap.add_argument("--tuning-dir", default=TUNING_DEFAULT)
    ap.add_argument("--final-dir", default=FINAL_DEFAULT)
    ap.add_argument("--out-dir", default=FORECAST_DEFAULT)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    k = args.k
    drop_group = None if args.variant == "with_nfhs" else "nfhs"

    print(f"k={k} | variant {args.variant} | device {device} | "
          f"{args.seeds} seeds")

    cfg = best_config(args.tuning_dir, args.variant, k)
    n_ep = tuned_epochs(args.final_dir, args.variant, k)
    t_lag = int(cfg["t_lag"])

    print(f"  lag {t_lag}, hidden {int(cfg['hidden_dim'])}, "
          f"layers {int(cfg['num_hidden_layers'])}, "
          f"lr {float(cfg['lr']):.5f}")
    print(f"  training for {n_ep} epochs, the count the tuned runs settled on")

    path = os.path.join(args.prefit_dir, f"k{k}_tlag{t_lag}.gzip")
    if not os.path.exists(path):
        raise FileNotFoundError(f"no matrix at {path}")
    X = pd.read_parquet(path)

    tf = pd.read_parquet(os.path.join(
        args.prefit_dir, f"k{k}_tlag{t_lag}_cases_transform.gzip"))
    c_mu, c_sd = float(tf["cases_mean"][0]), float(tf["cases_std"][0])

    # ---- everything becomes training ----------------------------------------
    # The split column is overwritten rather than filtered, so feature
    # selection goes through the same path as the evaluation runs.
    Xf = X.copy()
    Xf["split"] = "train"

        # the same selection process_data makes, so the forecast columns are
    # exactly those the model was fitted on
    feat_cols = list(Xf.filter(regex=msf.FEATURE_REGEX).columns)
    if drop_group is not None:
        g = msf.feature_groups(feat_cols)
        feat_cols = [c for c in feat_cols if c not in set(g[drop_group])]
    # ---- the forecast origins ------------------------------------------------
    # Every row carries lags from t-k backwards, all observed, and a target k
    # biweeks ahead. Taking the last n_origins biweeks as origins gives a
    # prediction for each of the n_origins biweeks after the data ends.
    last_t = int(Xf["t_global"].max())
    origins = sorted(Xf["t_global"].unique())[-args.n_origins:]
    fc = Xf[Xf["t_global"].isin(origins)].copy()

    print(f"  forecasting from t_global {min(origins)} to {last_t}, "
          f"{len(origins)} origins, {fc['uid'].nunique()} districts")
    print(f"  targets land at t_global {min(origins) + k} to {last_t + k}")

    miss = fc[feat_cols].isna().sum().sum()
    if miss:
        print(f"  {miss} missing feature values in the forecast rows")
        return

    # ---- fit and predict, several times --------------------------------------
    preds = []
    for seed in range(args.seeds):
        torch.manual_seed(seed)
        np.random.seed(seed)

        data, n_feat, _ = msf.process_data(Xf, drop_group=drop_group,
                                           verbose=(seed == 0))
        model = msf.NeuralNetwork(
            input_dim=n_feat,
            hidden_dim=int(cfg["hidden_dim"]),
            output_dim=1,
            num_hidden_layers=int(cfg["num_hidden_layers"])).to(device)
        opt = torch.optim.Adam(model.parameters(), lr=float(cfg["lr"]),
                               weight_decay=float(cfg["weight_decay"]))
        loaders = msf.get_dataloaders(data, batch_size=32)
        loss_fn = nn.MSELoss()

        for _ in range(n_ep):
            msf.train(model, device, loaders["train"], opt, loss_fn)

        model.eval()
        with torch.no_grad():
            xt = torch.tensor(fc[feat_cols].values).float().to(device)
            preds.append(model(xt).cpu().numpy().flatten())
        print(f"    seed {seed} done")

    P = np.vstack(preds)
    mean_trans, sd_trans = P.mean(axis=0), P.std(axis=0)

    # ---- assemble -------------------------------------------------------------
    out = pd.DataFrame({
        "state": fc["state"].values,
        "district": fc["district"].values,
        "uid": fc["uid"].values,
        "k": k,
        "variant": args.variant,
        "from_t_global": fc["t_global"].values,
        "target_t_global": fc["t_global"].values + k,
        "predicted_trans": mean_trans,
        "sd_across_seeds_trans": sd_trans,
        "predicted_cases": msf.back_transform(mean_trans, c_mu, c_sd),
    })
    if "persistence" in fc.columns:
        out["persistence_cases"] = msf.back_transform(
            fc["persistence"].to_numpy(), c_mu, c_sd)

    # the year and biweek the target falls in, for the caption
    out["target_year"] = 2008 + (out["target_t_global"] - 1) // 26
    out["target_biweek"] = ((out["target_t_global"] - 1) % 26) + 1

    f = os.path.join(args.out_dir, f"forecast_k{k}_{args.variant}.csv")
    out.to_csv(f, index=False)

    print(f"\n  wrote {os.path.basename(f)}, {len(out)} rows")
    print(f"  targets span {out['target_year'].min()} to "
          f"{out['target_year'].max()}")

    # ---- per district, summed over the forecast year --------------------------
    by_d = out.groupby(["state", "district", "uid"], as_index=False).agg(
        biweeks=("predicted_cases", "size"),
        predicted_total=("predicted_cases", "sum"),
        predicted_mean=("predicted_cases", "mean"))
    fs = os.path.join(args.out_dir, f"forecast_summary_k{k}_{args.variant}.csv")
    by_d.to_csv(fs, index=False)

    print("\n  by state:")
    print("  " + by_d.groupby("state").agg(
        districts=("district", "size"),
        total=("predicted_total", "sum"),
        median_per_district=("predicted_total", "median")
    ).round(1).to_string().replace("\n", "\n  "))

    print("\nThis is a forecast beyond the end of the data. There is no observed")
    print("value to compare it against, and the spread reported is over")
    print("initialisations rather than a forecast interval. The period forecast")
    print("is April 2020 onward, when surveillance was disrupted, so it")
    print("describes what the reporting of previous years would have implied.")


if __name__ == "__main__":
    main()
