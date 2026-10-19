#!/usr/bin/env python3
# =============================================================================
# 11_sfnn_raytune.py
#
# Hyperparameter search, one study per horizon and variant.
#
# The search runs twice per horizon, once with the static covariates and once
# without, each writing to its own folder. Tuning separately rather than tuning
# once and ablating afterwards is the fair comparison, since the model without
# the covariates then gets its own best architecture rather than inheriting one
# chosen while they were present.
#
# Trials are ranked on the validation split, so the test set plays no part in
# choosing either the architecture or the stopping epoch. The test error is
# recorded alongside for reference only.
#
# Search space
#   t_lag              {26, 52}, restricted to windows at least as long as k
#   hidden_dim         {32, 64, 128, 240}
#   num_hidden_layers  {1, 2}
#   lr                 log uniform, 1e-4 to 1e-2
#   weight_decay       uniform, 1e-4 to 1e-1
#
# The lag window is searched rather than fixed. Lags run from k to t_lag, so at
# k=26 a 26 window leaves a single lag and the feature blocks nearly vanish,
# while 52 gives 27. At the longer horizons the wider window is what makes them
# usable at all, so which wins is a real question.
#
# Depth is one or two layers: in earlier work every configuration that worked
# used a single layer, and three appeared only where the model failed.
#
# Usage
#   python 11_sfnn_raytune.py --k 2
#   python 11_sfnn_raytune.py --k 2 --drop-group nfhs
# =============================================================================
import argparse
import os

import pandas as pd
import torch
from torch import nn

import ray
from ray import tune
from ray.tune.schedulers import ASHAScheduler

import importlib
msf = importlib.import_module("10_sfnn_functions")

PREFIT_DEFAULT  = "/.../.../.../Derived/sfnn/prefit"
RESULTS_DEFAULT = "/.../.../.../Derived/sfnn/tuning"

T_LAG_OPTIONS = [26, 52]


def train_with_tuning(config):
    """One trial. Reports validation each epoch so the scheduler can stop it."""
    path = os.path.join(config["prefit"],
                        f"k{config['k']}_tlag{config['t_lag']}.gzip")
    cases = pd.read_parquet(path)

    data, num_features, _ = msf.process_data(
        cases, drop_group=config.get("drop_group"), verbose=False)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = msf.NeuralNetwork(
        input_dim=num_features,
        hidden_dim=config["hidden_dim"],
        output_dim=1,
        num_hidden_layers=config["num_hidden_layers"]).to(device)

    optimizer = torch.optim.Adam(model.parameters(),
                                 lr=config["lr"],
                                 weight_decay=config["weight_decay"])
    loaders = msf.get_dataloaders(data, batch_size=config["batch_size"])
    loss_fn = nn.MSELoss()

    has_val = "val" in loaders and len(loaders["val"].dataset) > 0

    for epoch in range(config["num_epochs"]):
        tr = msf.train(model, device, loaders["train"], optimizer, loss_fn)
        va = (msf.test(model, device, loaders["val"], loss_fn)
              if has_val else float("nan"))
        te = msf.test(model, device, loaders["test"], loss_fn)
        # the scheduler minimises val_mse; test_mse rides along so the gap can
        # be inspected later and is never used to rank or to stop
        tune.report({"val_mse": va, "test_mse": te,
                     "train_mse": tr, "epoch": epoch + 1})


def main():
    ap = argparse.ArgumentParser(description="SFNN hyperparameter search")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--num-samples", type=int, default=20,
                    help="draws of lr and weight decay per grid point")
    ap.add_argument("--max-num-epochs", type=int, default=150)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--drop-group", type=str, default=None,
                    help="ablation, e.g. nfhs, spatial, campaign")
    ap.add_argument("--cpus-per-trial", type=int, default=2)
    ap.add_argument("--gpus-per-trial", type=float, default=0.25)
    ap.add_argument("--max-cpus", type=int, default=32)
    ap.add_argument("--max-gpus", type=int, default=2)
    ap.add_argument("--prefit-dir", default=PREFIT_DEFAULT)
    ap.add_argument("--results-dir", default=RESULTS_DEFAULT)
    args = ap.parse_args()

    tag = "with_nfhs" if args.drop_group is None else f"no_{args.drop_group}"
    out_dir = os.path.join(args.results_dir, tag)
    os.makedirs(out_dir, exist_ok=True)

    print(f"horizon k={args.k} | variant {tag}")
    if args.drop_group:
        print(f"  the '{args.drop_group}' group is dropped before training")

    # A window shorter than the horizon would leave the most recent available
    # observation unused, so only windows at or above k are searched.
    opts = [t for t in T_LAG_OPTIONS if t >= args.k]
    if not opts:
        print(f"no window in {T_LAG_OPTIONS} reaches k={args.k}")
        return

    missing = []
    print("prefit files:")
    for t in opts:
        p = os.path.join(args.prefit_dir, f"k{args.k}_tlag{t}.gzip")
        ok = os.path.exists(p)
        print(f"  {os.path.basename(p)}  {'found' if ok else 'MISSING'}")
        if not ok:
            missing.append(p)
    if missing:
        print("\nrun the loader first for the missing windows")
        return

    if not ray.is_initialized():
        ray.init(num_cpus=args.max_cpus, num_gpus=args.max_gpus)

    config = {
        "k":                 args.k,
        "t_lag":             tune.grid_search(opts),
        "prefit":            args.prefit_dir,
        "drop_group":        args.drop_group,
        "num_epochs":        args.max_num_epochs,
        "batch_size":        args.batch_size,
        "hidden_dim":        tune.grid_search([32, 64, 128, 240]),
        "num_hidden_layers": tune.grid_search([1, 2]),
        "lr":                tune.loguniform(1e-4, 1e-2),
        "weight_decay":      tune.uniform(1e-4, 1e-1),
    }

    n_grid = len(opts) * 4 * 2
    print(f"\ngrid {n_grid} | samples {args.num_samples} | "
          f"trials {n_grid * args.num_samples} | epochs up to {args.max_num_epochs}")
    print(f"per trial: {args.cpus_per_trial} cpu, {args.gpus_per_trial} gpu")

    # stops weak trials early; the grace period lets every trial run five
    # epochs first, since a good configuration can start slowly
    scheduler = ASHAScheduler(metric="val_mse", mode="min",
                              max_t=args.max_num_epochs,
                              grace_period=5, reduction_factor=4)

    result = tune.run(
        train_with_tuning,
        resources_per_trial={"cpu": args.cpus_per_trial,
                             "gpu": args.gpus_per_trial},
        config=config,
        num_samples=args.num_samples,
        scheduler=scheduler,
        verbose=1,
        raise_on_failed_trial=False)

    rows = []
    for trial in result.trials:
        if trial.last_result and "val_mse" in trial.last_result:
            rows.append({
                "k":                 args.k,
                "variant":           tag,
                "trial_id":          trial.trial_id,
                "val_mse":           trial.last_result["val_mse"],
                "test_mse":          trial.last_result.get("test_mse"),
                "train_mse":         trial.last_result.get("train_mse"),
                "epoch":             trial.last_result.get("epoch"),
                "t_lag":             trial.config.get("t_lag"),
                "hidden_dim":        trial.config.get("hidden_dim"),
                "num_hidden_layers": trial.config.get("num_hidden_layers"),
                "lr":                trial.config.get("lr"),
                "weight_decay":      trial.config.get("weight_decay"),
            })
    if not rows:
        raise RuntimeError("no trial completed, see ~/ray_results for the logs")

    # ranked on validation, since ranking on test would select using the set
    # the error is meant to be held out from
    df = pd.DataFrame(rows).sort_values("val_mse").reset_index(drop=True)
    path = os.path.join(out_dir, f"tuning_k{args.k}.csv")
    df.to_csv(path, index=False)

    best = df.iloc[0]
    print(f"\n{'=' * 55}")
    print(f"best for k={args.k}, {tag}")
    print(f"{'=' * 55}")
    print(f"  validation mse      {best['val_mse']:.6f}   (what was ranked)")
    print(f"  test mse            {best['test_mse']:.6f}   (held out)")
    print(f"  train mse           {best['train_mse']:.6f}")
    print(f"  hidden dimension    {int(best['hidden_dim'])}")
    print(f"  hidden layers       {int(best['num_hidden_layers'])}")
    print(f"  learning rate       {best['lr']:.6f}")
    print(f"  weight decay        {best['weight_decay']:.6f}")
    print(f"  lag window          {int(best['t_lag'])}")
    print(f"  epochs completed    {int(best['epoch'])}")
    print(f"  train minus validation {best['train_mse'] - best['val_mse']:+.6f}")

    if best["epoch"] > 0.9 * args.max_num_epochs:
        print("\n  stopping near the cap, so raising it may improve the model")

    # which window won, and by how much, since that is a question in itself
    by_lag = df.groupby("t_lag")["val_mse"].min()
    print("\nbest validation mse by lag window:")
    for t, v in by_lag.items():
        print(f"  {int(t)} biweeks   {v:.6f}")

    print(f"\n{len(df)} trials, written to {path}")


if __name__ == "__main__":
    main()
