#!/usr/bin/env python3
# =============================================================================
# 12_sfnn_train.py
#
# Trains the final model at one horizon using the configuration the search
# settled on, for both variants, and writes the predictions.
#
# Both variants are trained in one run, so they are fitted on the same data
# with the same seeds and differ only in whether the covariates were present.
# Each takes its own tuned architecture, since the search was run separately.
#
# Training stops at the epoch that fits the validation year best and keeps that
# state, rather than running a fixed budget and keeping whatever it ended in.
#
# Two skill scores are reported. Against always predicting the training mean,
# and against persistence, meaning carrying the last observed value forward.
# The second matters more: in a series where most biweeks repeat their
# predecessor, persistence is what a person would do without a model.
#
# Several seeds, because a single fit gives one number and a small difference
# between variants can be initialisation noise. The spread is what tells them
# apart.
#
# Usage
#   python 12_sfnn_train.py --k 4 --seeds 3
# =============================================================================
import argparse
import os

import numpy as np
import pandas as pd
import torch
from torch import nn

import importlib
msf = importlib.import_module("10_sfnn_functions")

PREFIT_DEFAULT  = "/.../.../.../Derived/sfnn/prefit"
TUNING_DEFAULT  = "/.../.../.../Derived/sfnn/tuning"
RESULTS_DEFAULT = "/.../.../.../Derived/sfnn/final"

# the tag each variant writes under, and the group it drops
VARIANTS = [("with_nfhs", None), ("no_nfhs", "nfhs")]


def best_config(tuning_dir, variant, k):
    """The lowest validation error configuration from that study.

    Ranking on test would select the architecture using the set the reported
    error is meant to be held out from.
    """
    path = os.path.join(tuning_dir, variant, f"tuning_k{k}.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"no tuning results at {path}")
    d = pd.read_csv(path)
    if not len(d):
        raise ValueError(f"{path} is empty")
    key = "val_mse" if "val_mse" in d.columns else "test_mse"
    if key == "test_mse":
        print(f"  {path} predates the validation split, ranking on test")
    return d.sort_values(key).iloc[0].to_dict()


def fit_once(cases, cfg, drop_group, seed, device, max_epochs, patience):
    torch.manual_seed(seed)
    np.random.seed(seed)

    data, n_feat, ids = msf.process_data(
        cases, drop_group=drop_group, verbose=(seed == 0))

    model = msf.NeuralNetwork(
        input_dim=n_feat,
        hidden_dim=int(cfg["hidden_dim"]),
        output_dim=1,
        num_hidden_layers=int(cfg["num_hidden_layers"])).to(device)

    opt = torch.optim.Adam(model.parameters(), lr=float(cfg["lr"]),
                           weight_decay=float(cfg["weight_decay"]))
    loaders = msf.get_dataloaders(data, batch_size=32)

    model, history, best_epoch, best_val = msf.fit_with_early_stopping(
        model, device, loaders, opt, max_epochs=max_epochs,
        patience=patience, verbose=(seed == 0))

    model.eval()
    preds = {}
    with torch.no_grad():
        for nm, d in data.items():
            if len(d) == 0:
                continue
            preds[nm] = (ids[nm], d.y.numpy().flatten(),
                         model(d.X.to(device)).cpu().numpy().flatten())

    return model, n_feat, history, best_epoch, best_val, preds


def main():
    ap = argparse.ArgumentParser(description="SFNN final training")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--epochs", type=int, default=150,
                    help="the cap; early stopping usually ends sooner")
    ap.add_argument("--patience", type=int, default=15)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--prefit-dir", default=PREFIT_DEFAULT)
    ap.add_argument("--tuning-dir", default=TUNING_DEFAULT)
    ap.add_argument("--results-dir", default=RESULTS_DEFAULT)
    args = ap.parse_args()

    os.makedirs(args.results_dir, exist_ok=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"horizon k={args.k} | device {device} | {args.seeds} seeds\n")

    all_pred, all_summary, all_hist, all_dist = [], [], [], []

    for variant, drop_group in VARIANTS:
        cfg = best_config(args.tuning_dir, variant, args.k)

        # the window is whichever the search chose, so the two variants may
        # differ in it and each reads its own matrix
        t_lag = int(cfg["t_lag"])
        path = os.path.join(args.prefit_dir, f"k{args.k}_tlag{t_lag}.gzip")
        if not os.path.exists(path):
            raise FileNotFoundError(f"no matrix at {path}")
        cases = pd.read_parquet(path)

        tf = pd.read_parquet(os.path.join(
            args.prefit_dir, f"k{args.k}_tlag{t_lag}_cases_transform.gzip"))
        c_mu, c_sd = float(tf["cases_mean"][0]), float(tf["cases_std"][0])

        print("=" * 58)
        print(f"{variant}  |  lag {t_lag}, hidden {int(cfg['hidden_dim'])}, "
              f"layers {int(cfg['num_hidden_layers'])}, "
              f"lr {float(cfg['lr']):.5f}")
        print("=" * 58)

        seed_rows = []
        for seed in range(args.seeds):
            model, n_feat, hist, best_epoch, best_val, preds = fit_once(
                cases, cfg, drop_group, seed, device, args.epochs, args.patience)

            id_te, y_te, p_te = preds["test"]
            id_tr, y_tr, p_tr = preds["train"]

            mse_te = float(np.mean((p_te - y_te) ** 2))
            mse_tr = float(np.mean((p_tr - y_tr) ** 2))

            obs = msf.back_transform(y_te, c_mu, c_sd)
            prd = msf.back_transform(p_te, c_mu, c_sd)
            rmse_counts = float(np.sqrt(np.mean((prd - obs) ** 2)))
            mae_counts = float(np.mean(np.abs(prd - obs)))

            # Has it collapsed towards a constant? A prediction varying far
            # less than the data can post a respectable error while saying
            # nothing about when cases rise or fall.
            sd_ratio = float(np.std(p_te) / (np.std(y_te) or 1.0))
            corr = (float(np.corrcoef(p_te, y_te)[0, 1])
                    if np.std(p_te) > 0 and np.std(y_te) > 0 else np.nan)

            naive = float(np.mean((y_te - np.mean(y_tr)) ** 2))
            skill = 1.0 - mse_te / naive if naive > 0 else np.nan

            if "persistence" in id_te.columns:
                pers = id_te["persistence"].to_numpy()
                ok = np.isfinite(pers)
                pers_mse = float(np.mean((pers[ok] - y_te[ok]) ** 2))
                skill_pers = (1.0 - float(np.mean((p_te[ok] - y_te[ok]) ** 2))
                              / pers_mse) if pers_mse > 0 else np.nan
            else:
                pers_mse, skill_pers = np.nan, np.nan

            seed_rows.append(dict(
                seed=seed, test_mse=mse_te, train_mse=mse_tr,
                val_mse=best_val, best_epoch=best_epoch,
                rmse_counts=rmse_counts, mae_counts=mae_counts,
                sd_ratio=sd_ratio, correlation=corr,
                naive_mse=naive, skill=skill,
                persistence_mse=pers_mse, skill_vs_persistence=skill_pers))

            print(f"  seed {seed}: stopped {best_epoch:3d} | val {best_val:.4f} "
                  f"| test {mse_te:.4f} | rmse {rmse_counts:6.1f} "
                  f"| corr {corr:+.3f} | spread {sd_ratio:.3f} "
                  f"| skill vs mean {skill:+.3f} vs persistence {skill_pers:+.3f}")

            hist["seed"] = seed
            hist["variant"] = variant
            hist["k"] = args.k
            all_hist.append(hist)

            # the first seed's model is kept, since the explain script needs
            # a fitted model and refitting there would give a different one
            if seed == 0:
                mdir = os.path.join(args.results_dir, "models")
                os.makedirs(mdir, exist_ok=True)
                torch.save(model.state_dict(),
                           os.path.join(mdir, f"k{args.k}_{variant}.pt"))

                pr = id_te.copy()
                pr["variant"] = variant
                pr["k"] = args.k
                pr["observed_trans"] = y_te
                pr["predicted_trans"] = p_te
                pr["observed"] = obs
                pr["predicted"] = prd
                if "persistence" in pr.columns:
                    pr["persistence_counts"] = msf.back_transform(
                        pr["persistence"].to_numpy(), c_mu, c_sd)
                all_pred.append(pr)

                d = pr.groupby(["state", "district"]).apply(
                    lambda g: pd.Series({
                        "rmse": np.sqrt(np.mean((g["predicted"] - g["observed"]) ** 2)),
                        "mae": np.mean(np.abs(g["predicted"] - g["observed"])),
                        "mean_observed": g["observed"].mean(),
                        "correlation": (np.corrcoef(g["predicted"], g["observed"])[0, 1]
                                        if g["predicted"].std() > 0 and
                                           g["observed"].std() > 0 else np.nan),
                    }), include_groups=False).reset_index()
                d["variant"] = variant
                d["k"] = args.k
                all_dist.append(d)

        s = pd.DataFrame(seed_rows)
        summ = dict(
            k=args.k, variant=variant, n_features=n_feat, t_lag=t_lag,
            drop_group=("" if drop_group is None else drop_group),
            hidden_dim=int(cfg["hidden_dim"]),
            num_hidden_layers=int(cfg["num_hidden_layers"]),
            lr=float(cfg["lr"]), weight_decay=float(cfg["weight_decay"]),
            test_mse_mean=s["test_mse"].mean(), test_mse_sd=s["test_mse"].std(),
            train_mse_mean=s["train_mse"].mean(),
            val_mse_mean=s["val_mse"].mean(),
            best_epoch_mean=s["best_epoch"].mean(),
            rmse_counts_mean=s["rmse_counts"].mean(),
            mae_counts_mean=s["mae_counts"].mean(),
            correlation_mean=s["correlation"].mean(),
            sd_ratio_mean=s["sd_ratio"].mean(),
            skill_mean=s["skill"].mean(),
            persistence_mse_mean=s["persistence_mse"].mean(),
            skill_vs_persistence_mean=s["skill_vs_persistence"].mean())
        all_summary.append(summ)

        print(f"\n  across {args.seeds} seeds: test mse "
              f"{summ['test_mse_mean']:.4f} (sd {summ['test_mse_sd']:.4f}), "
              f"stopping at epoch {summ['best_epoch_mean']:.0f} on average")
        if summ["sd_ratio_mean"] < 0.3:
            print("  the predictions vary far less than the data, so this is")
            print("  close to a constant and the error is not informative")
        if summ["skill_mean"] < 0:
            print("  worse than always predicting the training mean")
        sp = summ["skill_vs_persistence_mean"]
        if np.isfinite(sp) and sp < 0:
            print("  and worse than carrying the last value forward")
        print()

    # ---- the comparison ------------------------------------------------------
    su = pd.DataFrame(all_summary)
    print("=" * 58)
    print(f"with and without the static covariates, k={args.k}")
    print("=" * 58)
    cols = ["variant", "t_lag", "n_features", "val_mse_mean", "test_mse_mean",
            "test_mse_sd", "rmse_counts_mean", "correlation_mean",
            "sd_ratio_mean", "skill_mean", "skill_vs_persistence_mean"]
    print(su[cols].round(4).to_string(index=False))

    w = su[su["variant"] == "with_nfhs"].iloc[0]
    n = su[su["variant"] == "no_nfhs"].iloc[0]
    diff = 100 * (n["test_mse_mean"] - w["test_mse_mean"]) / n["test_mse_mean"]
    print(f"\nthe covariates change the test error by {diff:+.2f} percent")
    if abs(w["test_mse_mean"] - n["test_mse_mean"]) < 2 * max(
            w["test_mse_sd"], n["test_mse_sd"]):
        print("that is within the variation between seeds, so they cannot be")
        print("said to help or hinder at this horizon")

    # ---- write ---------------------------------------------------------------
    tag = f"k{args.k}"
    pd.concat(all_pred).to_csv(
        os.path.join(args.results_dir, f"predictions_{tag}.csv"), index=False)
    pd.concat(all_dist).to_csv(
        os.path.join(args.results_dir, f"district_errors_{tag}.csv"), index=False)
    pd.concat(all_hist).to_csv(
        os.path.join(args.results_dir, f"training_history_{tag}.csv"), index=False)
    su.to_csv(os.path.join(args.results_dir, f"summary_{tag}.csv"), index=False)

    print(f"\nwritten to {args.results_dir}")


if __name__ == "__main__":
    main()
