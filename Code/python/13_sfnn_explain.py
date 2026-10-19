#!/usr/bin/env python3
# =============================================================================
# 13_sfnn_explain.py
#
# Grouped Shapley attribution for the trained SFNN.
#
# Every feature is assigned to a pathway, and Shapley values are computed over
# those groups rather than over individual columns. A single lag of a district's
# own incidence means little on its own; the twenty six of them together mean
# something. Grouping also handles the correlation within blocks, which would
# otherwise split one signal across many small attributions.
#
# The four covariate pathways are kept separate rather than pooled into one
# nfhs block. The ablation showed the block as a whole does not help, and the
# question this answers is whether any single pathway carries signal anyway.
#
# Run on both variants. The with_nfhs model is the one that has covariates to
# attribute; the no_nfhs model is a check that the other groups behave the same
# way without them.
#
# Attribution is over a subsample of the test rows by default. Shapley sampling
# costs roughly groups times samples forward passes per row, so the full test
# set is expensive and a few thousand rows give a stable mean.
#
# Usage
#   python 13_sfnn_explain.py --k 6
#   python 13_sfnn_explain.py --k 6 --n-rows 0        # every test row
# =============================================================================
import argparse
import os

import numpy as np
import pandas as pd
import torch
from captum.attr import ShapleyValueSampling

import importlib
msf = importlib.import_module("10_sfnn_functions")

PREFIT_DEFAULT  = "/home/brain/Msc_project/Manuscript/Derived/sfnn/prefit"
FINAL_DEFAULT   = "/home/brain/Msc_project/Manuscript/Derived/sfnn/final"
EXPLAIN_DEFAULT = "/home/brain/Msc_project/Manuscript/Derived/sfnn/explain"

VARIANTS = [("with_nfhs", None), ("no_nfhs", "nfhs")]


def assign_groups(cols):
    """Map every column to a pathway. Returns the names and an index per column.

    The four covariate pathways stay separate, and the two spatial sets are
    split, since a neighbour's incidence and a large district's incidence are
    different mechanisms.
    """
    order = ["own", "susc", "core", "nearest", "births", "pop", "vacc",
             "campaign", "state",
             "nfhs_supply", "nfhs_contact", "nfhs_access", "nfhs_nutrition"]
    idx = {g: i for i, g in enumerate(order)}

    out, unmatched = {}, []
    for c in cols:
        if c.startswith("cases_lag_"):
            g = "own"
        elif c.startswith("susc_lag_"):
            g = "susc"
        elif c.startswith("core"):
            g = "core"
        elif c.startswith("nbc"):
            g = "nearest"
        elif c.startswith("births_lag_"):
            g = "births"
        elif c.startswith("pop_lag_"):
            g = "pop"
        elif c.startswith(("v1_lag", "v2_lag")):
            g = "vacc"
        elif c.startswith("mr_"):
            g = "campaign"
        elif c.startswith("st_"):
            g = "state"
        elif c.startswith("nfhs_"):
            # the loader put the pathway in the name, as nfhs_supply_measles_zero
            parts = c.split("_")
            g = f"nfhs_{parts[1]}" if len(parts) > 2 and f"nfhs_{parts[1]}" in idx \
                else "nfhs_supply"
        else:
            unmatched.append(c)
            g = "own"
        out[c] = idx[g]

    if unmatched:
        raise ValueError(f"{len(unmatched)} columns match no group: "
                         f"{unmatched[:8]}")

    present = sorted(set(out.values()))
    names = [order[i] for i in present]
    remap = {old: new for new, old in enumerate(present)}
    return names, np.array([remap[out[c]] for c in cols])


def main():
    ap = argparse.ArgumentParser(description="grouped Shapley for the SFNN")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--n-samples", type=int, default=25,
                    help="Shapley samples per group; higher is steadier, slower")
    ap.add_argument("--n-rows", type=int, default=3000,
                    help="test rows to attribute over; 0 for all")
    ap.add_argument("--batch-rows", type=int, default=500,
                    help="rows per pass, to bound memory")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--prefit-dir", default=PREFIT_DEFAULT)
    ap.add_argument("--final-dir", default=FINAL_DEFAULT)
    ap.add_argument("--out-dir", default=EXPLAIN_DEFAULT)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    print(f"horizon k={args.k} | device {device}")

    summ_path = os.path.join(args.final_dir, f"summary_k{args.k}.csv")
    if not os.path.exists(summ_path):
        raise FileNotFoundError(f"no summary at {summ_path}, train first")
    summ = pd.read_csv(summ_path)

    all_long, all_wide = [], []

    for variant, drop_group in VARIANTS:
        row = summ[summ["variant"] == variant]
        if not len(row):
            print(f"  no {variant} in the summary, skipping")
            continue
        row = row.iloc[0]

        t_lag = int(row["t_lag"])
        path = os.path.join(args.prefit_dir, f"k{args.k}_tlag{t_lag}.gzip")
        cases = pd.read_parquet(path)

        data, n_feat, ids = msf.process_data(cases, drop_group=drop_group,
                                             verbose=False)

        # the column order has to match what the model was fitted on, so it is
        # taken the same way rather than reconstructed
        Xte = cases[cases["split"] == "test"].filter(regex=msf.FEATURE_REGEX)
        if drop_group is not None:
            g = msf.feature_groups(Xte.columns)
            Xte = Xte.drop(columns=g[drop_group])
        Xte = Xte.replace([np.inf, -np.inf], 0).fillna(0)
        assert Xte.shape[1] == n_feat, \
            f"columns {Xte.shape[1]} against the model's {n_feat}"

        model_path = os.path.join(args.final_dir, "models",
                                  f"k{args.k}_{variant}.pt")
        if not os.path.exists(model_path):
            print(f"  no model at {model_path}, skipping {variant}")
            continue

        model = msf.NeuralNetwork(
            input_dim=n_feat,
            hidden_dim=int(row["hidden_dim"]),
            output_dim=1,
            num_hidden_layers=int(row["num_hidden_layers"])).to(device)
        model.load_state_dict(torch.load(model_path, map_location=device))
        model.eval()

        names, garr = assign_groups(list(Xte.columns))
        print(f"\n=== {variant} | lag {t_lag} | {n_feat} features "
              f"in {len(names)} groups ===")
        for i, g in enumerate(names):
            print(f"  {g:<16s} {int((garr == i).sum())}")

        # subsample, since Shapley costs groups times samples passes per row
        idt = ids["test"].reset_index(drop=True)
        n_all = len(Xte)
        if args.n_rows and args.n_rows < n_all:
            take = np.sort(rng.choice(n_all, args.n_rows, replace=False))
            print(f"  attributing over {len(take)} of {n_all} test rows")
        else:
            take = np.arange(n_all)
            print(f"  attributing over all {n_all} test rows")

        Xs = Xte.iloc[take]
        ids_s = idt.iloc[take].reset_index(drop=True)

        mask = torch.tensor(garr).to(device)
        svs = ShapleyValueSampling(model)

        parts = []
        for s in range(0, len(Xs), args.batch_rows):
            chunk = torch.from_numpy(
                Xs.iloc[s:s + args.batch_rows].to_numpy()).float().to(device)
            a = svs.attribute(chunk, feature_mask=mask,
                              n_samples=args.n_samples)
            # captum returns the attribution shaped like the input, but the
            # network's output of (n, 1) can leave it with a trailing axis, so
            # it is flattened back to rows by features
            a = a.detach().cpu().numpy().reshape(chunk.shape[0], -1)
            if a.shape[1] != n_feat:
                raise RuntimeError(
                    f"attribution has {a.shape[1]} columns against "
                    f"{n_feat} features; chunk was {tuple(chunk.shape)}")
            parts.append(a)
            print(f"    {min(s + args.batch_rows, len(Xs))} of {len(Xs)}")
        attr = np.vstack(parts)

        # every column in a group carries that group's value, so one is taken
        first = {g: int(np.where(garr == g)[0][0]) for g in range(len(names))}
        wide = pd.DataFrame({names[g]: attr[:, first[g]]
                             for g in range(len(names))})

        # The identifiers are prefixed because one of the groups is also called
        # state, and a bare column of that name would overwrite its attribution
        # and then be handed to abs() as text.
        for c in ("state", "district", "uid", "t_global", "time"):
            if c in ids_s.columns:
                wide[f"id_{c}"] = ids_s[c].values
        wide["k"] = args.k
        wide["variant"] = variant
        all_wide.append(wide)

        # ---- what it says ----------------------------------------------------
        # The mean absolute attribution is how much a group moves the
        # prediction, regardless of direction. The signed mean says which way.
        imp = pd.DataFrame({
            "group": names,
            "mean_abs": [np.abs(attr[:, first[g]]).mean() for g in range(len(names))],
            "mean_signed": [attr[:, first[g]].mean() for g in range(len(names))],
            "n_features": [int((garr == g).sum()) for g in range(len(names))]})
        imp["share"] = 100 * imp["mean_abs"] / imp["mean_abs"].sum()
        imp = imp.sort_values("mean_abs", ascending=False).reset_index(drop=True)
        imp["k"] = args.k
        imp["variant"] = variant
        all_long.append(imp)

        print(f"\n  importance, {variant}:")
        print("  " + imp[["group", "mean_abs", "mean_signed", "share",
                          "n_features"]].round(4).to_string(index=False)
              .replace("\n", "\n  "))

        # per feature attribution, so a large group is not credited simply for
        # being large
        per = imp.assign(per_feature=imp["mean_abs"] / imp["n_features"])
        print("\n  per feature, which corrects for group size:")
        print("  " + per.sort_values("per_feature", ascending=False)[
            ["group", "per_feature"]].round(5).to_string(index=False)
            .replace("\n", "\n  "))

        # by state, since the two may lean on different things
        if "id_state" in wide.columns:
            bys = wide.groupby("id_state")[names].apply(
                lambda d: d.abs().mean()).round(4)
            print("\n  mean absolute attribution by state:")
            print("  " + bys.to_string().replace("\n", "\n  "))

    if not all_long:
        print("\nnothing attributed")
        return

    imp_all = pd.concat(all_long, ignore_index=True)
    imp_all.to_csv(os.path.join(args.out_dir, f"shap_importance_k{args.k}.csv"),
                   index=False)
    pd.concat(all_wide, ignore_index=True).to_csv(
        os.path.join(args.out_dir, f"shap_attributions_k{args.k}.csv"),
        index=False)

    # ---- do the covariate pathways carry anything? ---------------------------
    # The ablation said the block as a whole does not help. If a single pathway
    # still shows a substantial attribution, that is worth knowing; if all four
    # sit at the bottom, the two findings agree.
    w = imp_all[imp_all["variant"] == "with_nfhs"]
    if len(w):
        nf = w[w["group"].str.startswith("nfhs_")]
        rest = w[~w["group"].str.startswith("nfhs_")]
        print(f"\n{'=' * 55}")
        print("the covariate pathways against everything else")
        print(f"{'=' * 55}")
        print(f"  the four pathways together: "
              f"{nf['share'].sum():.1f} percent of attribution")
        print(f"  the largest single pathway: {nf.iloc[0]['group']} at "
              f"{nf.iloc[0]['share']:.1f} percent" if len(nf) else "")
        print(f"  the smallest other group:   {rest.iloc[-1]['group']} at "
              f"{rest.iloc[-1]['share']:.1f} percent" if len(rest) else "")
        if len(nf) and len(rest) and nf["share"].max() < rest["share"].min():
            print("\n  every covariate pathway attributes less than every other")
            print("  group, which agrees with the ablation")

    print(f"\nwritten to {args.out_dir}")


if __name__ == "__main__":
    main()
