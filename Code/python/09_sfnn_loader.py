#!/usr/bin/env python3
# =============================================================================
# 09_sfnn_loader.py
#
# Per horizon feature matrix, built from sfnn_master.csv.
#
# Follows Madden et al. (2024), adapted to two states of Indian districts.
# What differs from the multi-state version this is descended from:
#
#   two states rather than ten, so the state block is one indicator
#   no IMI, so that group is gone
#   thirteen static covariates rather than five, grouped by pathway
#   neighbours found within a state, not pooled
#   a three way split, so training can stop without touching the test set
#
# On the neighbours. The reference pooled the nearest ten across all districts
# and flagged the result as an artefact, since two non adjacent states produce
# false neighbours hundreds of kilometres apart. With only two states that
# artefact would dominate, so the nearest set is drawn within a state. A
# district's neighbours are then its actual neighbours.
#
# The split is by fiscal year, cut on t_global rather than the calendar year,
# since a fiscal year begins at biweek 7 and cutting on the year would split
# it across two.
#
#   train       April 2008 to March 2016
#   validation  April 2016 to March 2017
#   test        April 2017 to March 2020
#
# Validation exists so training can stop at the epoch that fits best without
# the test set influencing the choice. Stopping on the test loss, as the
# reference does by tuning there, means the test error is no longer held out.
#
# Direct k-step design: one matrix per horizon.
#
# Usage
#   python 09_sfnn_loader.py --k 2 --t_lag 26
# =============================================================================
import argparse
import os

import numpy as np
import pandas as pd

BASE_DEFAULT = "/.../.../.../Derived"
SAVE_DEFAULT = "/.../.../.../Derived/sfnn/prefit"

BIWEEKS_PER_YEAR = 26
BASE_YEAR = 2008          # must match the R config

# The fiscal year begins at biweek 7. t_global = (year - 2008) * 26 + biweek,
# so April 2016 is (2016 - 2008) * 26 + 7 and April 2017 one year later.
VAL_FROM  = (2016 - BASE_YEAR) * BIWEEKS_PER_YEAR + 7
TEST_FROM = (2017 - BASE_YEAR) * BIWEEKS_PER_YEAR + 7

N_NEAR_DEFAULT = 5
N_CORE_DEFAULT = 4

# Grouped by pathway rather than listed flat. Thirteen individual attributions
# would each be small and none clearly important, because the indicators are
# correlated within themes; the three vaccination measures in particular
# correlate above 0.7 with each other. Summing within a group answers whether
# the pathway matters rather than whether one indicator does.
NFHS_GROUPS = {
    "supply":    ["measles_zero", "penta_zero", "full_imm"],
    "contact":   ["hhsize", "urbanization"],
    "access":    ["women10plus", "poorest", "institutional_birth",
                  "maternal_immunisation", "sc_st"],
    "nutrition": ["stunting", "vitamin_a", "early_breastfeeding"],
}


def parse_args():
    p = argparse.ArgumentParser(description="SFNN feature loader, two states")
    p.add_argument("--base", default=BASE_DEFAULT)
    p.add_argument("--save_dir", default=SAVE_DEFAULT)
    p.add_argument("--k", type=int, default=2, help="horizon in biweeks")
    p.add_argument("--t_lag", type=int, default=26, help="how far the lags reach")
    p.add_argument("--n_near", type=int, default=N_NEAR_DEFAULT,
                   help="nearest districts, found within a state")
    p.add_argument("--n_core", type=int, default=N_CORE_DEFAULT,
                   help="largest districts by population, within a state")
    return p.parse_args()


def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def main():
    a = parse_args()
    k, t_lag = a.k, a.t_lag
    os.makedirs(a.save_dir, exist_ok=True)

    df = pd.read_csv(os.path.join(a.base, "sfnn_master.csv"))
    df = df.sort_values(["uid", "t_global"]).reset_index(drop=True)
    ids = df["uid"].drop_duplicates().tolist()
    n_dist = len(ids)
    row_of = {d: i for i, d in enumerate(ids)}
    states = sorted(df["state"].unique())

    print(f"master: {len(df)} rows, {n_dist} districts, {len(states)} states")
    print(f"k={k} t_lag={t_lag} | nearest-{a.n_near} and top-{a.n_core} within state")

    # ---- guards -------------------------------------------------------------
    # The lags shift within a district by position, so every district must share
    # the same complete grid or a shift would pull the wrong time point.
    grid = df.groupby("uid")["t_global"].agg(["min", "max", "nunique"])
    full = df["t_global"].nunique()
    assert (grid["nunique"] == full).all(), \
        "districts have unequal grids, the lag shifts would misalign"
    coord_na = df.groupby("uid")[["lat", "lon"]].first().isna().any(axis=1).sum()
    assert coord_na == 0, f"{coord_na} districts have no coordinates"
    print(f"  guards passed: a complete {full} step grid for every district")

    # ---- split --------------------------------------------------------------
    df["split"] = np.select(
        [df["t_global"] < VAL_FROM, df["t_global"] < TEST_FROM],
        ["train", "val"], default="test")
    tr_mask = df["split"].values == "train"

    print(f"  split on t_global: train < {VAL_FROM}, "
          f"validation < {TEST_FROM}, test onward")
    for nm in ("train", "val", "test"):
        n = int((df["split"] == nm).sum())
        print(f"    {nm:<6s} {n:6d} rows, {100 * n / len(df):.0f} percent")

    # ---- transform, then lag ------------------------------------------------
    # Cases and susceptibles are logged and standardised on the training rows
    # first, and the lags taken of the transformed series, so the target at t
    # and the features at t-k share a scale. Fitting the transform on training
    # alone keeps the later periods out of the earlier features.
    def make_trans(col):
        raw = np.log1p(df[col].values)
        mu_ = raw[tr_mask].mean()
        sd_ = raw[tr_mask].std() or 1.0
        return (raw - mu_) / sd_, mu_, sd_

    df["cases_trans"], cases_mu, cases_sd = make_trans("cases")
    df["susc_trans"], _, _ = make_trans("susc")

    def lag_of(col, prefix):
        tmp = pd.DataFrame({"d": df["uid"].values, "_v": df[col].values})
        g = tmp.groupby("d")["_v"]
        return pd.DataFrame(
            {f"{prefix}_lag_{i}": g.shift(i).values for i in range(k, t_lag + 1)},
            index=df.index)

    def lag_raw(col, prefix):
        tmp = pd.DataFrame({"d": df["uid"].values, "_v": df[col].values})
        g = tmp.groupby("d")["_v"]
        return pd.DataFrame(
            {f"{prefix}_lag_{i}": g.shift(i).values for i in range(k, t_lag + 1)},
            index=df.index)

    blocks = [df[["state", "district", "uid", "t_global", "time", "year",
                  "biweek", "split"]].copy()]

    blocks.append(lag_of("cases_trans", "cases"))       # own incidence
    blocks.append(lag_of("susc_trans", "susc"))         # susceptibles
    blocks.append(lag_raw("births_eff", "births"))      # two dose corrected
    blocks.append(lag_raw("population", "pop"))
    blocks.append(lag_raw("v1_lag", "v1"))              # cohort shifted already
    blocks.append(lag_raw("v2_lag", "v2"))

    # ---- geography, within state --------------------------------------------
    geo = (df.groupby("uid")
             .agg(state=("state", "first"), lat=("lat", "first"),
                  lon=("lon", "first"), pop=("population", "mean"))
             .reindex(ids).reset_index())
    lat, lon = geo["lat"].values, geo["lon"].values

    # Two distance matrices. D has an infinite diagonal so a district is never
    # its own nearest neighbour. D_self0 keeps the true zero, since a core
    # district's distance to itself is genuinely zero and an infinity there
    # would poison the column when it is standardised.
    D = np.vstack([haversine(lat[i], lon[i], lat, lon) for i in range(n_dist)])
    D_self0 = D.copy()
    np.fill_diagonal(D_self0, 0.0)
    np.fill_diagonal(D, np.inf)

    # A district in the other state is not a neighbour. Setting the cross state
    # distances to infinity confines both the nearest and the core sets within
    # a state without changing the code that finds them.
    same_state = (geo["state"].values[:, None] == geo["state"].values[None, :])
    D[~same_state] = np.inf

    inc = (df.pivot_table(index="t_global", columns="uid", values="cases_trans")
             .reindex(columns=ids))

    # core: the largest by population, within each state
    core_by_state = {}
    for s in states:
        idx = np.where(geo["state"].values == s)[0]
        order = idx[np.argsort(-geo["pop"].values[idx])][:a.n_core]
        core_by_state[s] = [ids[i] for i in order]
        print(f"  core, {s}: " + ", ".join(
            d.split("__")[1] for d in core_by_state[s]))

    # Each district looks at its own state's core, so slot ci means the ci-th
    # largest district of that state rather than a fixed district everywhere.
    state_of_row = df["state"].values
    for ci in range(a.n_core):
        vals = np.full(len(df), np.nan)
        dist = np.full(len(df), np.nan)
        for s in states:
            if ci >= len(core_by_state[s]):
                continue
            cid = core_by_state[s][ci]
            m = state_of_row == s
            vals[m] = df.loc[m, "t_global"].map(inc[cid]).values
            dist[m] = D_self0[df.loc[m, "uid"].map(row_of).values, row_of[cid]]

        tmp = pd.DataFrame({"d": df["uid"].values, "_v": vals})
        g = tmp.groupby("d")["_v"]
        cols = {f"core{ci}_lag_{i}": g.shift(i).values
                for i in range(k, t_lag + 1)}
        cols[f"core{ci}_dist"] = dist
        blocks.append(pd.DataFrame(cols, index=df.index))

    # nearest, within state by construction since cross state is infinite
    near_idx = np.argsort(D, axis=1)[:, :a.n_near]
    long_inc = inc.reset_index().melt(id_vars="t_global",
                                      var_name="neigh", value_name="ninc")

    for slot in range(a.n_near):
        neigh_for = {ids[h]: ids[near_idx[h, slot]] for h in range(n_dist)}
        sub = pd.DataFrame({"t_global": df["t_global"].values,
                            "neigh": df["uid"].map(neigh_for).values,
                            "_row": df.index.values})
        merged = (sub.merge(long_inc, on=["t_global", "neigh"], how="left")
                     .set_index("_row").reindex(df.index))
        tmp = pd.DataFrame({"d": df["uid"].values, "_v": merged["ninc"].values})
        g = tmp.groupby("d")["_v"]
        cols = {f"nbc{slot}_lag_{i}": g.shift(i).values
                for i in range(k, t_lag + 1)}
        home = D[np.arange(n_dist), near_idx[:, slot]]
        cols[f"nbc{slot}_dist"] = df["uid"].map(
            {ids[h]: home[h] for h in range(n_dist)}).values
        blocks.append(pd.DataFrame(cols, index=df.index))

    # ---- the campaign --------------------------------------------------------
    # mu_mr is the share of the pool the campaign removed in that biweek. West
    # Bengal has none, so these are zero throughout for its districts. The
    # biweeks since counter is lagged by k so it is known at prediction time,
    # and districts before their campaign take a finite sentinel rather than an
    # infinity, which would poison the standardisation.
    act = (df["mu_mr"].values > 0).astype(float)
    ta = pd.DataFrame({"d": df["uid"].values, "_a": act})
    tm = pd.DataFrame({"d": df["uid"].values, "_m": df["mu_mr"].values})
    ga, gm = ta.groupby("d")["_a"], tm.groupby("d")["_m"]

    cols = {}
    for i in range(k, t_lag + 1):
        cols[f"mr_act_lag_{i}"] = ga.shift(i).values
        cols[f"mr_mu_lag_{i}"] = gm.shift(i).values

    NEVER = float(full)
    since = np.empty(len(df))
    since[:] = np.nan
    for d, idx in df.groupby("uid").indices.items():
        seq = act[idx]
        cnt = np.full(len(idx), np.nan)
        c = NEVER
        for j in range(len(idx)):
            c = 0.0 if seq[j] > 0 else min(c + 1.0, NEVER)
            cnt[j] = c
        since[idx] = cnt
    ts = pd.DataFrame({"d": df["uid"].values, "_s": since})
    cols["mr_biweeks_since"] = ts.groupby("d")["_s"].shift(k).values
    blocks.append(pd.DataFrame(cols, index=df.index))

    # ---- static covariates, prefixed by pathway ------------------------------
    # The group is in the name, so the trainer can drop or attribute a pathway
    # without a separate lookup.
    nfhs_cols = {}
    for grp, names in NFHS_GROUPS.items():
        have = [c for c in names if c in df.columns]
        if len(have) != len(names):
            print(f"  absent from {grp}: {', '.join(set(names) - set(have))}")
        for c in have:
            nfhs_cols[f"nfhs_{grp}_{c}"] = df[c].values
    blocks.append(pd.DataFrame(nfhs_cols, index=df.index))
    print(f"  static covariates: {len(nfhs_cols)} in "
          f"{len(NFHS_GROUPS)} pathway groups")

    # ---- the state indicator -------------------------------------------------
    # One column for two states. A dummy per state would be collinear with the
    # intercept; one indicator carries the same information.
    blocks.append(pd.DataFrame(
        {f"st_{states[1].lower().replace(' ', '_')}":
            (df["state"].values == states[1]).astype(float)},
        index=df.index))

    # ---- persistence ---------------------------------------------------------
    # The last value known at prediction time. Already present as cases_lag_k,
    # so it adds nothing as a feature, but it is written out separately so the
    # model can be scored against carrying that value forward. Skill against
    # the training mean is a weak claim where most biweeks are zero.
    tmp = pd.DataFrame({"d": df["uid"].values, "_v": df["cases_trans"].values})
    blocks.append(pd.DataFrame(
        {"persistence": tmp.groupby("d")["_v"].shift(k).values}, index=df.index))

    # ---- assemble -----------------------------------------------------------
    X = pd.concat(blocks, axis=1)
    X["cases_trans"] = df["cases_trans"].values

    meta = ["state", "district", "uid", "t_global", "time", "year", "biweek",
            "split", "cases_trans", "persistence"]
    feat_cols = [c for c in X.columns if c not in meta]

    before = len(X)
    X = X.dropna(subset=feat_cols + ["cases_trans"]).reset_index(drop=True)
    print(f"\nafter the warm up trim: {len(X)} rows of {before}, "
          f"{len(feat_cols)} features")

    # ---- standardise, on the training rows only -----------------------------
    no_std = [c for c in feat_cols if c.startswith(("cases_lag_", "susc_lag_"))]
    std_cols = [c for c in feat_cols if c not in no_std]
    tr = X["split"] == "train"
    mu = X.loc[tr, std_cols].mean()
    sd = X.loc[tr, std_cols].std().replace(0, 1.0)
    X[std_cols] = (X[std_cols] - mu) / sd
    print(f"  standardised {len(std_cols)}, left as they were {len(no_std)}")

    # ---- by group ------------------------------------------------------------
    cl = list(X.columns)
    counts = {
        "own incidence": sum(c.startswith("cases_lag_") for c in cl),
        "susceptibles": sum(c.startswith("susc_lag_") for c in cl),
        "births": sum(c.startswith("births_lag_") for c in cl),
        "population": sum(c.startswith("pop_lag_") for c in cl),
        "MCV1": sum(c.startswith("v1_lag_") for c in cl),
        "MCV2": sum(c.startswith("v2_lag_") for c in cl),
        "core districts": sum(c.startswith("core") for c in cl),
        "nearest districts": sum(c.startswith("nbc") for c in cl),
        "campaign": sum(c.startswith("mr_") for c in cl),
        "static covariates": sum(c.startswith("nfhs_") for c in cl),
        "state": sum(c.startswith("st_") for c in cl),
    }
    print("\nfeatures by group:")
    for nm, n in counts.items():
        if n:
            print(f"  {nm:<20s} {n}")
    print(f"  {'total':<20s} {len(feat_cols)}")

    # ---- how far the neighbours reach ---------------------------------------
    # Within state now, so these should be genuinely local. A large distance to
    # the nth neighbour would mean a district on a state edge with few close
    # ones, which is worth knowing rather than assuming away.
    finite = np.where(np.isinf(D), np.nan, D)
    diag = pd.DataFrame({
        "uid": ids,
        "state": geo["state"].values,
        "dist_to_nearest": np.nanmin(finite, axis=1),
        f"dist_to_{a.n_near}th": np.sort(finite, axis=1)[:, a.n_near - 1]})
    diag.to_csv(os.path.join(a.save_dir, "spatial_neighbour_diagnostic.csv"),
                index=False)
    print("\nneighbour distances, km:")
    for s in states:
        sub = diag[diag["state"] == s]
        print(f"  {s}: nearest median {sub['dist_to_nearest'].median():.0f}, "
              f"{a.n_near}th median {sub[f'dist_to_{a.n_near}th'].median():.0f}, "
              f"max {sub[f'dist_to_{a.n_near}th'].max():.0f}")

    # ---- what persistence alone scores ---------------------------------------
    # If the network cannot beat this the extra machinery is not earning its
    # place, and it is better to know now than after the fitting.
    print()
    mean_ref = X.loc[X["split"] == "train", "cases_trans"].mean()
    for nm in ("train", "val", "test"):
        m = X["split"] == nm
        if not m.any():
            continue
        pmse = ((X.loc[m, "persistence"] - X.loc[m, "cases_trans"]) ** 2).mean()
        naive = ((X.loc[m, "cases_trans"] - mean_ref) ** 2).mean()
        print(f"  persistence on {nm:<5s}: mse {pmse:.4f}   "
              f"against the training mean at {naive:.4f}   "
              f"skill {1 - pmse / naive:+.3f}")

    # ---- save ---------------------------------------------------------------
    out = os.path.join(a.save_dir, f"k{k}_tlag{t_lag}.gzip")
    X.to_parquet(out, index=False, compression="gzip")

    pd.DataFrame({"cases_mean": [cases_mu], "cases_std": [cases_sd]}).to_parquet(
        os.path.join(a.save_dir, f"k{k}_tlag{t_lag}_cases_transform.gzip"),
        index=False, compression="gzip")
    pd.DataFrame({"feature": std_cols, "mean": mu.values,
                  "sd": sd.values}).to_parquet(
        os.path.join(a.save_dir, f"k{k}_tlag{t_lag}_feature_transform.gzip"),
        index=False, compression="gzip")

    print(f"\nsaved {out}")
    for nm in ("train", "val", "test"):
        print(f"  {nm:<6s} {int((X['split'] == nm).sum())}")


if __name__ == "__main__":
    main()
