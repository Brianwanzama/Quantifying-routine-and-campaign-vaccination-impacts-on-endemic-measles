#!/usr/bin/env python3
# =============================================================================
# 10_sfnn_functions.py
#
# The dataset wrapper, the network, and the train and test loops. Imported by
# the tuning, training and explanation scripts so all three see the same
# feature selection and the same architecture.
#
# The Data class and the network follow Madden et al. (2024): a linear input
# layer, then hidden layers with rectified linear activation, then a linear
# output of dimension one.
#
# The feature groups matter for two things. Ablation, where a group is dropped
# from the same matrix so two runs differ only in what was removed. And grouped
# Shapley values, where attributions are summed within a group so the answer is
# about the pathway rather than an individual lag.
#
# The static covariates carry their pathway in the column name, so
# nfhs_supply_measles_zero belongs to the supply group. Thirteen individual
# attributions would each be small and none clearly important, because the
# indicators correlate within themes.
# =============================================================================
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader


class Data(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X.values).float()
        self.y = torch.from_numpy(y).float()
        self.len = self.X.shape[0]

    def __getitem__(self, i):
        return self.X[i], self.y[i]

    def __len__(self):
        return self.len


class NeuralNetwork(nn.Module):
    """Feed forward, rectified linear, a single linear output."""

    def __init__(self, input_dim, hidden_dim, output_dim, num_hidden_layers=3):
        super(NeuralNetwork, self).__init__()
        layers = [nn.Linear(input_dim, hidden_dim), nn.ReLU()]
        for _ in range(num_hidden_layers):
            layers += [nn.Linear(hidden_dim, hidden_dim), nn.ReLU()]
        layers.append(nn.Linear(hidden_dim, output_dim))
        self.linear_relu_stack = nn.Sequential(*layers)

    def forward(self, x):
        return self.linear_relu_stack(x)


# ---------------------------------------------------------- feature selection --
#   lag_   the lagged blocks, own incidence through to the campaign
#   dist_  the distances to the core and nearest districts
#   nfhs_  the static covariates
#   ^mr_   the campaign
#   ^st_   the state indicator
#
# persistence is deliberately absent. It is a benchmark on the target scale
# rather than a feature, and the loader lists it with the identifiers.
FEATURE_REGEX = r"lag_|dist_|nfhs_|^mr_|^st_"

# the pathway groups the loader prefixed onto the covariate names
NFHS_PATHWAYS = ["supply", "contact", "access", "nutrition"]


def feature_groups(columns):
    """Groups by pathway, for ablation and for grouped Shapley values.

    A lag of a district's own incidence and a lag of its neighbour's are
    different pathways even though both are incidence, so they are separated
    rather than pooled by variable type.
    """
    g = {
        "own":        [c for c in columns if c.startswith("cases_lag_")],
        "susc":       [c for c in columns if c.startswith("susc_lag_")],
        "spatial":    [c for c in columns if c.startswith(("core", "nbc"))],
        "core":       [c for c in columns if c.startswith("core")],
        "nearest":    [c for c in columns if c.startswith("nbc")],
        "demography": [c for c in columns if c.startswith(("births_lag_", "pop_lag_"))],
        "vacc":       [c for c in columns if c.startswith(("v1_lag_", "v2_lag_"))],
        "campaign":   [c for c in columns if c.startswith("mr_")],
        "nfhs":       [c for c in columns if c.startswith("nfhs_")],
        "state":      [c for c in columns if c.startswith("st_")],
    }
    # each covariate pathway on its own, so Shapley can attribute to one
    for p in NFHS_PATHWAYS:
        g[f"nfhs_{p}"] = [c for c in columns if c.startswith(f"nfhs_{p}_")]
    return g


# The groups reported and used for Shapley. spatial is the union of core and
# nearest, and nfhs the union of the four pathways, so listing the parts as
# well would double count.
REPORT_GROUPS = ["own", "susc", "spatial", "demography", "vacc",
                 "campaign", "nfhs", "state"]

# for the Shapley chart, where the covariate pathways are worth separating
SHAP_GROUPS = ["own", "susc", "core", "nearest", "demography", "vacc",
               "campaign", "state"] + [f"nfhs_{p}" for p in NFHS_PATHWAYS]


def process_data(cases, drop_group=None, verbose=True):
    """Split three ways, select features, wrap for the dataloader.

    The loader writes the split column, so the boundaries are set there.
    """
    if "split" not in cases.columns:
        raise ValueError("no split column, rebuild the matrix with the loader")

    parts = {nm: cases[cases["split"] == nm].copy()
             for nm in ("train", "val", "test")}

    id_cols = [c for c in ["state", "district", "uid", "t_global", "time",
                           "year", "biweek"] if c in cases.columns]

    ids = {nm: d[id_cols].reset_index(drop=True) for nm, d in parts.items()}

    # the persistence benchmark rides along with the identifiers, since it is
    # on the target scale and is scored against rather than fitted on
    for nm, d in parts.items():
        if "persistence" in d.columns:
            ids[nm]["persistence"] = d["persistence"].to_numpy()

    ys = {nm: d["cases_trans"].to_numpy().reshape((-1, 1))
          for nm, d in parts.items()}
    Xs = {nm: d.filter(regex=FEATURE_REGEX) for nm, d in parts.items()}

    # Ablation drops a group from every split alike, so two runs differ only in
    # what was removed rather than in how the matrix was built.
    if drop_group is not None:
        groups = feature_groups(Xs["train"].columns)
        if drop_group not in groups:
            raise ValueError(f"unknown group '{drop_group}', "
                             f"options: {sorted(groups)}")
        cols = groups[drop_group]
        if not cols:
            raise ValueError(f"group '{drop_group}' has no columns to drop")
        Xs = {nm: X.drop(columns=cols) for nm, X in Xs.items()}
        if verbose:
            print(f"  ablation: dropped '{drop_group}', {len(cols)} columns")

    ref = list(Xs["train"].columns)
    for nm, X in Xs.items():
        if len(X) and list(X.columns) != ref:
            raise ValueError(f"the {nm} columns differ from the training ones")

    Xs = {nm: X.replace([np.inf, -np.inf], 0).fillna(0) for nm, X in Xs.items()}

    if verbose:
        g = feature_groups(Xs["train"].columns)
        print(f"  train {len(Xs['train'])} | validation {len(Xs['val'])} | "
              f"test {len(Xs['test'])} rows | {Xs['train'].shape[1]} features")
        for nm in REPORT_GROUPS:
            if g[nm]:
                print(f"    {nm:<12s} {len(g[nm])}")
        total = sum(len(g[nm]) for nm in REPORT_GROUPS)
        assert total == Xs["train"].shape[1], \
            f"the groups account for {total} of {Xs['train'].shape[1]} features"

    data = {nm: Data(Xs[nm], ys[nm]) for nm in Xs}
    return data, Xs["train"].shape[1], ids


def get_dataloaders(data, batch_size):
    """One loader per split. Only training is shuffled."""
    return {nm: DataLoader(dataset=d, batch_size=batch_size,
                           shuffle=(nm == "train"))
            for nm, d in data.items() if len(d) > 0}


def train(model, device, loader, optimizer, loss_fn=None):
    if loss_fn is None:
        loss_fn = nn.MSELoss()
    model.train()
    total, n = 0.0, 0
    for X, y in loader:
        X, y = X.to(device), y.to(device)
        optimizer.zero_grad()
        loss = loss_fn(model(X), y)
        loss.backward()
        optimizer.step()
        total += loss.item() * X.size(0)
        n += X.size(0)
    return total / max(n, 1)


def test(model, device, loader, loss_fn=None):
    if loss_fn is None:
        loss_fn = nn.MSELoss()
    model.eval()
    total, n = 0.0, 0
    with torch.no_grad():
        for X, y in loader:
            X, y = X.to(device), y.to(device)
            total += loss_fn(model(X), y).item() * X.size(0)
            n += X.size(0)
    return total / max(n, 1)


def fit_with_early_stopping(model, device, loaders, optimizer, max_epochs,
                            patience=15, loss_fn=None, verbose=False):
    """Fit until the validation loss stops improving, keeping the best state.

    Training a fixed number of epochs and keeping the last one keeps whatever
    state the model ended in, which is usually not the best it reached. The
    stopping epoch is chosen on validation so the test set plays no part in
    selecting the model.
    """
    if loss_fn is None:
        loss_fn = nn.MSELoss()

    has_val = "val" in loaders and len(loaders["val"].dataset) > 0
    if not has_val and verbose:
        print("  no validation split, running the full epoch budget")

    best_loss, best_epoch, best_state, since_best = float("inf"), 0, None, 0
    history = []

    for epoch in range(1, max_epochs + 1):
        tr = train(model, device, loaders["train"], optimizer, loss_fn)
        va = test(model, device, loaders["val"], loss_fn) if has_val else float("nan")
        te = (test(model, device, loaders["test"], loss_fn)
              if "test" in loaders else float("nan"))
        history.append({"epoch": epoch, "train_mse": tr,
                        "val_mse": va, "test_mse": te})

        if has_val:
            if va < best_loss - 1e-6:
                best_loss, best_epoch, since_best = va, epoch, 0
                best_state = {k: v.detach().clone()
                              for k, v in model.state_dict().items()}
            else:
                since_best += 1
                if since_best >= patience:
                    if verbose:
                        print(f"  stopped at {epoch}, best was {best_epoch}")
                    break

    if best_state is not None:
        model.load_state_dict(best_state)
    else:
        best_epoch = len(history)

    return model, pd.DataFrame(history), best_epoch, best_loss


def back_transform(y, cases_mean, cases_std):
    """From the modelled scale to counts, floored at zero."""
    return np.maximum(np.expm1(np.asarray(y).flatten() * cases_std + cases_mean), 0.0)
