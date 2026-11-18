#!/usr/bin/env python
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORKSPACE_DIR = Path(__file__).resolve().parents[1]
FIGURE_DIR = WORKSPACE_DIR / "benchmark_outputs" / "figures"
RESULTS_DIR = WORKSPACE_DIR / "benchmark_outputs" / "results"

CHRONOPACK_RANKS = FIGURE_DIR / "chronopack_zscore_composite_ligand_ranks.csv"
HARMONIZED_ALL = RESULTS_DIR / "harmonized_all_methods.csv"
LIANA_INDIVIDUAL = RESULTS_DIR / "harmonized_liana_individual_methods.csv"


def load_other_method_ligand_ranks() -> pd.DataFrame:
    """Return best ligand rank percentile per non-ChronoPack method/run."""
    all_methods = pd.read_csv(HARMONIZED_ALL)
    all_methods = all_methods.loc[~all_methods["method"].astype(str).str.startswith("ChronoPack")].copy()

    individual = pd.read_csv(LIANA_INDIVIDUAL)
    shared_cols = [col for col in all_methods.columns if col in individual.columns]
    other = pd.concat([all_methods, individual[shared_cols]], ignore_index=True, sort=False)
    other = other.dropna(subset=["method", "run_id", "ligand", "rank"]).copy()

    # Collapse receptor-resolved methods to ligand-level ranks within each method/run.
    ligand_rank = (
        other.groupby(["method", "run_id", "ligand"], as_index=False)
        .agg(method_ligand_rank=("rank", "min"))
    )
    max_rank = (
        ligand_rank.groupby(["method", "run_id"], as_index=False)
        .agg(method_max_ligand_rank=("method_ligand_rank", "max"))
    )
    ligand_rank = ligand_rank.merge(max_rank, on=["method", "run_id"], how="left")
    ligand_rank["method_rank_pct"] = (
        ligand_rank["method_ligand_rank"] / ligand_rank["method_max_ligand_rank"]
    )
    return ligand_rank


def build_specificity_table() -> pd.DataFrame:
    chrono = pd.read_csv(CHRONOPACK_RANKS).copy()
    chrono = chrono.rename(columns={"zscore_composite_rank": "chronopack_ligand_rank"})
    chrono["chronopack_rank_pct"] = chrono["chronopack_ligand_rank"] / chrono[
        "chronopack_ligand_rank"
    ].max()

    other = load_other_method_ligand_ranks()
    best_other = (
        other.sort_values(["ligand", "method_rank_pct", "method_ligand_rank"])
        .groupby("ligand", as_index=False)
        .first()
        .rename(
            columns={
                "method": "best_other_method",
                "run_id": "best_other_run_id",
                "method_ligand_rank": "best_other_ligand_rank",
                "method_max_ligand_rank": "best_other_max_ligand_rank",
                "method_rank_pct": "best_other_rank_pct",
            }
        )
    )

    table = chrono.merge(best_other, on="ligand", how="left")
    missing_other_pct = 1.05
    table["best_other_rank_pct_filled"] = table["best_other_rank_pct"].fillna(missing_other_pct)
    table["best_other_label"] = np.where(
        table["best_other_method"].notna(),
        table["best_other_method"].astype(str) + " / " + table["best_other_run_id"].astype(str),
        "not recovered by comparator methods",
    )
    table["chronopack_specificity"] = (
        table["best_other_rank_pct_filled"] - table["chronopack_rank_pct"]
    )
    table["log2_rank_advantage"] = np.log2(
        (table["best_other_rank_pct_filled"] + 1e-3) / (table["chronopack_rank_pct"] + 1e-3)
    )
    return table.sort_values("chronopack_specificity", ascending=False)


def plot_specificity(table: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    plot_df = table.sort_values("chronopack_ligand_rank").copy()
    top_specific = table.head(14).copy()

    fig, (ax_scatter, ax_bar) = plt.subplots(
        1,
        2,
        figsize=(11.2, 5.8),
        dpi=220,
        gridspec_kw={"width_ratios": [1.2, 0.95], "wspace": 0.35},
    )

    colors = np.where(plot_df["chronopack_specificity"] > 0, "#218380", "#9aa0a6")
    ax_scatter.scatter(
        plot_df["chronopack_rank_pct"],
        plot_df["best_other_rank_pct_filled"],
        c=colors,
        s=np.clip(78 * (plot_df["chronopack_specificity"].clip(lower=0) + 0.12), 18, 82),
        alpha=0.88,
        edgecolor="white",
        linewidth=0.45,
    )
    ax_scatter.plot([0, 1.05], [0, 1.05], color="0.58", lw=1.0, ls="--")
    ax_scatter.set_xlim(0, 1.07)
    ax_scatter.set_ylim(0, 1.07)
    ax_scatter.invert_xaxis()
    ax_scatter.invert_yaxis()
    ax_scatter.set_xlabel("ChronoPack ligand rank percentile (better left)", fontsize=9)
    ax_scatter.set_ylabel("Best comparator ligand rank percentile (worse downward)", fontsize=9)
    ax_scatter.set_title("ChronoPack-specific ligand recovery", loc="left", fontsize=11)
    ax_scatter.grid(True, color="#e9ecef", lw=0.75)
    ax_scatter.spines["top"].set_visible(False)
    ax_scatter.spines["right"].set_visible(False)

    for row in top_specific.itertuples(index=False):
        ax_scatter.annotate(
            row.ligand,
            (row.chronopack_rank_pct, row.best_other_rank_pct_filled),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=7.2,
            color="#1f2933",
        )

    bar_df = top_specific.sort_values("chronopack_specificity", ascending=True)
    y = np.arange(len(bar_df))
    ax_bar.barh(y, bar_df["chronopack_specificity"], color="#218380", alpha=0.9)
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(bar_df["ligand"], fontsize=8)
    ax_bar.set_xlabel("Best comparator percentile - ChronoPack percentile", fontsize=9)
    ax_bar.set_title("Largest ChronoPack advantage", loc="left", fontsize=11)
    ax_bar.grid(True, axis="x", color="#e9ecef", lw=0.75)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.spines["left"].set_visible(False)
    ax_bar.tick_params(axis="y", length=0)

    fig.suptitle(
        "Ligands ranked much earlier by ChronoPack than by LIANA, CellChat-style methods, NicheNet, or TraSig",
        x=0.02,
        ha="left",
        fontsize=12,
    )
    fig.text(
        0.02,
        0.02,
        "Comparator rank is the best ligand-level percentile across all non-ChronoPack method/run combinations; missing ligands are assigned 105th percentile.",
        fontsize=7.5,
        color="#4b5563",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.95))

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(FIGURE_DIR / f"chronopack_specific_ligand_recovery.{ext}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    table = build_specificity_table()
    out_csv = FIGURE_DIR / "chronopack_specific_ligand_recovery.csv"
    table.to_csv(out_csv, index=False)
    plot_specificity(table)
    print("Top ChronoPack-specific ligands:")
    cols = [
        "ligand",
        "chronopack_ligand_rank",
        "best_pair",
        "best_other_method",
        "best_other_run_id",
        "best_other_ligand_rank",
        "chronopack_specificity",
    ]
    print(table[cols].head(16).to_string(index=False))
    print("Wrote:")
    print(out_csv)
    for ext in ["png", "pdf", "svg"]:
        print(FIGURE_DIR / f"chronopack_specific_ligand_recovery.{ext}")


if __name__ == "__main__":
    main()
