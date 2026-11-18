# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import count_top_predicted_tfs as c


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"
FIG_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "regulator_prioritization"
RECALIBRATION_SCORES = TABLE_DIR / "thp1_k562_lfdr_recalibration_scores.csv"

METHOD_ORDER = [
    "ChronoPack_K562",
    "GENIE3",
    "SINCERITIES",
    "GRNBOOST2",
    "PPCOR",
    "Replogle2022_K562",
    "ChronoPack_THP1",
    "Yao2023_THP1",
]
METHOD_ORDER = fig5_paths.insert_prior(METHOD_ORDER)  # Fig5_GRN_v3

BUDGET_FRACTIONS = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]

METHOD_LABELS = {
    "ChronoPack_K562": "Chrono K562",
    "GENIE3": "GENIE3",
    "SINCERITIES": "SINCERITIES",
    "GRNBOOST2": "GRNBOOST2",
    "PPCOR": "PPCOR",
    "Replogle2022_K562": "Replogle",
    "ChronoPack_THP1": "Chrono THP1",
    "Yao2023_THP1": "Yao",
    **fig5_paths.PRIOR_LABELS,
}


def load_chronopack_beta_lfdr(
    name: str,
    dataset: str,
    omnipath_path: Path,
    valid_targets: set[str],
) -> pd.DataFrame:
    scores = pd.read_csv(RECALIBRATION_SCORES)
    genes = scores[scores["dataset"] == dataset][["gene", "lfdr_bum"]].copy()
    genes["target"] = genes["gene"].astype(str)
    genes["gene_score"] = 1 - pd.to_numeric(genes["lfdr_bum"], errors="coerce")
    genes = genes[genes["target"].isin(valid_targets)]
    genes = genes[genes["gene_score"] > 0].dropna(subset=["gene_score"])

    omni = pd.read_csv(omnipath_path)
    omni["TF"] = omni["TF"].astype(str).str.strip()
    omni["priority_act"] = pd.to_numeric(omni["priority_act"], errors="coerce")
    omni = omni[omni["priority_act"] > 0].dropna(subset=["priority_act"])
    omni["tf_score"] = omni["priority_act"] / omni["priority_act"].sum()

    edge_weight = np.outer(omni["tf_score"].to_numpy(float), genes["gene_score"].to_numpy(float))
    pred = (
        pd.DataFrame(edge_weight, index=omni["TF"], columns=genes["target"])
        .stack()
        .rename("EdgeWeight")
        .reset_index()
        .rename(columns={"level_0": "TF", "level_1": "target"})
    )
    return c.break_ties(c.normalize_edges(pred, name, valid_targets))


def load_recalibrated_methods(valid_targets: set[str]) -> dict[str, pd.DataFrame]:
    methods = {
        name: c.load_beeline_method(name, valid_targets)
        for name in ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]
    }
    methods["ChronoPack_K562"] = load_chronopack_beta_lfdr(
        "ChronoPack_K562", "K562", c.K562_OMNIPATH_PATH, valid_targets
    )
    methods["ChronoPack_THP1"] = load_chronopack_beta_lfdr(
        "ChronoPack_THP1", "THP1", c.THP1_OMNIPATH_PATH, valid_targets
    )
    methods["Replogle2022_K562"] = c.load_replogle(valid_targets)
    methods["Yao2023_THP1"] = c.load_yao(valid_targets)
    if fig5_paths.WITH_PRIOR:
        methods.update({name: c.load_prior_method(name, valid_targets) for name in fig5_paths.PRIOR_METHODS})
    return {method: methods[method] for method in METHOD_ORDER if method in methods}


def common_universe(methods: dict[str, pd.DataFrame], reference_tfs: set[str]) -> tuple[list[str], list[str]]:
    tf_sets = []
    target_sets = []
    for method in METHOD_ORDER:
        edges = methods[method]
        edges = edges[edges["TF"].isin(reference_tfs)]
        tf_sets.append(set(edges["TF"]))
        target_sets.append(set(edges["target"]))
    tfs = sorted(set.intersection(*tf_sets))
    targets = sorted(set.intersection(*target_sets))
    return tfs, targets


def method_universes(methods: dict[str, pd.DataFrame], reference_tfs: set[str]) -> dict[str, tuple[list[str], list[str]]]:
    universes = {}
    for method in METHOD_ORDER:
        edges = methods[method]
        edges = edges[edges["TF"].isin(reference_tfs)]
        universes[method] = (sorted(edges["TF"].unique()), sorted(edges["target"].unique()))
    return universes


def summarize_prioritization(
    methods: dict[str, pd.DataFrame],
    universes: dict[str, tuple[list[str], list[str]]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    rank_rows = []
    for method in METHOD_ORDER:
        universe_tfs, universe_targets = universes[method]
        universe_tfs_set = set(universe_tfs)
        universe_targets_set = set(universe_targets)

        pred = methods[method]
        pred = pred[pred["TF"].isin(universe_tfs_set) & pred["target"].isin(universe_targets_set)].copy()
        pred = (
            pred.sort_values("EdgeWeight", ascending=False, kind="mergesort")
            .drop_duplicates(["TF", "target"], keep="first")
            .reset_index(drop=True)
        )
        pred["rank"] = np.arange(1, len(pred) + 1)
        pred["percentile"] = 1 - (pred["rank"] - 1) / max(len(pred), 1)
        rank_rows.append(pred[["method", "TF", "target", "EdgeWeight", "rank", "percentile"]])
        ranked_edge_count = len(pred)
        possible_edge_count = len(universe_tfs) * len(universe_targets)
        budget_edges = {frac: max(1, int(round(ranked_edge_count * frac))) for frac in BUDGET_FRACTIONS}

        for frac, budget in budget_edges.items():
            top = pred.head(budget)
            counts = top.groupby("TF")["target"].nunique()
            for tf in universe_tfs:
                out_degree = int(counts.get(tf, 0))
                rows.append(
                    {
                        "method": method,
                        "method_label": METHOD_LABELS.get(method, method),
                        "TF": tf,
                        "budget_fraction": frac,
                        "budget_label": f"{100 * frac:g}%",
                        "budget_edges": budget,
                        "ranked_edges": ranked_edge_count,
                        "possible_edges": possible_edge_count,
                        "n_universe_tfs": len(universe_tfs),
                        "n_universe_targets": len(universe_targets),
                        "out_degree": out_degree,
                        "edge_share": out_degree / budget,
                    }
                )
    summary = pd.DataFrame(rows)
    ranked_edges = pd.concat(rank_rows, ignore_index=True)
    return summary, ranked_edges


def plot_heatmaps(summary: pd.DataFrame, value_col: str, stem: str, title: str, cmap: str) -> None:
    sns.set_theme(style="white", font_scale=0.86)
    budget_labels = [f"{100 * frac:g}%" for frac in BUDGET_FRACTIONS]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2), sharex=True, sharey=True)
    axes = axes.ravel()
    vmax = summary[value_col].max() if value_col == "edge_share" else None
    for ax, budget_label in zip(axes, budget_labels):
        sub = summary[summary["budget_label"] == budget_label].copy()
        top_tfs = (
            sub.groupby("TF")[value_col]
            .max()
            .sort_values(ascending=False)
            .head(24)
            .index
        )
        mat = (
            sub[sub["TF"].isin(top_tfs)]
            .pivot(index="method_label", columns="TF", values=value_col)
            .reindex(index=[METHOD_LABELS[m] for m in METHOD_ORDER], columns=list(top_tfs))
            .fillna(0)
        )
        sns.heatmap(
            mat,
            ax=ax,
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            cbar=budget_label == budget_labels[-1],
            linewidths=0.5,
            linecolor="white",
        )
        budget_edges = int(sub["budget_edges"].iloc[0])
        ax.set_title(f"Top {budget_label} ({budget_edges} edges)")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)
    fig.suptitle(title, y=1.02, fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_tf_share_lines(summary: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.95)
    top_tfs = (
        summary.groupby("TF")["edge_share"]
        .max()
        .sort_values(ascending=False)
        .head(12)
        .index
        .tolist()
    )
    fig, axes = plt.subplots(3, 4, figsize=(15.0, 9.5), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, tf in zip(axes, top_tfs):
        sub = summary[summary["TF"] == tf].copy()
        sub["budget_percent"] = 100 * sub["budget_fraction"]
        sns.lineplot(
            data=sub,
            x="budget_percent",
            y="edge_share",
            hue="method_label",
            hue_order=[METHOD_LABELS[m] for m in METHOD_ORDER],
            marker="o",
            linewidth=1.6,
            ax=ax,
            legend=False,
        )
        ax.set_xscale("log")
        ax.set_title(tf)
        ax.set_xlabel("Top edge budget (%)")
        ax.set_ylabel("Edge share")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for ax in axes[len(top_tfs):]:
        ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    if not handles:
        handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.suptitle("TF Edge Share Across Method-Specific Edge Budgets", y=1.02, fontsize=14)
    fig.tight_layout(rect=(0, 0, 0.86, 1))
    fig.savefig(FIG_DIR / "tf_edge_share_budget_lines.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "tf_edge_share_budget_lines.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    valid_targets = c.load_valid_targets()
    reference_tfs = c.load_reference_tfs(valid_targets)
    methods = load_recalibrated_methods(valid_targets)
    common_tfs, common_targets = common_universe(methods, reference_tfs)
    universes = method_universes(methods, reference_tfs)
    summary, ranked_edges = summarize_prioritization(methods, universes)

    universe_rows = []
    for method in METHOD_ORDER:
        tfs, targets = universes[method]
        universe_rows.append(
            {
                "method": method,
                "method_label": METHOD_LABELS.get(method, method),
                "n_tfs": len(tfs),
                "n_targets": len(targets),
                "n_edges": len(tfs) * len(targets),
                "tfs": ",".join(tfs),
                "budget_fractions": ",".join(str(x) for x in BUDGET_FRACTIONS),
            }
        )
    universe = pd.DataFrame(universe_rows)
    common_universe_df = pd.DataFrame(
        [
            {
                "n_common_tfs": len(common_tfs),
                "common_tfs": ",".join(common_tfs),
                "n_common_targets": len(common_targets),
                "n_common_edges": len(common_tfs) * len(common_targets),
            }
        ]
    )
    universe_out = TABLE_DIR / "regulator_prioritization_method_specific_universes.csv"
    common_universe_out = TABLE_DIR / "regulator_prioritization_common_universe_reference.csv"
    summary_out = TABLE_DIR / "regulator_prioritization_method_specific_tf_budget_summary.csv"
    ranked_out = TABLE_DIR / "regulator_prioritization_method_specific_ranked_edges.csv"
    universe.to_csv(universe_out, index=False)
    common_universe_df.to_csv(common_universe_out, index=False)
    summary.to_csv(summary_out, index=False)
    ranked_edges.to_csv(ranked_out, index=False)

    plot_heatmaps(
        summary,
        "edge_share",
        "method_specific_tf_edge_share_budget_heatmaps",
        "TF Edge Share in Top-Ranked Method-Specific Edges",
        "viridis",
    )
    plot_heatmaps(
        summary,
        "out_degree",
        "method_specific_tf_out_degree_budget_heatmaps",
        "TF Out-Degree in Top-Ranked Method-Specific Edges",
        "mako",
    )
    plot_tf_share_lines(summary)

    print(f"Wrote {universe_out}")
    print(f"Wrote {common_universe_out}")
    print(f"Wrote {summary_out}")
    print(f"Wrote {ranked_out}")
    print(f"Wrote {FIG_DIR / 'method_specific_tf_edge_share_budget_heatmaps.png'}")
    print(f"Wrote {FIG_DIR / 'method_specific_tf_out_degree_budget_heatmaps.png'}")
    print(f"Wrote {FIG_DIR / 'tf_edge_share_budget_lines.png'}")


if __name__ == "__main__":
    main()
