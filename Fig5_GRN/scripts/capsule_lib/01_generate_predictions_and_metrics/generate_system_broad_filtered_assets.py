# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import seaborn as sns

import benchmark_grn_systems as b


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "system_broad_filtered"
SCENARIO = "system_broad_all_positive_tfs_five_methods"
OUT_SCENARIO = "system_broad_filtered"

METHOD_ORDER = [
    "ChronoPack_K562",
    "GENIE3",
    "Replogle2022_K562",
    "ChronoPack_THP1",
    "Yao2023_THP1",
]
METHOD_SHORT = {
    "ChronoPack_K562": "CP K562",
    "GENIE3": "GENIE3",
    "Replogle2022_K562": "Replogle",
    "ChronoPack_THP1": "CP THP1",
    "Yao2023_THP1": "Yao THP1",
}
BENCHMARK_LABELS = {"TRRUST": "TRRUST", "K562_ChIP": "K562 ChIP"}
METRIC_LABELS = {
    "EPR": "EPR",
    "AUPRC_ratio": "AUPRCr",
    "positive_coverage": "positive coverage",
    "median_EPR": "median EPR",
    "median_AUPRC_ratio": "median AUPRCr",
}
SYSTEM_BLUE = "#2d4056"
CMAP = LinearSegmentedColormap.from_list(
    "reference_purple_orange",
    [
        "#4b116f",
        "#7a1f82",
        "#a12c89",
        "#c53d78",
        "#e25163",
        "#f36e54",
        "#fb8f61",
        "#f9b873",
        "#fee2a1",
    ],
)


def save(fig: plt.Figure, stem: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def tf_order_from_sets(tf_sets: pd.DataFrame, benchmark: str) -> list[str]:
    row = tf_sets[tf_sets["benchmark"] == benchmark].iloc[0]
    return [tf for tf in str(row["TFs"]).split(",") if tf]


def normalize_over_method_rows(raw: pd.DataFrame) -> pd.DataFrame:
    arr = raw.to_numpy(dtype=float)
    out = np.full_like(arr, np.nan, dtype=float)
    for j in range(arr.shape[1]):
        col = arr[:, j]
        mask = np.isfinite(col)
        if not mask.any():
            continue
        lo = np.nanmin(col[mask])
        hi = np.nanmax(col[mask])
        out[mask, j] = 0.5 if np.isclose(lo, hi) else (col[mask] - lo) / (hi - lo)
    return pd.DataFrame(out, index=raw.index, columns=raw.columns)


def format_annotations(raw: pd.DataFrame) -> pd.DataFrame:
    return raw.map(lambda x: "" if pd.isna(x) else f"{x:.2g}")


def add_row_label_bands(ax, labels: list[str], width: float = 0.62, fontsize: float = 12) -> None:
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    for i, label in enumerate(labels):
        if i % 2 == 0:
            ax.add_patch(
                patches.Rectangle(
                    (-width, i),
                    width,
                    1,
                    transform=trans,
                    facecolor="#e6e6e6",
                    edgecolor="none",
                    clip_on=False,
                    zorder=-1,
                )
            )
        ax.text(
            -0.06,
            i + 0.5,
            label,
            transform=trans,
            ha="right",
            va="center",
            fontsize=fontsize,
            color=SYSTEM_BLUE,
            clip_on=False,
        )


def style_heatmap_axis(ax, title: str) -> None:
    ax.set_title(title, fontsize=14, pad=18)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", labelrotation=0, labelsize=12, top=False, bottom=False, labeltop=True, labelbottom=False, pad=8)
    ax.tick_params(axis="y", left=False, right=False, labelleft=False)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.plot([0, 1], [1.015, 1.015], transform=ax.transAxes, color="black", lw=2.0, clip_on=False)


def apply_annotation_colors(ax, color: pd.DataFrame) -> None:
    vals = color.to_numpy(dtype=float).ravel()
    for text, val in zip(ax.texts, vals):
        text.set_color("white" if np.isfinite(val) and val < 0.48 else "#111111")


def add_low_high_legend(fig, x: float = 0.34, y: float = 0.03, width: float = 0.32, height: float = 0.02) -> None:
    ax = fig.add_axes([x, y, width, height])
    ax.imshow(np.linspace(0, 1, 12).reshape(1, -1), aspect="auto", cmap=CMAP, interpolation="nearest")
    ax.set_axis_off()
    fig.patches.append(
        patches.Rectangle(
            (x - 0.03, y - 0.028),
            width + 0.06,
            height + 0.065,
            transform=fig.transFigure,
            facecolor="#e6e6e6",
            edgecolor="none",
            zorder=-1,
        )
    )
    fig.text(x, y - 0.012, "Low/Poor", ha="left", va="top", fontsize=14)
    fig.text(x + width, y - 0.012, "High/Good", ha="right", va="top", fontsize=14)


def plot_panel(ax, raw: pd.DataFrame, title: str, show_y_label: bool = True) -> None:
    color = normalize_over_method_rows(raw)
    sns.heatmap(
        color,
        annot=format_annotations(raw),
        fmt="",
        cmap=CMAP,
        vmin=0,
        vmax=1,
        cbar=False,
        linewidths=1.25,
        linecolor="white",
        ax=ax,
        annot_kws={"fontsize": 12},
    )
    apply_annotation_colors(ax, color)
    if show_y_label:
        add_row_label_bands(ax, [str(x) for x in raw.index])
    style_heatmap_axis(ax, title)


def per_tf_matrix(metrics: pd.DataFrame, benchmark: str, metric: str, method_order: list[str], tf_order: list[str]) -> pd.DataFrame:
    sub = metrics[
        (metrics["scenario"] == SCENARIO)
        & (metrics["benchmark"] == benchmark)
        & (metrics["scope"] == "per_tf")
        & (metrics["method"].isin(method_order))
        & (metrics["TF"].isin(tf_order))
        & (metrics["k_pos"].gt(0))
    ]
    mat = sub.pivot_table(index="method", columns="TF", values=metric, aggfunc="first")
    mat = mat.reindex(index=method_order, columns=tf_order)
    return mat.rename(index=METHOD_SHORT)


def median_matrix(summary: pd.DataFrame, metric: str, method_order: list[str]) -> pd.DataFrame:
    sub = summary[
        (summary["scenario"] == SCENARIO)
        & (summary["benchmark"].isin(["TRRUST", "K562_ChIP"]))
        & (summary["method"].isin(method_order))
    ]
    mat = sub.pivot_table(index="method", columns="benchmark", values=metric, aggfunc="first")
    mat = mat.reindex(index=method_order, columns=["TRRUST", "K562_ChIP"])
    return mat.rename(index=METHOD_SHORT, columns=BENCHMARK_LABELS)


def grouped_median_matrix(summary: pd.DataFrame, method_order: list[str]) -> pd.DataFrame:
    mats = []
    for benchmark in ["TRRUST", "K562_ChIP"]:
        cols = []
        sub = summary[
            (summary["scenario"] == SCENARIO)
            & (summary["benchmark"] == benchmark)
            & (summary["method"].isin(method_order))
        ].set_index("method")
        for metric in ["median_EPR", "median_AUPRC_ratio"]:
            cols.append(sub[metric].reindex(method_order).rename((BENCHMARK_LABELS[benchmark], METRIC_LABELS[metric].replace("median ", ""))))
        mats.append(pd.concat(cols, axis=1))
    return pd.concat(mats, axis=1).rename(index=METHOD_SHORT)


def plot_grouped_panel(ax, raw: pd.DataFrame, group_label_y: float = 1.22) -> None:
    color = normalize_over_method_rows(raw)
    sns.heatmap(
        color,
        annot=format_annotations(raw),
        fmt="",
        cmap=CMAP,
        vmin=0,
        vmax=1,
        cbar=False,
        linewidths=1.25,
        linecolor="white",
        ax=ax,
        annot_kws={"fontsize": 12},
    )
    apply_annotation_colors(ax, color)
    add_row_label_bands(ax, [str(x) for x in raw.index])
    style_heatmap_axis(ax, "")
    ax.set_xticklabels(raw.columns.get_level_values(1), rotation=0, fontsize=12)
    for group in list(dict.fromkeys(raw.columns.get_level_values(0))):
        idxs = [i for i, val in enumerate(raw.columns.get_level_values(0)) if val == group]
        center = (min(idxs) + max(idxs) + 1) / 2
        ax.text(
            center / raw.shape[1],
            group_label_y,
            group,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=14,
            fontweight="bold",
            color=SYSTEM_BLUE,
        )
    for x in range(2, raw.shape[1], 2):
        ax.axvline(x, color="white", lw=5.0)


def plot_summary(summary: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    sub = summary[summary["scenario"] == SCENARIO].copy()
    sub["method"] = pd.Categorical(sub["method"], METHOD_ORDER, ordered=True)
    sub = sub.sort_values("method")
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, metric, title in zip(axes, ["median_EPR", "median_AUPRC_ratio"], ["Median EPR", "Median AUPRC ratio"]):
        sns.barplot(data=sub, y="method_label", x=metric, hue="benchmark", ax=ax, palette={"TRRUST": "#777777", "K562_ChIP": "#4C78A8"})
        ax.axvline(1, color="black", ls="--", lw=1, alpha=0.6)
        ax.set_title(title)
        ax.set_xlabel(metric)
        ax.set_ylabel("")
        if metric == "median_AUPRC_ratio" and ax.legend_:
            ax.legend_.remove()
    axes[0].legend(frameon=False, loc="lower right")
    fig.suptitle("system_broad_filtered: median per-TF performance")
    fig.tight_layout()
    save(fig, f"summary_{OUT_SCENARIO}")


def plot_simple_metric_heatmaps(metrics: pd.DataFrame, tf_sets: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.85)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        tf_order = tf_order_from_sets(tf_sets, benchmark)
        sub = metrics[
            (metrics["scenario"] == SCENARIO)
            & (metrics["benchmark"] == benchmark)
            & (metrics["scope"] == "per_tf")
            & (metrics["method"].isin(METHOD_ORDER))
            & (metrics["TF"].isin(tf_order))
        ].copy()
        sub["method"] = pd.Categorical(sub["method"], METHOD_ORDER, ordered=True)
        for metric in ["EPR", "AUPRC_ratio", "positive_coverage"]:
            mat = sub.pivot_table(index="method", columns="TF", values=metric, aggfunc="first")
            mat = mat.reindex(index=METHOD_ORDER, columns=tf_order).rename(index=METHOD_SHORT)
            fig, ax = plt.subplots(figsize=(8, 3.3))
            sns.heatmap(mat, annot=True, fmt=".2g", cmap="viridis", ax=ax, cbar_kws={"label": metric})
            ax.set_title(f"system_broad_filtered {BENCHMARK_LABELS[benchmark]}: per-TF {metric}")
            ax.set_xlabel("")
            ax.set_ylabel("")
            fig.tight_layout()
            save(fig, f"heatmap_{OUT_SCENARIO}_{benchmark}_{metric}")


def plot_normalized_per_tf_heatmaps(metrics: pd.DataFrame, tf_sets: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.9)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        tf_order = tf_order_from_sets(tf_sets, benchmark)
        fig_width = 7.6 if len(tf_order) <= 4 else 10.0
        fig, axes = plt.subplots(1, 2, figsize=(fig_width, 3.2), sharey=True)
        for i, (ax, metric) in enumerate(zip(axes, ["EPR", "AUPRC_ratio"])):
            mat = per_tf_matrix(metrics, benchmark, metric, METHOD_ORDER, tf_order)
            plot_panel(ax, mat, METRIC_LABELS[metric], show_y_label=(i == 0))
        fig.suptitle(
            f"system_broad_filtered {BENCHMARK_LABELS[benchmark]} per-TF performance\n"
            "Filtered TFs match five_method_positive_tfs macro PR curves; colors normalized within each TF column",
            fontsize=12,
            y=1.05,
        )
        add_low_high_legend(fig)
        fig.tight_layout(rect=[0, 0.1, 1, 0.98])
        save(fig, f"heatmap_{OUT_SCENARIO}_{benchmark}_per_tf_EPR_AUPRCr_method_normalized")


def plot_median_heatmaps(summary: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.9)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for i, (ax, metric) in enumerate(zip(axes, ["median_EPR", "median_AUPRC_ratio"])):
        mat = median_matrix(summary, metric, METHOD_ORDER)
        plot_panel(ax, mat, METRIC_LABELS[metric], show_y_label=(i == 0))
    fig.suptitle(
        "system_broad_filtered median per-TF performance\n"
        "Filtered TFs match five_method_positive_tfs macro PR curves",
        fontsize=12,
        y=1.09,
    )
    add_low_high_legend(fig, x=0.31, y=0.03, width=0.38, height=0.025)
    fig.tight_layout(rect=[0, 0.12, 1, 0.98])
    save(fig, f"heatmap_{OUT_SCENARIO}_median_EPR_AUPRCr_method_normalized")

    fig, ax = plt.subplots(1, 1, figsize=(7.6, 3.9))
    plot_grouped_panel(ax, grouped_median_matrix(summary, METHOD_ORDER))
    fig.suptitle("system_broad_filtered median metrics grouped by benchmark", fontsize=12, color=SYSTEM_BLUE, y=0.98)
    add_low_high_legend(fig, x=0.31, y=0.035, width=0.38, height=0.025)
    fig.subplots_adjust(left=0.25, right=0.98, top=0.58, bottom=0.23)
    save(fig, f"heatmap_{OUT_SCENARIO}_median_grouped_by_benchmark")


def plot_network_pr_curves(curves: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.9)
    sub = curves[(curves["scenario"] == SCENARIO) & (curves["scope"] == "network")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=False)
    for ax, benchmark in zip(axes, ["TRRUST", "K562_ChIP"]):
        panel = sub[sub["benchmark"] == benchmark]
        for method in METHOD_ORDER:
            d = panel[panel["method"] == method].sort_values("k")
            if d.empty:
                continue
            obs = d[~d["is_expectation"]]
            tail = d[d["is_expectation"]]
            color = b.METHOD_COLORS[method]
            if len(obs):
                ax.plot(obs["recall"], obs["precision"], drawstyle="steps-post", lw=1.8, color=color, label=b.METHOD_LABELS[method])
                ax.scatter(obs.iloc[-1]["recall"], obs.iloc[-1]["precision"], s=18, color=color, zorder=4)
            if len(tail):
                start = pd.concat([obs.tail(1), tail.head(1)]) if len(obs) else tail.head(1)
                tail_plot = pd.concat([start, tail]).drop_duplicates(["k", "recall", "precision"])
                ax.plot(tail_plot["recall"], tail_plot["precision"], lw=1.0, ls="--", color=color, alpha=0.55)
        ax.set_title(f"{BENCHMARK_LABELS[benchmark]} - filtered network PR")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_xlim(0, 1)
        ymax = panel["precision"].replace([np.inf, -np.inf], np.nan).dropna().quantile(0.995) if len(panel) else 1
        ax.set_ylim(0, min(1, max(0.05, ymax * 1.1)))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.suptitle("Precision-recall curves: system_broad_filtered (network)")
    fig.tight_layout(rect=[0, 0, 0.78, 0.94])
    save(fig, f"pr_curves_{OUT_SCENARIO}_network")


def plot_macro_curves(macro: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.9)
    for agg in ["mean_equal_TF", "median_equal_TF"]:
        sub = macro[(macro["scenario"] == SCENARIO) & (macro["aggregation"] == agg)].copy()
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=False)
        for ax, benchmark in zip(axes, ["TRRUST", "K562_ChIP"]):
            panel = sub[sub["benchmark"] == benchmark]
            for method in METHOD_ORDER:
                d = panel[panel["method"] == method].sort_values("recall")
                if d.empty:
                    continue
                ax.plot(d["recall"], d["precision"], lw=1.8, color=b.METHOD_COLORS[method], label=b.METHOD_LABELS[method])
            tf_label = panel["TFs"].iloc[0] if len(panel) else ""
            ax.set_title(f"{BENCHMARK_LABELS[benchmark]}: {agg.replace('_', ' ')} macro PR\nTFs: {tf_label}")
            ax.set_xlabel("Recall")
            ax.set_ylabel("Precision")
            ax.set_xlim(0, 1)
            ymax = panel["precision"].replace([np.inf, -np.inf], np.nan).dropna().quantile(0.995) if len(panel) else 1
            ax.set_ylim(0, min(1, max(0.05, ymax * 1.1)))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
        fig.tight_layout(rect=[0, 0, 0.78, 1])
        save(fig, f"{OUT_SCENARIO}_{agg}_macro_pr")


def compute_per_tf_curves(tf_sets: pd.DataFrame) -> pd.DataFrame:
    valid_targets = b.load_valid_targets()
    gts = b.load_ground_truths(valid_targets)
    methods = b.load_methods(valid_targets)
    rows = []
    for benchmark in ["TRRUST", "K562_ChIP"]:
        tf_order = tf_order_from_sets(tf_sets, benchmark)
        for tf in tf_order:
            for method in METHOD_ORDER:
                _, curve = b.evaluate_one(methods[method].edges, gts[benchmark], {tf}, valid_targets, curve_points=500)
                c = curve.copy()
                c["scenario"] = SCENARIO
                c["benchmark"] = benchmark
                c["scope"] = "per_tf"
                c["TF"] = tf
                c["method"] = method
                rows.append(c)
    return pd.concat(rows, ignore_index=True)


def plot_per_tf_pr_curves(per_tf_curves: pd.DataFrame, tf_sets: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.85)
    max_tfs = max(len(tf_order_from_sets(tf_sets, bench)) for bench in ["TRRUST", "K562_ChIP"])
    fig, axes = plt.subplots(max_tfs, 2, figsize=(9, 2.8 * max_tfs), squeeze=False)
    for j, benchmark in enumerate(["TRRUST", "K562_ChIP"]):
        tf_order = tf_order_from_sets(tf_sets, benchmark)
        for i in range(max_tfs):
            ax = axes[i, j]
            if i >= len(tf_order):
                ax.set_axis_off()
                continue
            tf = tf_order[i]
            panel = per_tf_curves[(per_tf_curves["benchmark"] == benchmark) & (per_tf_curves["TF"] == tf) & (~per_tf_curves["is_expectation"])]
            for method in METHOD_ORDER:
                d = panel[panel["method"] == method].sort_values("k")
                if d.empty:
                    continue
                ax.plot(d["recall"], d["precision"], drawstyle="steps-post", lw=1.2, color=b.METHOD_COLORS[method], label=b.METHOD_LABELS[method])
            ax.set_title(f"{BENCHMARK_LABELS[benchmark]} - {tf}")
            ax.set_xlabel("Recall")
            ax.set_ylabel("Precision")
            ax.set_xlim(0, 1)
            ymax = panel["precision"].replace([np.inf, -np.inf], np.nan).dropna().quantile(0.995) if len(panel) else 1
            ax.set_ylim(0, min(1, max(0.05, ymax * 1.1)))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.suptitle("Precision-recall curves: system_broad_filtered (per TF)")
    fig.tight_layout(rect=[0, 0, 0.84, 0.97])
    save(fig, f"pr_curves_{OUT_SCENARIO}_per_tf")


def main() -> None:
    metrics = pd.read_csv(TABLE_DIR / "metrics_five_method_positive_tfs.csv")
    summary = pd.read_csv(TABLE_DIR / "summary_five_method_positive_tfs.csv")
    macro = pd.read_csv(TABLE_DIR / "macro_pr_curves_five_method_positive_tfs.csv")
    network_curves = pd.read_csv(TABLE_DIR / "micro_pr_curves_five_method_positive_tfs.csv")
    tf_sets = pd.read_csv(TABLE_DIR / "five_method_positive_tf_sets.csv")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_summary(summary)
    plot_simple_metric_heatmaps(metrics, tf_sets)
    plot_normalized_per_tf_heatmaps(metrics, tf_sets)
    plot_median_heatmaps(summary)
    plot_network_pr_curves(network_curves)
    plot_macro_curves(macro)
    per_tf_curves = compute_per_tf_curves(tf_sets)
    per_tf_curves.to_csv(TABLE_DIR / "system_broad_filtered_per_tf_pr_curves.csv", index=False)
    plot_per_tf_pr_curves(per_tf_curves, tf_sets)

    print("Wrote system_broad_filtered assets to", OUT_DIR)


if __name__ == "__main__":
    main()
