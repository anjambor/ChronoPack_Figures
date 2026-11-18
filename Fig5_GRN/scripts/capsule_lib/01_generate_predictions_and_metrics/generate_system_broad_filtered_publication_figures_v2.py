# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from pathlib import Path

import fig5_paths
import pickle

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns

import benchmark_grn_systems as b
from generate_system_broad_filtered_assets import (
    BENCHMARK_LABELS,
    CMAP,
    SYSTEM_BLUE,
    add_low_high_legend,
    apply_annotation_colors,
    format_annotations,
    normalize_over_method_rows,
    style_heatmap_axis,
)
from generate_system_broad_filtered_publication_figures import (
    PR_METHOD_ORDER,
    VIOLIN_METHOD_LABELS,
    add_cell_type_row_groups,
    add_compact_row_label_bands,
    add_trrust_tf_cell_type_row_groups,
)


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"
OUT_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "legacy_publication"
V2_TABLE_DIR = TABLE_DIR / "publication_figures"
SCENARIO = "system_broad_filtered_publication_v2"
ALL_METHOD_ORDER = PR_METHOD_ORDER

METHOD_CONTEXT = {
    "ChronoPack_K562": "K562",
    "GENIE3": "K562",
    "Replogle2022_K562": "K562",
    "ChronoPack_THP1": "THP1",
    "Yao2023_THP1": "THP1",
    **{m: "K562" for m in fig5_paths.PRIOR_METHODS},
}

BASE_TFS_BY_BENCHMARK = {
    "TRRUST": ["RELA", "NFKB1", "JUN", "ELK1", "CREB1", "ATF1", "FOS", "ATF2"],
    "K562_ChIP": ["RELA", "JUN", "ELK1", "CREB1", "ATF1", "FOS", "ATF2"],
}
THP1_EXTRA_TFS = ["REL", "IRF3"]
EPR_CAP = 25.0
RANK_COVERAGE_CUTOFFS = ["Benchmark top-K", "Top 5%", "Top 1000"]

MILD_VIRIDIS = LinearSegmentedColormap.from_list(
    "mild_viridis_v2",
    plt.get_cmap("viridis")(np.linspace(0.0, 0.90, 256)),
)


def save(fig: plt.Figure, stem: str, bbox_inches: str | None = "tight") -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches=bbox_inches)
    fig.savefig(OUT_DIR / f"{stem}.pdf", dpi=300, bbox_inches=bbox_inches)
    plt.close(fig)


def save_existing_figure(fig: plt.Figure, stems: list[str], bbox_inches: str | None = "tight") -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for stem in stems:
        fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches=bbox_inches)
        fig.savefig(OUT_DIR / f"{stem}.pdf", dpi=300, bbox_inches=bbox_inches)
    plt.close(fig)


def method_requested_tfs(benchmark: str, method: str) -> list[str]:
    tfs = list(BASE_TFS_BY_BENCHMARK[benchmark])
    if METHOD_CONTEXT[method] == "THP1":
        tfs.extend(THP1_EXTRA_TFS)
    return list(dict.fromkeys(tfs))


def ordered_tfs_for_benchmark(benchmark: str) -> list[str]:
    out = []
    for method in ALL_METHOD_ORDER:
        out.extend(method_requested_tfs(benchmark, method))
    return list(dict.fromkeys(out))


def load_reference_ground_truths(valid_targets: set[str]) -> dict[str, pd.DataFrame]:
    trrust = pd.read_csv(b.TRRUST_PATH, sep="\t", header=None)
    trrust.columns = ["TF", "target", "interaction_type", "ids"]
    trrust = trrust[["TF", "target"]].drop_duplicates()

    chip_mat = pd.read_csv(b.CHIP_ENCODE_GT_PATH, sep="\t", index_col=0)
    chip = (
        chip_mat.stack()
        .loc[lambda s: s.eq(1)]
        .reset_index()[["level_1", "level_0"]]
        .rename(columns={"level_1": "TF", "level_0": "target"})
        .drop_duplicates()
    )

    out = {}
    for name, gt in {"TRRUST": trrust, "K562_ChIP": chip}.items():
        gt = gt.copy()
        gt["TF"] = gt["TF"].astype(str).str.strip()
        gt["target"] = gt["target"].astype(str).str.strip()
        gt = gt[gt["target"].isin(valid_targets)]
        out[name] = gt.drop_duplicates(["TF", "target"]).reset_index(drop=True)
    return out


def load_ranked_edges_for_v2(valid_targets: set[str], tfs: set[str]) -> dict[str, pd.DataFrame]:
    rows = []
    ranked_path = TABLE_DIR / "regulator_prioritization_method_specific_ranked_edges.csv"
    for chunk in pd.read_csv(ranked_path, chunksize=1_000_000):
        chunk["TF"] = chunk["TF"].astype(str).str.strip()
        sub = chunk[(chunk["method"].isin(ALL_METHOD_ORDER)) & (chunk["TF"].isin(tfs))].copy()
        if len(sub):
            sub["target"] = sub["target"].astype(str).str.strip()
            sub = sub[sub["target"].isin(valid_targets)]
            rows.append(sub[["method", "TF", "target", "EdgeWeight"]])
    edges = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["method", "TF", "target", "EdgeWeight"])

    yao_irf3 = fig5_paths.CAPSULE_INPUT_DIR / "perturb_de" / "yao_irf3" / "TFs_vs_non_targeting_DE_IRF3.csv"
    if yao_irf3.exists() and "IRF3" in tfs:
        tmp = pd.read_csv(yao_irf3)
        tmp = tmp.rename(columns={"Gene name": "target", "P-value": "p_value"})
        tmp["method"] = "Yao2023_THP1"
        tmp["TF"] = "IRF3"
        tmp["target"] = tmp["target"].astype(str).str.strip()
        tmp = tmp[tmp["target"].isin(valid_targets)]
        tmp["EdgeWeight"] = -np.log10(pd.to_numeric(tmp["p_value"], errors="coerce").clip(lower=np.nextafter(0, 1)))
        edges = pd.concat([edges, tmp[["method", "TF", "target", "EdgeWeight"]]], ignore_index=True)

    replogle_fos = fig5_paths.CAPSULE_INPUT_DIR / "perturb_de" / "replogle_fos" / "deg_tables" / "ENSG00000170345.findmarkers.all.txt"
    if replogle_fos.exists() and "FOS" in tfs:
        with open(b.EN2SYM_PATH, "rb") as f:
            en2sym = pickle.load(f)
        tmp = pd.read_csv(replogle_fos, sep="\t")
        target_col = tmp.columns[0]
        tmp = tmp.rename(columns={target_col: "target", "p_val": "p_value"})
        tmp["method"] = "Replogle2022_K562"
        tmp["TF"] = "FOS"
        tmp["target"] = tmp["target"].astype(str).str.strip().map(en2sym)
        tmp = tmp.dropna(subset=["target"])
        tmp = tmp[tmp["target"].isin(valid_targets)]
        tmp["EdgeWeight"] = -np.log10(pd.to_numeric(tmp["p_value"], errors="coerce").clip(lower=np.nextafter(0, 1)))
        edges = pd.concat([edges, tmp[["method", "TF", "target", "EdgeWeight"]]], ignore_index=True)

    edges["EdgeWeight"] = pd.to_numeric(edges["EdgeWeight"], errors="coerce")
    edges = edges.dropna(subset=["EdgeWeight"])
    methods = {}
    for method in ALL_METHOD_ORDER:
        sub = edges[edges["method"] == method].copy()
        sub = (
            sub.sort_values("EdgeWeight", ascending=False, kind="mergesort")
            .drop_duplicates(["TF", "target"], keep="first")
            .reset_index(drop=True)
        )
        methods[method] = sub[["TF", "target", "EdgeWeight"]]
    return methods


def build_v2_metrics(method_edges: dict[str, pd.DataFrame], gts: dict[str, pd.DataFrame], valid_targets: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    curve_rows = []
    availability_rows = []

    for benchmark, gt in gts.items():
        gt_tfs = set(gt["TF"].unique())
        for method in ALL_METHOD_ORDER:
            edges = method_edges[method]
            available_tfs = set(edges["TF"].unique())
            for tf in method_requested_tfs(benchmark, method):
                has_reference = tf in gt_tfs and len(gt[gt["TF"] == tf]) > 0
                has_prediction = tf in available_tfs
                availability_rows.append({
                    "scenario": SCENARIO,
                    "benchmark": benchmark,
                    "method": method,
                    "method_label": VIOLIN_METHOD_LABELS[method],
                    "TF": tf,
                    "requested": True,
                    "has_reference_positives": has_reference,
                    "has_method_predictions": has_prediction,
                    "included_in_method_summary": has_reference and has_prediction,
                    "method_context": METHOD_CONTEXT[method],
                })
                if not (has_reference and has_prediction):
                    continue
                row, curve = b.evaluate_one(edges, gt, {tf}, valid_targets, curve_points=500)
                row["EPR_raw"] = row["EPR"]
                row["EPR"] = min(float(row["EPR"]), EPR_CAP) if pd.notna(row["EPR"]) else row["EPR"]
                row.update({
                    "scenario": SCENARIO,
                    "benchmark": benchmark,
                    "strategy": "method_specific_requested_tfs",
                    "scope": "per_tf",
                    "TF": tf,
                    "method": method,
                    "method_label": VIOLIN_METHOD_LABELS[method],
                    "method_context": METHOD_CONTEXT[method],
                    "n_tfs_eval": 1,
                    "n_targets_eval": len(valid_targets),
                })
                metric_rows.append(row)
                c = curve.copy()
                c["scenario"] = SCENARIO
                c["benchmark"] = benchmark
                c["scope"] = "per_tf"
                c["TF"] = tf
                c["method"] = method
                curve_rows.append(c)

    metrics = pd.DataFrame(metric_rows)
    curves = pd.concat(curve_rows, ignore_index=True) if curve_rows else pd.DataFrame()
    availability = pd.DataFrame(availability_rows)
    return metrics, curves, availability


def add_summary_tables(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, sub in metrics.groupby(["scenario", "benchmark", "method", "method_label"], dropna=False):
        rows.append({
            "scenario": keys[0],
            "benchmark": keys[1],
            "method": keys[2],
            "method_label": keys[3],
            "n_TFs": sub["TF"].nunique(),
            "TFs": ",".join(sorted(sub["TF"].unique())),
            "median_EPR": sub["EPR"].median(),
            "median_AUPRC_ratio": sub["AUPRC_ratio"].median(),
            "mean_EPR": sub["EPR"].mean(),
            "mean_AUPRC_ratio": sub["AUPRC_ratio"].mean(),
        })
    return pd.DataFrame(rows)


def publication_metric_matrix(summary: pd.DataFrame) -> pd.DataFrame:
    cols = []
    for benchmark in ["TRRUST", "K562_ChIP"]:
        sub = summary[summary["benchmark"] == benchmark].set_index("method")
        cols.append(sub["median_EPR"].reindex(ALL_METHOD_ORDER).rename((BENCHMARK_LABELS[benchmark], "EPR")))
        cols.append(sub["median_AUPRC_ratio"].reindex(ALL_METHOD_ORDER).rename((BENCHMARK_LABELS[benchmark], "AUPRCr")))
    mat = pd.concat(cols, axis=1).rename(index={m: VIOLIN_METHOD_LABELS[m] for m in ALL_METHOD_ORDER})
    mat.columns = pd.MultiIndex.from_tuples(mat.columns)
    return mat


def plot_publication_heatmap(summary: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.9)
    mat = publication_metric_matrix(summary)
    mat.to_csv(V2_TABLE_DIR / "system_broad_filtered_publication_v2_benchmark_metrics.csv")
    color = normalize_over_method_rows(mat)

    fig, ax = plt.subplots(1, 1, figsize=(6.75, 5.0))
    sns.heatmap(
        color,
        annot=format_annotations(mat),
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
    add_compact_row_label_bands(ax, [str(x) for x in mat.index])
    add_cell_type_row_groups(ax)
    style_heatmap_axis(ax, "")
    ax.set_xticklabels(mat.columns.get_level_values(1), rotation=0, fontsize=12)
    for group in list(dict.fromkeys(mat.columns.get_level_values(0))):
        idxs = [i for i, val in enumerate(mat.columns.get_level_values(0)) if val == group]
        center = (min(idxs) + max(idxs) + 1) / 2
        ax.text(center / mat.shape[1], 1.22, group, transform=ax.transAxes, ha="center", va="bottom", fontsize=14, fontweight="bold", color=SYSTEM_BLUE)
    for x in range(2, mat.shape[1], 2):
        ax.axvline(x, color="white", lw=5.0)
    add_low_high_legend(fig, x=0.31, y=0.075, width=0.38, height=0.025)
    fig.subplots_adjust(left=0.25, right=0.98, top=0.79, bottom=0.29)
    save(fig, "publication_heatmap_system_broad_filtered_EPR_AUPRCr", bbox_inches=None)


def per_tf_metric_matrix(metrics: pd.DataFrame, benchmark: str) -> pd.DataFrame:
    tf_order = ordered_tfs_for_benchmark(benchmark)
    cols = []
    sub = metrics[metrics["benchmark"] == benchmark]
    for tf in tf_order:
        tf_sub = sub[sub["TF"] == tf].set_index("method")
        cols.append(tf_sub["EPR"].reindex(ALL_METHOD_ORDER).rename((tf, "EPR")))
        cols.append(tf_sub["AUPRC_ratio"].reindex(ALL_METHOD_ORDER).rename((tf, "AUPRCr")))
    mat = pd.concat(cols, axis=1).rename(index={m: VIOLIN_METHOD_LABELS[m] for m in ALL_METHOD_ORDER})
    mat.columns = pd.MultiIndex.from_tuples(mat.columns)
    return mat


def plot_expanded_tf_heatmap(metrics: pd.DataFrame, benchmark: str) -> None:
    mat = per_tf_metric_matrix(metrics, benchmark)
    bench_slug = benchmark.lower()
    mat.to_csv(V2_TABLE_DIR / f"system_broad_filtered_publication_v2_{bench_slug}_tf_metrics.csv")
    color = normalize_over_method_rows(mat)
    fig, ax = plt.subplots(figsize=(12.5, 4.0))
    sns.heatmap(color, annot=format_annotations(mat), fmt="", cmap=MILD_VIRIDIS, vmin=0, vmax=1, cbar=False, linewidths=1.25, linecolor="white", ax=ax, annot_kws={"fontsize": 8.5})
    apply_annotation_colors(ax, color)
    add_compact_row_label_bands(ax, [str(x) for x in mat.index], width=0.18, fontsize=10.5)
    add_trrust_tf_cell_type_row_groups(ax)
    style_heatmap_axis(ax, "")
    ax.set_xticklabels(mat.columns.get_level_values(1), rotation=0, fontsize=9)
    for tf in list(dict.fromkeys(mat.columns.get_level_values(0))):
        idxs = [i for i, val in enumerate(mat.columns.get_level_values(0)) if val == tf]
        center = (min(idxs) + max(idxs) + 1) / 2
        ax.text(center / mat.shape[1], 1.20, tf, transform=ax.transAxes, ha="center", va="bottom", fontsize=10.5, fontweight="bold", color=SYSTEM_BLUE)
    for x in range(2, mat.shape[1], 2):
        ax.axvline(x, color="white", lw=4.0)
    fig.subplots_adjust(left=0.18, right=0.99, top=0.62, bottom=0.08)
    stems = [f"publication_heatmap_{benchmark}_v2_TFs_EPR_AUPRCr"]
    if benchmark == "TRRUST":
        stems.append("publication_heatmap_TRRUST_RELA_NFKB1_JUN_EPR_AUPRCr")
    save_existing_figure(fig, stems, bbox_inches=None)


def violin_source(metrics: pd.DataFrame) -> pd.DataFrame:
    out = metrics[["benchmark", "method", "method_label", "TF", "EPR", "AUPRC_ratio"]].copy()
    out["benchmark_label"] = out["benchmark"].map(BENCHMARK_LABELS)
    out.to_csv(V2_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_violin_metrics.csv", index=False)
    return out


def plot_violin_benchmark(metrics: pd.DataFrame, benchmark: str) -> None:
    data = violin_source(metrics)
    data = data[data["benchmark"] == benchmark].copy()
    method_labels = [VIOLIN_METHOD_LABELS[m] for m in ALL_METHOD_ORDER]
    tf_order = ordered_tfs_for_benchmark(benchmark)
    tf_offsets = {tf: offset for tf, offset in zip(tf_order, np.linspace(-0.32, 0.32, len(tf_order)))}
    palette = {VIOLIN_METHOD_LABELS[m]: b.METHOD_COLORS[m] for m in ALL_METHOD_ORDER}

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.8), sharey=True)
    for ax, metric, label in zip(axes, ["EPR", "AUPRC_ratio"], ["EPR", "AUPRCr"]):
        sns.violinplot(data=data, y="method_label", x=metric, hue="method_label", order=method_labels, hue_order=method_labels, palette=palette, inner=None, cut=0, linewidth=1.0, saturation=0.65, density_norm="width", dodge=False, legend=False, ax=ax)
        xmax = float(data[metric].replace([np.inf, -np.inf], np.nan).max())
        xmax = max(1.2, xmax * 1.25)
        ax.set_xlim(0, xmax)
        ax.axvline(1, color="#333333", lw=0.8, ls="--", alpha=0.45)
        for method_i, method in enumerate(ALL_METHOD_ORDER):
            for _, row in data[data["method"] == method].iterrows():
                y = method_i + tf_offsets[row["TF"]]
                x = float(row[metric])
                ax.scatter(x, y, s=16, color=b.METHOD_COLORS[method], edgecolor="white", linewidth=0.35, zorder=4)
                ax.text(x + xmax * 0.012, y, row["TF"], ha="left", va="center", fontsize=5.4, color="#222222", clip_on=True, zorder=5)
        ax.set_title(label, fontsize=12, pad=8)
        ax.set_xlabel(label)
        ax.set_ylabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_yticks(range(len(method_labels)))
    axes[0].set_yticklabels(method_labels, fontsize=9)
    fig.suptitle(BENCHMARK_LABELS[benchmark], fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, f"publication_violin_{benchmark}_per_tf_EPR_AUPRCr")


def compute_rank_coverage(method_edges: dict[str, pd.DataFrame], availability: pd.DataFrame, gts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for benchmark in ["TRRUST", "K562_ChIP"]:
        for method in ALL_METHOD_ORDER:
            relevant = availability[
                (availability["benchmark"] == benchmark)
                & (availability["method"] == method)
                & (availability["included_in_method_summary"])
            ]["TF"].tolist()
            if not relevant:
                continue
            pred = method_edges[method]
            pred = pred[pred["TF"].isin(relevant)].sort_values("EdgeWeight", ascending=False, kind="mergesort").reset_index(drop=True)
            benchmark_k = int(gts[benchmark][gts[benchmark]["TF"].isin(relevant)].drop_duplicates(["TF", "target"]).shape[0])
            requested_ks = {
                "Benchmark top-K": benchmark_k,
                "Top 5%": int(np.ceil(0.05 * len(pred))) if len(pred) else 0,
                "Top 1000": 1000,
            }
            for cutoff_label, requested_k in requested_ks.items():
                effective_k = min(requested_k, len(pred)) if requested_k > 0 else 0
                counts = pred.head(effective_k)["TF"].value_counts()
                for tf in relevant:
                    rows.append({
                        "benchmark": benchmark,
                        "benchmark_label": BENCHMARK_LABELS[benchmark],
                        "cutoff_label": cutoff_label,
                        "requested_k": requested_k,
                        "effective_emitted_k": effective_k,
                        "method": method,
                        "method_label": VIOLIN_METHOD_LABELS[method],
                        "TF": tf,
                        "rank_coverage": counts.get(tf, 0) / requested_k if requested_k > 0 else np.nan,
                        "n_reference_tfs": len(relevant),
                    })
    out = pd.DataFrame(rows)
    out.to_csv(V2_TABLE_DIR / "system_broad_filtered_publication_v2_tf_rank_coverage.csv", index=False)
    return out


def plot_tf_rank_coverage(source: pd.DataFrame, benchmark: str) -> None:
    data = source[source["benchmark"] == benchmark].copy()
    method_labels = [VIOLIN_METHOD_LABELS[m] for m in ALL_METHOD_ORDER]
    data["method_label"] = pd.Categorical(data["method_label"], method_labels, ordered=True)
    data["cutoff_label"] = pd.Categorical(data["cutoff_label"], RANK_COVERAGE_CUTOFFS, ordered=True)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 5.0), sharey=True)
    for ax, cutoff_label in zip(axes, RANK_COVERAGE_CUTOFFS):
        panel = data[data["cutoff_label"] == cutoff_label]
        xmax = panel["rank_coverage"].replace([np.inf, -np.inf], np.nan).max()
        xmax = 0.25 if not np.isfinite(xmax) else max(0.25, float(xmax) * 1.22)
        for method_i, method in enumerate(ALL_METHOD_ORDER):
            method_data = panel[panel["method"] == method]
            if method_data.empty:
                continue
            tf_order = method_data["TF"].tolist()
            offsets = {tf: offset for tf, offset in zip(tf_order, np.linspace(-0.34, 0.34, len(tf_order)))}
            for _, row in method_data.iterrows():
                y = method_i + offsets[row["TF"]]
                x = float(row["rank_coverage"])
                ax.scatter(x, y, s=28, marker="o", color=b.METHOD_COLORS[method], edgecolor="white", linewidth=0.5, alpha=0.98, zorder=4)
                ax.text(x + xmax * 0.015, y, row["TF"], ha="left", va="center", fontsize=5.3, color="#222222", clip_on=True, zorder=5)
        ax.set_xlim(0, xmax)
        ax.set_ylim(len(method_labels) - 0.5, -0.5)
        ax.set_title(cutoff_label, fontsize=12, pad=8)
        ax.set_xlabel("Fraction of top-ranked edges")
        ax.set_ylabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_yticks(range(len(method_labels)))
    axes[0].set_yticklabels(method_labels, fontsize=9)
    fig.suptitle(BENCHMARK_LABELS[benchmark], fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, f"publication_tf_rank_coverage_{benchmark}")


def plot_tf_pr_curves(per_tf_curves: pd.DataFrame, benchmark: str) -> None:
    for tf in ordered_tfs_for_benchmark(benchmark):
        panel = per_tf_curves[(per_tf_curves["benchmark"] == benchmark) & (per_tf_curves["TF"] == tf) & (per_tf_curves["method"].isin(PR_METHOD_ORDER))]
        if panel.empty:
            continue
        fig, ax = plt.subplots(figsize=(5.0, 4.0))
        for method in PR_METHOD_ORDER:
            d = panel[panel["method"] == method].sort_values("k")
            if d.empty:
                continue
            observed = d[~d["is_expectation"]].dropna(subset=["precision"])
            expected = d[d["is_expectation"]].dropna(subset=["precision"])
            color = b.METHOD_COLORS[method]
            if len(observed):
                ax.plot(observed["recall"], observed["precision"], drawstyle="steps-post", lw=2.2, color=color, label=b.METHOD_LABELS[method])
            if len(expected):
                start = observed.tail(1) if len(observed) else expected.head(1)
                tail = pd.concat([start, expected]).drop_duplicates(["recall", "precision"])
                ax.plot(tail["recall"], tail["precision"], lw=1.9, ls="--", color=color, alpha=0.75)
        ax.set_xlim(0, 0.25)
        ax.set_ylim(0, min(1, max(0.05, panel["precision"].replace([np.inf, -np.inf], np.nan).quantile(0.995) * 1.1)))
        ax.set_box_aspect(1)
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_title(tf, fontsize=12, pad=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if tf == "RELA":
            ax.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.6, borderaxespad=0.2)
        fig.tight_layout()
        save(fig, f"publication_pr_{benchmark}_{tf}")


def plot_macro_pr_for_benchmark(curves: pd.DataFrame, benchmark: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    recall_grid = np.linspace(0, 1, 250)
    for method in PR_METHOD_ORDER:
        method_curves = []
        panel = curves[(curves["benchmark"] == benchmark) & (curves["method"] == method)]
        for tf in sorted(panel["TF"].unique()):
            d = panel[panel["TF"] == tf].sort_values(["recall", "is_expectation"])
            if d.empty:
                continue
            d = d.dropna(subset=["recall", "precision"]).drop_duplicates("recall", keep="last")
            x = d["recall"].to_numpy(float)
            y = d["precision"].to_numpy(float)
            if len(x) == 0:
                continue
            order = np.argsort(x)
            method_curves.append(np.interp(recall_grid, x[order], y[order], left=y[order][0], right=y[order][-1]))
        if not method_curves:
            continue
        mean_precision = np.nanmean(np.vstack(method_curves), axis=0)
        ax.plot(recall_grid, mean_precision, lw=2.4, color=b.METHOD_COLORS[method], label=b.METHOD_LABELS[method])
    ax.set_xlim(0, 0.25 if benchmark == "TRRUST" else 1)
    ax.set_ylim(0, 0.5 if benchmark == "TRRUST" else 1)
    ax.set_box_aspect(1)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(BENCHMARK_LABELS[benchmark], fontsize=12, pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if benchmark == "TRRUST":
        ax.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.6, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, f"publication_mean_macro_pr_{benchmark}")


def main() -> None:
    V2_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    valid_targets = b.load_valid_targets()
    gts = load_reference_ground_truths(valid_targets)
    all_requested_tfs = set()
    for benchmark in ["TRRUST", "K562_ChIP"]:
        for method in ALL_METHOD_ORDER:
            all_requested_tfs.update(method_requested_tfs(benchmark, method))

    method_edges = load_ranked_edges_for_v2(valid_targets, all_requested_tfs)
    metrics, curves, availability = build_v2_metrics(method_edges, gts, valid_targets)
    summary = add_summary_tables(metrics)

    metrics.to_csv(V2_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_metrics.csv", index=False)
    curves.to_csv(V2_TABLE_DIR / "system_broad_filtered_publication_v2_per_tf_pr_curves.csv", index=False)
    availability.to_csv(V2_TABLE_DIR / "system_broad_filtered_publication_v2_tf_availability.csv", index=False)
    summary.to_csv(V2_TABLE_DIR / "system_broad_filtered_publication_v2_summary_by_method.csv", index=False)

    plot_publication_heatmap(summary)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_expanded_tf_heatmap(metrics, benchmark)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_violin_benchmark(metrics, benchmark)
    rank_coverage = compute_rank_coverage(method_edges, availability, gts)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_tf_rank_coverage(rank_coverage, benchmark)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_tf_pr_curves(curves, benchmark)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_macro_pr_for_benchmark(curves, benchmark)

    print(f"Wrote legacy publication figures to {OUT_DIR}")
    print(f"Wrote publication tables to {V2_TABLE_DIR}")


if __name__ == "__main__":
    main()
