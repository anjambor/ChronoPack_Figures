from __future__ import annotations

import glob
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns

import benchmark_grn_systems as b
import fig5_paths
from generate_system_broad_filtered_assets import (
    BENCHMARK_LABELS,
    OUT_DIR,
    SCENARIO,
    TABLE_DIR,
    add_low_high_legend,
    apply_annotation_colors,
    format_annotations,
    normalize_over_method_rows,
    style_heatmap_axis,
    CMAP,
    SYSTEM_BLUE,
)


SOURCE_TABLE_OUT = TABLE_DIR / "system_broad_filtered_publication_benchmark_metrics.csv"
TRRUST_TF_TABLE_OUT = TABLE_DIR / "system_broad_filtered_publication_trrust_core_tf_metrics.csv"
VIOLIN_TABLE_OUT = TABLE_DIR / "system_broad_filtered_publication_per_tf_violin_metrics.csv"
TF_RANK_COVERAGE_OUT = TABLE_DIR / "system_broad_filtered_publication_tf_rank_coverage.csv"
PUB_OUT_DIR = OUT_DIR / "publication_figures"
ALL_METHOD_ORDER = b.METHOD_ORDER
ALL_METHOD_SHORT = {
    "ChronoPack_K562": "Chrono",
    "GENIE3": "GENIE3",
    "SINCERITIES": "SINCERITIES",
    "GRNBOOST2": "GRNBOOST2",
    "PPCOR": "PPCOR",
    "Replogle2022_K562": "Replogle",
    "ChronoPack_THP1": "Chrono",
    "Yao2023_THP1": "Yao",
    **fig5_paths.PRIOR_LABELS,
}
VIOLIN_METHOD_LABELS = {
    "ChronoPack_K562": "Chrono\nK562",
    "GENIE3": "GENIE3\nK562",
    "SINCERITIES": "SINCERITIES\nK562",
    "GRNBOOST2": "GRNBOOST2\nK562",
    "PPCOR": "PPCOR\nK562",
    "Replogle2022_K562": "Replogle\nK562",
    "ChronoPack_THP1": "Chrono\nTHP1",
    "Yao2023_THP1": "Yao\nTHP1",
    **{m: f"{label}\nK562" for m, label in fig5_paths.PRIOR_LABELS.items()},
}
PR_METHOD_ORDER = [
    "ChronoPack_K562",
    "GENIE3",
    "Replogle2022_K562",
    "ChronoPack_THP1",
    "Yao2023_THP1",
]
# Fig5_GRN_v3: PR panels show GENIE3 for the BEELINE methods; with_prior adds the motif-prior methods after it
PR_METHOD_ORDER = fig5_paths.insert_prior(PR_METHOD_ORDER, after="GENIE3")
PR_METHOD_LABELS = {
    "ChronoPack_K562": "Chrono",
    "GENIE3": "GENIE3",
    "Replogle2022_K562": "Replogle",
    "ChronoPack_THP1": "Chrono",
    "Yao2023_THP1": "Yao",
    **fig5_paths.PRIOR_LABELS,
}
PR_CURVE_LABELS = {
    "ChronoPack_K562": "Chrono (K562)",
    "GENIE3": "GENIE3 (K562)",
    "Replogle2022_K562": "Replogle (K562)",
    "ChronoPack_THP1": "Chrono (THP1)",
    "Yao2023_THP1": "Yao (THP1)",
    **{m: f"{label} (K562)" for m, label in fig5_paths.PRIOR_LABELS.items()},
}
TRRUST_CORE_TFS = ["RELA", "NFKB1", "JUN"]
RANK_COVERAGE_CUTOFFS = ["Benchmark top-K", "Top 5%", "Top 1000"]
MILD_VIRIDIS = LinearSegmentedColormap.from_list(
    "mild_viridis",
    plt.get_cmap("viridis")(np.linspace(0.0, 0.90, 256)),
)


def save(fig: plt.Figure, stem: str, bbox_inches: str | None = "tight") -> None:
    PUB_OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PUB_OUT_DIR / f"{stem}.png", dpi=300, bbox_inches=bbox_inches)
    fig.savefig(PUB_OUT_DIR / f"{stem}.pdf", dpi=300, bbox_inches=bbox_inches)
    plt.close(fig)


def filtered_tfs(tf_sets: pd.DataFrame, benchmark: str) -> list[str]:
    row = tf_sets[tf_sets["benchmark"] == benchmark].iloc[0]
    return [tf for tf in str(row["TFs"]).split(",") if tf]


def publication_metric_matrix(metrics: pd.DataFrame, tf_sets: pd.DataFrame) -> pd.DataFrame:
    cols = []
    for benchmark in ["TRRUST", "K562_ChIP"]:
        tfs = filtered_tfs(tf_sets, benchmark)
        sub = metrics[
            (metrics["scenario"] == "system_broad_all")
            & (metrics["benchmark"] == benchmark)
            & (metrics["scope"] == "per_tf")
            & (metrics["method"].isin(ALL_METHOD_ORDER))
            & (metrics["TF"].isin(tfs))
        ]
        med = sub.groupby("method")[["EPR", "AUPRC_ratio"]].median()
        cols.append(med["EPR"].reindex(ALL_METHOD_ORDER).rename((BENCHMARK_LABELS[benchmark], "EPR")))
        cols.append(med["AUPRC_ratio"].reindex(ALL_METHOD_ORDER).rename((BENCHMARK_LABELS[benchmark], "AUPRCr")))
    mat = pd.concat(cols, axis=1).rename(index=ALL_METHOD_SHORT)
    mat.columns = pd.MultiIndex.from_tuples(mat.columns)
    return mat


def plot_publication_heatmap(metrics: pd.DataFrame, tf_sets: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.9)
    mat = publication_metric_matrix(metrics, tf_sets)
    mat.to_csv(SOURCE_TABLE_OUT)

    fig, ax = plt.subplots(1, 1, figsize=(6.75, 5.0))
    plot_compact_grouped_panel(ax, mat, group_label_y=1.22)
    add_cell_type_row_groups(ax)
    add_low_high_legend(fig, x=0.31, y=0.075, width=0.38, height=0.025)
    fig.subplots_adjust(left=0.25, right=0.98, top=0.79, bottom=0.29)
    save(fig, "publication_heatmap_system_broad_filtered_EPR_AUPRCr", bbox_inches=None)


def plot_compact_grouped_panel(ax, raw: pd.DataFrame, group_label_y: float = 1.22) -> None:
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
    add_compact_row_label_bands(ax, [str(x) for x in raw.index])
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


def add_compact_row_label_bands(ax, labels: list[str], width: float = 0.24, fontsize: float = 12) -> None:
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    for i, label in enumerate(labels):
        if i % 2 == 0:
            ax.add_patch(
                plt.Rectangle(
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
            -0.018,
            i + 0.5,
            label,
            transform=trans,
            ha="right",
            va="center",
            fontsize=fontsize,
            color=SYSTEM_BLUE,
            clip_on=False,
        )


def add_cell_type_row_groups(ax) -> None:
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    n_k562 = len([m for m in b.METHOD_ORDER if m not in b.THP1_METHODS])  # Fig5_GRN_v3: 6 in v2
    for label, y0, y1 in [("K562", 0, n_k562), ("THP1", n_k562, len(b.METHOD_ORDER))]:
        ax.text(
            -0.31,
            (y0 + y1) / 2,
            label,
            transform=trans,
            ha="center",
            va="center",
            rotation=90,
            fontsize=13,
            fontweight="bold",
            color="#2d4056",
            clip_on=False,
        )
        ax.plot([-0.255, -0.255], [y0 + 0.08, y1 - 0.08], transform=trans, color="#2d4056", lw=1.5, clip_on=False)
        ax.plot([-0.275, -0.255], [y0 + 0.08, y0 + 0.08], transform=trans, color="#2d4056", lw=1.5, clip_on=False)
        ax.plot([-0.275, -0.255], [y1 - 0.08, y1 - 0.08], transform=trans, color="#2d4056", lw=1.5, clip_on=False)


def trrust_core_tf_matrix(metrics: pd.DataFrame) -> pd.DataFrame:
    sub = metrics[
        (metrics["scenario"] == SCENARIO)
        & (metrics["benchmark"] == "TRRUST")
        & (metrics["scope"] == "per_tf")
        & (metrics["method"].isin(PR_METHOD_ORDER))
        & (metrics["TF"].isin(TRRUST_CORE_TFS))
    ]
    cols = []
    for tf in TRRUST_CORE_TFS:
        tf_sub = sub[sub["TF"] == tf].set_index("method")
        cols.append(tf_sub["EPR"].reindex(PR_METHOD_ORDER).rename((tf, "EPR")))
        cols.append(tf_sub["AUPRC_ratio"].reindex(PR_METHOD_ORDER).rename((tf, "AUPRCr")))
    mat = pd.concat(cols, axis=1).rename(index=PR_METHOD_LABELS)
    mat.columns = pd.MultiIndex.from_tuples(mat.columns)
    return mat


def plot_trrust_core_tf_heatmap(metrics: pd.DataFrame) -> None:
    sns.set_theme(style="white", font_scale=0.9)
    mat = trrust_core_tf_matrix(metrics)
    mat.to_csv(TRRUST_TF_TABLE_OUT)
    color = normalize_over_method_rows(mat)

    fig, ax = plt.subplots(figsize=(8.0, 4.0))
    sns.heatmap(
        color,
        annot=format_annotations(mat),
        fmt="",
        cmap=MILD_VIRIDIS,
        vmin=0,
        vmax=1,
        cbar=False,
        linewidths=1.25,
        linecolor="white",
        ax=ax,
        annot_kws={"fontsize": 10.5},
    )
    apply_annotation_colors(ax, color)
    add_compact_row_label_bands(ax, [str(x) for x in mat.index], width=0.18, fontsize=10.5)
    add_trrust_tf_cell_type_row_groups(ax)
    style_heatmap_axis(ax, "")
    ax.set_xticklabels(mat.columns.get_level_values(1), rotation=0, fontsize=10.5)
    for tf in TRRUST_CORE_TFS:
        idxs = [i for i, val in enumerate(mat.columns.get_level_values(0)) if val == tf]
        center = (min(idxs) + max(idxs) + 1) / 2
        ax.text(
            center / mat.shape[1],
            1.20,
            tf,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
            color=SYSTEM_BLUE,
        )
    for x in range(2, mat.shape[1], 2):
        ax.axvline(x, color="white", lw=5.0)
    fig.subplots_adjust(left=0.28, right=0.98, top=0.62, bottom=0.08)
    save(fig, "publication_heatmap_TRRUST_RELA_NFKB1_JUN_EPR_AUPRCr", bbox_inches=None)


def add_trrust_tf_cell_type_row_groups(ax) -> None:
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    n_k562 = len([m for m in PR_METHOD_ORDER if m not in b.THP1_METHODS])  # Fig5_GRN_v3: 3 in v2
    for label, y0, y1 in [("K562", 0, n_k562), ("THP1", n_k562, len(PR_METHOD_ORDER))]:
        ax.text(
            -0.27,
            (y0 + y1) / 2,
            label,
            transform=trans,
            ha="center",
            va="center",
            rotation=90,
            fontsize=11.5,
            fontweight="bold",
            color="#2d4056",
            clip_on=False,
        )
        ax.plot([-0.215, -0.215], [y0 + 0.08, y1 - 0.08], transform=trans, color="#2d4056", lw=1.4, clip_on=False)
        ax.plot([-0.235, -0.215], [y0 + 0.08, y0 + 0.08], transform=trans, color="#2d4056", lw=1.4, clip_on=False)
        ax.plot([-0.235, -0.215], [y1 - 0.08, y1 - 0.08], transform=trans, color="#2d4056", lw=1.4, clip_on=False)


def violin_source(metrics: pd.DataFrame, tf_sets: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for benchmark in ["TRRUST", "K562_ChIP"]:
        tfs = filtered_tfs(tf_sets, benchmark)
        sub = metrics[
            (metrics["scenario"] == "system_broad_all")
            & (metrics["benchmark"] == benchmark)
            & (metrics["scope"] == "per_tf")
            & (metrics["method"].isin(ALL_METHOD_ORDER))
            & (metrics["TF"].isin(tfs))
        ].copy()
        sub["benchmark_label"] = BENCHMARK_LABELS[benchmark]
        sub["method_label"] = sub["method"].map(VIOLIN_METHOD_LABELS)
        rows.append(sub[["benchmark", "benchmark_label", "method", "method_label", "TF", "EPR", "AUPRC_ratio"]])
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(VIOLIN_TABLE_OUT, index=False)
    return out


def plot_violin_benchmark(metrics: pd.DataFrame, tf_sets: pd.DataFrame, benchmark: str) -> None:
    source = violin_source(metrics, tf_sets)
    data = source[source["benchmark"] == benchmark].copy()
    method_labels = [VIOLIN_METHOD_LABELS[m] for m in ALL_METHOD_ORDER]
    tf_order = filtered_tfs(tf_sets, benchmark)
    tf_offsets = {
        tf: offset
        for tf, offset in zip(tf_order, np.linspace(-0.30, 0.30, len(tf_order)))
    }
    palette = {VIOLIN_METHOD_LABELS[m]: b.METHOD_COLORS[m] for m in ALL_METHOD_ORDER}

    sns.set_theme(style="white", font_scale=0.9)
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.8), sharey=True)
    for ax, metric, label in zip(axes, ["EPR", "AUPRC_ratio"], ["EPR", "AUPRCr"]):
        sns.violinplot(
            data=data,
            y="method_label",
            x=metric,
            hue="method_label",
            order=method_labels,
            hue_order=method_labels,
            palette=palette,
            inner=None,
            cut=0,
            linewidth=1.0,
            saturation=0.65,
            density_norm="width",
            dodge=False,
            legend=False,
            ax=ax,
        )
        xmax = float(data[metric].replace([np.inf, -np.inf], np.nan).max())
        xmax = max(1.2, xmax * 1.25)
        ax.set_xlim(0, xmax)
        ax.axvline(1, color="#333333", lw=0.8, ls="--", alpha=0.45)
        for method_i, method in enumerate(ALL_METHOD_ORDER):
            method_data = data[data["method"] == method]
            for _, row in method_data.iterrows():
                y = method_i + tf_offsets[row["TF"]]
                x = float(row[metric])
                color = b.METHOD_COLORS[method]
                ax.scatter(x, y, s=16, color=color, edgecolor="white", linewidth=0.35, zorder=4)
                ax.text(
                    x + xmax * 0.012,
                    y,
                    row["TF"],
                    ha="left",
                    va="center",
                    fontsize=5.7,
                    color="#222222",
                    clip_on=True,
                    zorder=5,
                )
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


def benchmark_reference_ground_truths(valid_targets: set[str]) -> dict[str, pd.DataFrame]:
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


def normalize_edges_for_reference_tfs(
    df: pd.DataFrame,
    method: str,
    valid_targets: set[str],
    reference_tfs: set[str],
) -> pd.DataFrame:
    df = df[["TF", "target", "EdgeWeight"]].copy()
    df["TF"] = df["TF"].astype(str).str.strip()
    df["target"] = df["target"].astype(str).str.strip()
    df["EdgeWeight"] = pd.to_numeric(df["EdgeWeight"], errors="coerce")
    df = df.dropna(subset=["TF", "target", "EdgeWeight"])
    df = df[df["TF"].isin(reference_tfs)]
    df = df[df["target"].isin(valid_targets)]
    df = df[df["TF"] != df["target"]]
    df = (
        df.sort_values("EdgeWeight", ascending=False, kind="mergesort")
        .drop_duplicates(["TF", "target"], keep="first")
        .reset_index(drop=True)
    )
    df["method"] = method
    return df


def load_beeline_reference_tf_method(name: str, valid_targets: set[str], reference_tfs: set[str]) -> pd.DataFrame:
    if name == "PPCOR":
        path = b.BEELINE_RESULTS_DIR / name / "outFile.txt"
        score_col = "corVal"
    else:
        path = b.BEELINE_RESULTS_DIR / name / "rankedEdges.csv"
        score_col = "EdgeWeight"

    rows = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=1_000_000):
        chunk = chunk.rename(columns={"Gene1": "TF", "Gene2": "target"})
        chunk["TF"] = chunk["TF"].astype(str).str.strip()
        chunk["target"] = chunk["target"].astype(str).str.strip()
        chunk = chunk[chunk["TF"].isin(reference_tfs)]
        chunk = chunk[chunk["target"].isin(valid_targets)]
        chunk = chunk[chunk["TF"] != chunk["target"]]
        if len(chunk):
            if name == "PPCOR":
                chunk["EdgeWeight"] = pd.to_numeric(chunk[score_col], errors="coerce").abs()
            rows.append(chunk[["TF", "target", "EdgeWeight"]])
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["TF", "target", "EdgeWeight"])
    return normalize_edges_for_reference_tfs(df, name, valid_targets, reference_tfs)


def load_chronopack_reference_tf_method(
    name: str,
    asso_path: Path,
    omnipath_path: Path,
    valid_targets: set[str],
    reference_tfs: set[str],
) -> pd.DataFrame:
    asso = pd.read_csv(asso_path)
    gene_col = "X" if "X" in asso.columns else asso.columns[0]
    asso = asso.rename(columns={gene_col: "target"})
    asso["target"] = asso["target"].astype(str)
    asso = asso[asso["target"].isin(valid_targets)].copy()
    asso["gene_score"] = 1 - pd.to_numeric(asso["lfdr"], errors="coerce")
    asso = asso.dropna(subset=["gene_score"])

    omni = pd.read_csv(omnipath_path)
    omni["TF"] = omni["TF"].astype(str).str.strip()
    omni = omni[omni["TF"].isin(reference_tfs)].copy()
    omni["priority_act"] = pd.to_numeric(omni["priority_act"], errors="coerce")
    omni = omni[omni["priority_act"] > 0].dropna(subset=["priority_act"])
    omni["tf_score"] = omni["priority_act"] / omni["priority_act"].sum()

    scores = np.outer(omni["tf_score"].to_numpy(float), asso["gene_score"].to_numpy(float))
    pred = (
        pd.DataFrame(scores, index=omni["TF"], columns=asso["target"])
        .stack()
        .rename("EdgeWeight")
        .reset_index()
        .rename(columns={"level_0": "TF", "level_1": "target"})
    )
    pred = normalize_edges_for_reference_tfs(pred, name, valid_targets, reference_tfs)
    return b.break_ties(pred)


def load_replogle_reference_tf_method(valid_targets: set[str], reference_tfs: set[str]) -> pd.DataFrame:
    with open(b.EN2SYM_PATH, "rb") as f:
        en2sym = pickle.load(f)
    rows = []
    for fname in sorted(glob.glob(str(b.REPLOGLIE_DEG_DIR / "*"))):
        tf_ensg = Path(fname).name.split(".")[0]
        tf = en2sym.get(tf_ensg)
        if tf is None or tf not in reference_tfs:
            continue
        tmp = pd.read_csv(fname, sep="\t")
        target_col = tmp.columns[0]
        tmp["TF"] = tf
        tmp["target"] = tmp[target_col].astype(str).map(en2sym)
        tmp = tmp.dropna(subset=["target"])
        tmp["EdgeWeight"] = -np.log10(pd.to_numeric(tmp["p_val"], errors="coerce").clip(lower=np.nextafter(0, 1)))
        rows.append(tmp[["TF", "target", "EdgeWeight"]])
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["TF", "target", "EdgeWeight"])
    return normalize_edges_for_reference_tfs(df, "Replogle2022_K562", valid_targets, reference_tfs)


def load_yao_reference_tf_method(valid_targets: set[str], reference_tfs: set[str]) -> pd.DataFrame:
    df = pd.read_csv(b.YAO2023_DEGS_PATH)
    df = df.rename(columns={"Gene name": "target"})
    df["EdgeWeight"] = -np.log10(pd.to_numeric(df["P-value"], errors="coerce").clip(lower=np.nextafter(0, 1)))
    return normalize_edges_for_reference_tfs(df[["TF", "target", "EdgeWeight"]], "Yao2023_THP1", valid_targets, reference_tfs)


def load_reference_tf_methods(valid_targets: set[str], reference_tfs: set[str]) -> dict[str, pd.DataFrame]:
    methods = {
        name: load_beeline_reference_tf_method(name, valid_targets, reference_tfs)
        for name in ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]
    }
    methods["ChronoPack_K562"] = load_chronopack_reference_tf_method(
        "ChronoPack_K562", b.K562_TRADESEQ_PATH, b.K562_OMNIPATH_PATH, valid_targets, reference_tfs
    )
    methods["ChronoPack_THP1"] = load_chronopack_reference_tf_method(
        "ChronoPack_THP1", b.THP1_TRADESEQ_PATH, b.THP1_OMNIPATH_PATH, valid_targets, reference_tfs
    )
    methods["Replogle2022_K562"] = load_replogle_reference_tf_method(valid_targets, reference_tfs)
    methods["Yao2023_THP1"] = load_yao_reference_tf_method(valid_targets, reference_tfs)
    return {m: methods[m] for m in ALL_METHOD_ORDER if m in methods}


def rank_coverage_source(tf_sets: pd.DataFrame) -> pd.DataFrame:
    valid_targets = b.load_valid_targets()
    gts = benchmark_reference_ground_truths(valid_targets)
    all_reference_tfs = set().union(*(set(gt["TF"].unique()) for gt in gts.values()))
    methods = load_reference_tf_methods(valid_targets, all_reference_tfs)
    rows = []

    for benchmark in ["TRRUST", "K562_ChIP"]:
        relevant_tfs = set(filtered_tfs(tf_sets, benchmark))
        gt = gts[benchmark]
        benchmark_tfs = set(gt["TF"].unique())

        for method in ALL_METHOD_ORDER:
            if method not in methods:
                continue
            pred = methods[method]
            method_tfs = benchmark_tfs & set(pred["TF"].unique())
            pred = pred[pred["TF"].isin(method_tfs)].copy()
            pred = (
                pred.sort_values("EdgeWeight", ascending=False, kind="mergesort")
                .drop_duplicates(["TF", "target"], keep="first")
                .reset_index(drop=True)
            )
            benchmark_k = int(gt[gt["TF"].isin(method_tfs)].shape[0])
            requested_ks = {
                "Benchmark top-K": benchmark_k,
                "Top 5%": int(np.ceil(0.05 * len(pred))) if len(pred) else 0,
                "Top 1000": 1000,
            }
            for cutoff_label in RANK_COVERAGE_CUTOFFS:
                requested_k = requested_ks[cutoff_label]
                effective_k = min(requested_k, len(pred)) if requested_k > 0 else 0
                counts = pred.head(effective_k)["TF"].value_counts()
                for tf in sorted(method_tfs):
                    rows.append({
                        "benchmark": benchmark,
                        "benchmark_label": BENCHMARK_LABELS[benchmark],
                        "cutoff_label": cutoff_label,
                        "requested_k": requested_k,
                        "effective_emitted_k": effective_k,
                        "method": method,
                        "method_label": VIOLIN_METHOD_LABELS[method],
                        "TF": tf,
                        "tf_class": "Relevant" if tf in relevant_tfs else "Background",
                        "n_reference_tfs": len(benchmark_tfs),
                        "n_predicted_reference_tfs": len(method_tfs),
                        "rank_coverage": counts.get(tf, 0) / requested_k if requested_k > 0 else np.nan,
                    })

    out = pd.DataFrame(rows)
    out.to_csv(TF_RANK_COVERAGE_OUT, index=False)
    return out


def plot_tf_rank_coverage(source: pd.DataFrame, benchmark: str) -> None:
    data = source[source["benchmark"] == benchmark].copy()
    method_labels = [VIOLIN_METHOD_LABELS[m] for m in ALL_METHOD_ORDER]
    data["method_label"] = pd.Categorical(data["method_label"], method_labels, ordered=True)
    data["cutoff_label"] = pd.Categorical(data["cutoff_label"], RANK_COVERAGE_CUTOFFS, ordered=True)

    sns.set_theme(style="white", font_scale=0.9)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 5.0), sharey=True)
    for ax, cutoff_label in zip(axes, RANK_COVERAGE_CUTOFFS):
        panel = data[data["cutoff_label"] == cutoff_label]
        xmax = panel["rank_coverage"].replace([np.inf, -np.inf], np.nan).max()
        xmax = 0.25 if not np.isfinite(xmax) else max(0.25, float(xmax) * 1.22)

        for method_i, method in enumerate(ALL_METHOD_ORDER):
            method_data = panel[panel["method"] == method]
            if method_data.empty:
                continue
            tf_order = sorted(method_data["TF"].unique())
            offsets = {tf: offset for tf, offset in zip(tf_order, np.linspace(-0.34, 0.34, len(tf_order)))}
            expected_share = 1 / float(method_data["n_predicted_reference_tfs"].iloc[0])
            ax.plot(
                [expected_share, expected_share],
                [method_i - 0.38, method_i + 0.38],
                color="#333333",
                lw=0.65,
                ls="--",
                alpha=0.28,
                zorder=1,
            )
            background = method_data[method_data["tf_class"] != "Relevant"]
            if len(background):
                ax.scatter(
                    background["rank_coverage"],
                    [method_i + offsets[tf] for tf in background["TF"]],
                    s=9,
                    marker="o",
                    color="#8f969e",
                    edgecolor="none",
                    alpha=0.25,
                    zorder=2,
                )
            for _, row in method_data.iterrows():
                if row["tf_class"] != "Relevant":
                    continue
                y = method_i + offsets[row["TF"]]
                x = float(row["rank_coverage"])
                ax.scatter(
                    x,
                    y,
                    s=28,
                    marker="o",
                    color=b.METHOD_COLORS[method],
                    edgecolor="white",
                    linewidth=0.5,
                    alpha=0.98,
                    zorder=4,
                )
                ax.text(
                    x + xmax * 0.015,
                    y,
                    row["TF"],
                    ha="left",
                    va="center",
                    fontsize=5.5,
                    color="#222222",
                    clip_on=True,
                    zorder=5,
                )

        ax.set_xlim(0, xmax)
        ax.set_ylim(len(method_labels) - 0.5, -0.5)
        ax.set_title(cutoff_label, fontsize=12, pad=8)
        ax.set_xlabel("Fraction of top-ranked edges")
        ax.set_ylabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_yticks(range(len(method_labels)))
    axes[0].set_yticklabels(method_labels, fontsize=9)
    handles = [
        Line2D([0], [0], marker="o", color="none", label="Relevant TF", markerfacecolor="#555555", markeredgecolor="white", markersize=5.5),
        Line2D([0], [0], marker="o", color="none", label="Background TF", markerfacecolor="#8f969e", markeredgecolor="none", markersize=4.0, alpha=0.35),
        Line2D([0], [0], color="#333333", lw=0.8, ls="--", alpha=0.35, label="Equal predicted-TF share"),
    ]
    axes[-1].legend(handles=handles, frameon=False, fontsize=7, loc="lower right", handlelength=1.5)
    fig.suptitle(BENCHMARK_LABELS[benchmark], fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, f"publication_tf_rank_coverage_{benchmark}")


def plot_trrust_tf_pr_curve(per_tf_curves: pd.DataFrame, tf: str) -> None:
    sns.set_theme(style="white", font_scale=0.95)
    panel = per_tf_curves[
        (per_tf_curves["scenario"] == SCENARIO)
        & (per_tf_curves["benchmark"] == "TRRUST")
        & (per_tf_curves["scope"] == "per_tf")
        & (per_tf_curves["TF"] == tf)
        & (per_tf_curves["method"].isin(PR_METHOD_ORDER))
    ].copy()

    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    for method in PR_METHOD_ORDER:
        d = panel[panel["method"] == method].sort_values("k")
        if d.empty:
            continue
        observed = d[~d["is_expectation"]].dropna(subset=["precision"])
        expected = d[d["is_expectation"]].dropna(subset=["precision"])
        color = b.METHOD_COLORS[method]
        label = PR_CURVE_LABELS[method]
        if len(observed):
            ax.plot(observed["recall"], observed["precision"], drawstyle="steps-post", lw=2.2, color=color, label=label)
        if len(expected):
            start = observed.tail(1) if len(observed) else expected.head(1)
            tail = pd.concat([start, expected]).drop_duplicates(["recall", "precision"])
            ax.plot(tail["recall"], tail["precision"], lw=1.9, ls="--", color=color, alpha=0.75)

    xmax = 0.10 if tf == "JUN" else 0.25
    ax.set_xlim(0, xmax)
    ymax = 0.40 if tf == "JUN" else 1
    ax.set_ylim(0, ymax)
    ax.set_box_aspect(1)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(tf, fontsize=12, pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if tf == "RELA":
        ax.legend(frameon=False, fontsize=7, loc="upper right", handlelength=1.6, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, f"publication_pr_TRRUST_{tf}")


def plot_mean_macro_pr_for_benchmark(curves: pd.DataFrame, benchmark: str) -> None:
    sns.set_theme(style="white", font_scale=0.95)
    panel = curves[
        (curves["scenario"] == SCENARIO)
        & (curves["benchmark"] == benchmark)
        & (curves["method"].isin(PR_METHOD_ORDER))
    ].copy()

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for method in PR_METHOD_ORDER:
        d = panel[panel["method"] == method].sort_values("recall")
        if d.empty:
            continue
        observed = d[d["segment"] == "observed_all_TFs"].dropna(subset=["precision"])
        expected = d[d["segment"] == "expected_completion"].dropna(subset=["precision"])
        color = b.METHOD_COLORS[method]
        label = b.METHOD_LABELS[method]
        if benchmark == "TRRUST" and method == "Yao2023_THP1":
            label = "Yao2023 | THP1 LPS perturb"
        if len(observed):
            ax.plot(
                observed["recall"],
                observed["precision"],
                lw=2.6,
                color=color,
                label=label,
            )
        if len(expected):
            ax.plot(
                expected["recall"],
                expected["precision"],
                lw=2.21,
                ls="--",
                color=color,
                alpha=0.75,
            )

    if benchmark == "TRRUST":
        xmin = 0
        xmax = 0.25
        ymin = 0
        ymax = 0.5
    else:
        xmin = 0
        xmax = 1
        ymin = 0.10
        ymax = panel["precision"].replace([np.inf, -np.inf], np.nan).dropna().quantile(0.995)
        ymax = 1 if not np.isfinite(ymax) else min(1, max(0.05, ymax * 1.12))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_box_aspect(1)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(BENCHMARK_LABELS[benchmark], fontsize=12, pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if benchmark == "TRRUST":
        ax.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.6, borderaxespad=0.2)
        fig.tight_layout()
    else:
        fig.tight_layout(rect=[0, 0, 1, 1])
    save(fig, f"publication_mean_macro_pr_{benchmark}")


def main() -> None:
    metrics = pd.read_csv(TABLE_DIR / "metrics_long.csv")
    filtered_metrics = pd.read_csv(TABLE_DIR / "metrics_five_method_positive_tfs.csv")
    tf_sets = pd.read_csv(TABLE_DIR / "five_method_positive_tf_sets.csv")
    curves = pd.read_csv(TABLE_DIR / "mean_macro_pr_observed_expected_positive_tfs.csv")
    per_tf_curves = pd.read_csv(TABLE_DIR / "system_broad_filtered_per_tf_pr_curves.csv")

    plot_publication_heatmap(metrics, tf_sets)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_violin_benchmark(metrics, tf_sets, benchmark)
    rank_coverage = rank_coverage_source(tf_sets)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_tf_rank_coverage(rank_coverage, benchmark)
    plot_trrust_core_tf_heatmap(filtered_metrics)
    for tf in TRRUST_CORE_TFS:
        plot_trrust_tf_pr_curve(per_tf_curves, tf)
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_mean_macro_pr_for_benchmark(curves, benchmark)

    print(f"Wrote {PUB_OUT_DIR / 'publication_heatmap_system_broad_filtered_EPR_AUPRCr.pdf'}")
    print(f"Wrote {PUB_OUT_DIR / 'publication_mean_macro_pr_TRRUST.pdf'}")
    print(f"Wrote {PUB_OUT_DIR / 'publication_mean_macro_pr_K562_ChIP.pdf'}")
    print(f"Wrote {PUB_OUT_DIR / 'publication_heatmap_TRRUST_RELA_NFKB1_JUN_EPR_AUPRCr.pdf'}")
    print(f"Wrote {PUB_OUT_DIR / 'publication_violin_TRRUST_per_tf_EPR_AUPRCr.pdf'}")
    print(f"Wrote {PUB_OUT_DIR / 'publication_violin_K562_ChIP_per_tf_EPR_AUPRCr.pdf'}")
    print(f"Wrote {PUB_OUT_DIR / 'publication_tf_rank_coverage_TRRUST.pdf'}")
    print(f"Wrote {PUB_OUT_DIR / 'publication_tf_rank_coverage_K562_ChIP.pdf'}")
    print(f"Wrote {SOURCE_TABLE_OUT}")
    print(f"Wrote {TRRUST_TF_TABLE_OUT}")
    print(f"Wrote {VIOLIN_TABLE_OUT}")
    print(f"Wrote {TF_RANK_COVERAGE_OUT}")


if __name__ == "__main__":
    main()
