# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

from collections import defaultdict
import pickle
from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import benchmark_grn_systems as b


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"
FIG_DIR = BASE_DIR / "benchmark_outputs" / "figures" / "regulator_prioritization_v2"
RANKED_EDGES_PATH = TABLE_DIR / "regulator_prioritization_method_specific_ranked_edges.csv"

METHOD_ORDER = [
    "ChronoPack_K562",
    "GENIE3",
    "SINCERITIES",
    "GRNBOOST2",
    "PPCOR",
    "ChronoPack_THP1",
]
METHOD_ORDER = fig5_paths.insert_prior(METHOD_ORDER)  # Fig5_GRN_v3

METHOD_LABELS = {
    "ChronoPack_K562": "Chrono\nK562",
    "GENIE3": "GENIE3",
    "SINCERITIES": "SINCERITIES",
    "GRNBOOST2": "GRNBOOST2",
    "PPCOR": "PPCOR",
    "ChronoPack_THP1": "Chrono\nTHP1",
    **fig5_paths.PRIOR_LABELS,
}

K562_EVAL_TFS = ["RELA", "NFKB1", "JUN", "ELK1", "CREB1", "ATF1", "FOS", "ATF2"]
THP1_EVAL_TFS = K562_EVAL_TFS + ["REL", "IRF3"]
TOP_KS = [1000, 5000, 10000]
SUPPLEMENTED_METHODS: set[str] = set()

TF_COLORS = {
    "RELA": "#D55E00",
    "NFKB1": "#0072B2",
    "JUN": "#009E73",
    "ELK1": "#8B5A2B",
    "CREB1": "#6A3D9A",
    "ATF1": "#E69F00",
    "FOS": "#CC79A7",
    "ATF2": "#56B4E9",
    "REL": "#A6CEE3",
    "IRF3": "#009E9E",
    "IRF7": "#F8E7BE",
    "IRF5": "#A5AA99",
    "STAT3": "#D9CBE6",
    "NFKB2": "#B3A2C7",
    "E2F1": "#D6B656",
    "Other": "#D4D4D4",
}


def requested_tfs(method: str) -> list[str]:
    return THP1_EVAL_TFS if method.endswith("_THP1") else K562_EVAL_TFS


def all_requested_tfs() -> list[str]:
    return list(dict.fromkeys(K562_EVAL_TFS + THP1_EVAL_TFS))


def enrichment_tfs(method: str) -> list[str]:
    if method == "ChronoPack_K562":
        return THP1_EVAL_TFS
    return requested_tfs(method)


def supplemental_perturb_edges() -> pd.DataFrame:
    rows = []

    yao_irf3 = fig5_paths.CAPSULE_INPUT_DIR / "perturb_de" / "yao_irf3" / "TFs_vs_non_targeting_DE_IRF3.csv"
    if yao_irf3.exists():
        tmp = pd.read_csv(yao_irf3)
        tmp = tmp.rename(columns={"Gene name": "target", "P-value": "p_value"})
        tmp["method"] = "Yao2023_THP1"
        tmp["TF"] = "IRF3"
        tmp["EdgeWeight"] = -np.log10(pd.to_numeric(tmp["p_value"], errors="coerce").clip(lower=np.nextafter(0, 1)))
        rows.append(tmp[["method", "TF", "target", "EdgeWeight"]])

    replogle_fos = fig5_paths.CAPSULE_INPUT_DIR / "perturb_de" / "replogle_fos" / "deg_tables" / "ENSG00000170345.findmarkers.all.txt"
    if replogle_fos.exists():
        with open(b.EN2SYM_PATH, "rb") as f:
            en2sym = pickle.load(f)
        tmp = pd.read_csv(replogle_fos, sep="\t")
        target_col = tmp.columns[0]
        tmp = tmp.rename(columns={target_col: "target", "p_val": "p_value"})
        tmp["method"] = "Replogle2022_K562"
        tmp["TF"] = "FOS"
        tmp["target"] = tmp["target"].astype(str).str.strip().map(en2sym)
        tmp = tmp.dropna(subset=["target"])
        tmp["EdgeWeight"] = -np.log10(pd.to_numeric(tmp["p_value"], errors="coerce").clip(lower=np.nextafter(0, 1)))
        rows.append(tmp[["method", "TF", "target", "EdgeWeight"]])

    if not rows:
        return pd.DataFrame(columns=["method", "TF", "target", "EdgeWeight"])
    out = pd.concat(rows, ignore_index=True)
    out = out.dropna(subset=["target"])
    out["target"] = out["target"].astype(str).str.strip()
    out["TF"] = out["TF"].astype(str).str.strip()
    out["EdgeWeight"] = pd.to_numeric(out["EdgeWeight"], errors="coerce")
    return out.dropna(subset=["target", "EdgeWeight"])


def load_supplemented_method_edges() -> dict[str, pd.DataFrame]:
    rows = []
    usecols = ["method", "TF", "target", "EdgeWeight"]
    for chunk in pd.read_csv(RANKED_EDGES_PATH, usecols=usecols, chunksize=1_000_000):
        sub = chunk[chunk["method"].isin(SUPPLEMENTED_METHODS)].copy()
        if len(sub):
            rows.append(sub)
    base = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=usecols)
    edges = pd.concat([base, supplemental_perturb_edges()], ignore_index=True)
    edges["TF"] = edges["TF"].astype(str).str.strip()
    edges["target"] = edges["target"].astype(str).str.strip()
    edges["EdgeWeight"] = pd.to_numeric(edges["EdgeWeight"], errors="coerce")
    edges = edges.dropna(subset=["EdgeWeight"])

    out = {}
    for method in SUPPLEMENTED_METHODS:
        sub = edges[edges["method"] == method].copy()
        sub = (
            sub.sort_values("EdgeWeight", ascending=False, kind="mergesort")
            .drop_duplicates(["TF", "target"], keep="first")
            .reset_index(drop=True)
        )
        sub["rank"] = np.arange(1, len(sub) + 1)
        sub["percentile"] = 1.0 - (sub["rank"] - 1) / max(len(sub), 1)
        out[method] = sub[["method", "TF", "target", "EdgeWeight", "rank", "percentile"]]
    return out


def ranked_edge_chunks(usecols: list[str]):
    for chunk in pd.read_csv(RANKED_EDGES_PATH, usecols=usecols, chunksize=1_000_000):
        chunk = chunk[chunk["method"].isin(METHOD_ORDER) & ~chunk["method"].isin(SUPPLEMENTED_METHODS)]
        if len(chunk):
            yield chunk
    if not SUPPLEMENTED_METHODS:
        return
    supplemented = load_supplemented_method_edges()
    for method in METHOD_ORDER:
        if method in supplemented:
            yield supplemented[method][usecols]


def load_topk_counts() -> pd.DataFrame:
    counts: dict[tuple[str, int], defaultdict[str, int]] = {
        (method, k): defaultdict(int) for method in METHOD_ORDER for k in TOP_KS
    }
    usecols = ["method", "TF", "rank"]
    for chunk in ranked_edge_chunks(usecols):
        chunk = chunk[chunk["method"].isin(METHOD_ORDER) & chunk["rank"].le(max(TOP_KS))]
        if chunk.empty:
            continue
        for k in TOP_KS:
            sub = chunk[chunk["rank"].le(k)]
            grouped = sub.groupby(["method", "TF"]).size()
            for (method, tf), n in grouped.items():
                counts[(method, k)][tf] += int(n)

    rows = []
    for method in METHOD_ORDER:
        for k in TOP_KS:
            total = sum(counts[(method, k)].values())
            for tf, count in counts[(method, k)].items():
                rows.append({
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "top_k": k,
                    "TF": tf,
                    "count": count,
                    "edge_share": count / total if total else 0.0,
                    "total_edges": total,
                    "is_evaluated_for_method": tf in requested_tfs(method),
                })
    return pd.DataFrame(rows)


def selected_tfs(data: pd.DataFrame) -> list[str]:
    selected = set(all_requested_tfs())
    for (_, _), sub in data.groupby(["method", "top_k"], sort=False):
        selected.update(sub.sort_values("edge_share", ascending=False).head(3)["TF"])
    max_share = data.groupby("TF")["edge_share"].max()
    return max_share.reindex(list(selected)).fillna(0).sort_values(ascending=False).index.tolist()


def evaluated_signal_share(data: pd.DataFrame, top_k: int = 10000) -> pd.Series:
    sub = data[data["top_k"] == top_k].copy()
    rows = {}
    for method in METHOD_ORDER:
        method_sub = sub[sub["method"] == method]
        rows[method] = method_sub[method_sub["TF"].isin(requested_tfs(method))]["edge_share"].sum()
    return pd.Series(rows).reindex(METHOD_ORDER)


def plot_stacked_topks(data: pd.DataFrame, keep_tfs: list[str]) -> None:
    order = keep_tfs + ["Other"]
    method_order = evaluated_signal_share(data, top_k=10000).sort_values(ascending=False).index.tolist()
    fig, axes = plt.subplots(1, len(TOP_KS), figsize=(18.2, 6.2), sharey=True)
    for ax, k in zip(axes, TOP_KS):
        sub = data[data["top_k"] == k].copy()
        kept = sub[sub["TF"].isin(keep_tfs)]
        other = (
            sub[~sub["TF"].isin(keep_tfs)]
            .groupby(["method", "method_label", "top_k", "total_edges"], as_index=False)
            .agg(count=("count", "sum"), edge_share=("edge_share", "sum"))
        )
        other["TF"] = "Other"
        plot = pd.concat([kept, other], ignore_index=True)
        pivot = (
            plot.pivot_table(index="method", columns="TF", values="edge_share", aggfunc="sum", fill_value=0)
            .reindex(index=method_order, columns=order, fill_value=0)
        )
        x = np.arange(len(method_order))
        bottom = np.zeros(len(method_order))
        for tf in order:
            vals = pivot[tf].to_numpy(float)
            color = TF_COLORS.get(tf, "#BDBDBD")
            bars = ax.bar(x, vals, bottom=bottom, width=0.72, color=color, edgecolor="white", linewidth=0.35)
            for i, (bar, val) in enumerate(zip(bars, vals)):
                method = method_order[i]
                should_label = tf != "Other" and val >= 0.025 and (tf in requested_tfs(method) or val >= 0.06)
                if should_label:
                    alpha_note = "" if tf in requested_tfs(method) else "*"
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bottom[i] + val / 2,
                        f"{tf}{alpha_note}",
                        ha="center",
                        va="center",
                        rotation=0,
                        fontsize=6.5,
                        color="white",
                    )
            bottom += vals
        ax.set_title(f"Top {k:,}")
        ax.set_xticks(x)
        ax.set_xticklabels([METHOD_LABELS[m] for m in method_order], rotation=0)
        ax.set_ylim(0, 1)
        ax.grid(axis="y", color="#e8e8e8")
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("TF proportion")
    fig.suptitle("Regulator Proportions in Fixed Top-K Predictions (v2 evaluated TFs)", y=1.02, fontsize=14)
    legend_tfs = [tf for tf in order if tf != "Other" and tf in TF_COLORS] + ["Other"]
    handles = [plt.Rectangle((0, 0), 1, 1, color=TF_COLORS[tf]) for tf in legend_tfs]
    fig.legend(handles, legend_tfs, frameon=False, bbox_to_anchor=(1.005, 0.92), loc="upper left", fontsize=8)
    fig.tight_layout(rect=(0, 0, 0.90, 1))
    fig.savefig(FIG_DIR / "fixed_topk_tf_proportions_stacked_v2.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "fixed_topk_tf_proportions_stacked_v2.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_top10000_evaluated_stacked(data: pd.DataFrame, keep_tfs: list[str]) -> None:
    top_k = 10000
    sub = data[data["top_k"] == top_k].copy()
    method_order = evaluated_signal_share(data, top_k=top_k).sort_values(ascending=False).index.tolist()
    keep_tfs = [tf for tf in keep_tfs if tf in all_requested_tfs()] + [
        tf for tf in keep_tfs if tf not in all_requested_tfs()
    ]
    order = keep_tfs + ["Other"]
    kept = sub[sub["TF"].isin(keep_tfs)]
    other = (
        sub[~sub["TF"].isin(keep_tfs)]
        .groupby(["method", "method_label", "top_k", "total_edges"], as_index=False)
        .agg(count=("count", "sum"), edge_share=("edge_share", "sum"))
    )
    other["TF"] = "Other"
    plot = pd.concat([kept, other], ignore_index=True)
    pivot = (
        plot.pivot_table(index="method", columns="TF", values="edge_share", aggfunc="sum", fill_value=0)
        .reindex(index=method_order, columns=order, fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(12.2, 6.4))
    x = np.arange(len(method_order))
    bottom = np.zeros(len(method_order))
    for tf in order:
        vals = pivot[tf].to_numpy(float)
        bars = ax.bar(x, vals, bottom=bottom, width=0.72, color=TF_COLORS.get(tf, "#BDBDBD"), edgecolor="white", linewidth=0.35)
        for i, (bar, val) in enumerate(zip(bars, vals)):
            method = method_order[i]
            if tf != "Other" and val >= 0.018 and (tf in requested_tfs(method) or val >= 0.06):
                ax.text(bar.get_x() + bar.get_width() / 2, bottom[i] + val / 2, tf, ha="center", va="center", fontsize=6.5, color="white")
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([METHOD_LABELS[m] for m in method_order])
    ax.set_ylim(0, 1)
    ax.set_ylabel(f"TF proportion in top {top_k:,} predictions")
    ax.set_title(f"Regulator Proportions in Top {top_k:,} Predictions (v2 evaluated TFs)")
    ax.grid(axis="y", color="#e8e8e8")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    legend_tfs = [tf for tf in order if tf != "Other" and tf in TF_COLORS] + ["Other"]
    handles = [plt.Rectangle((0, 0), 1, 1, color=TF_COLORS[tf]) for tf in legend_tfs]
    ax.legend(handles, legend_tfs, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "top10000_tf_proportions_stacked_v2.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "top10000_tf_proportions_stacked_v2.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_top10000_evaluated_only_stacked(data: pd.DataFrame) -> None:
    top_k = 10000
    sub = data[data["top_k"] == top_k].copy()
    method_order = evaluated_signal_share(data, top_k=top_k).sort_values(ascending=False).index.tolist()
    eval_tfs = all_requested_tfs()
    order = eval_tfs + ["Other"]

    sub["plot_TF"] = np.where(sub["is_evaluated_for_method"], sub["TF"], "Other")
    plot = (
        sub.groupby(["method", "method_label", "plot_TF"], as_index=False)
        .agg(count=("count", "sum"), edge_share=("edge_share", "sum"))
        .rename(columns={"plot_TF": "TF"})
    )
    pivot = (
        plot.pivot_table(index="method", columns="TF", values="edge_share", aggfunc="sum", fill_value=0)
        .reindex(index=method_order, columns=order, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(12.2, 6.4))
    x = np.arange(len(method_order)) * 0.616
    bottom = np.zeros(len(method_order))
    for tf in order:
        vals = pivot[tf].to_numpy(float)
        bars = ax.bar(
            x,
            vals,
            bottom=bottom,
            width=0.36,
            color=TF_COLORS.get(tf, "#BDBDBD") if tf != "Other" else TF_COLORS["Other"],
            edgecolor="white",
            linewidth=0.35,
        )
        for i, (bar, val) in enumerate(zip(bars, vals)):
            method = method_order[i]
            if tf != "Other" and tf in requested_tfs(method) and val >= 0.018:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bottom[i] + val / 2,
                    tf,
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color="white",
                )
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([METHOD_LABELS[m] for m in method_order])
    ax.set_ylim(0, 0.6)
    ax.set_ylabel(f"TF proportion in top {top_k:,} predictions")
    ax.set_title(f"Evaluated/Signaling Regulator Proportions in Top {top_k:,} Predictions")
    ax.grid(axis="y", color="#e8e8e8")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    legend_tfs = [tf for tf in order if tf != "Other" and pivot[tf].sum() > 0] + ["Other"]
    handles = [plt.Rectangle((0, 0), 1, 1, color=TF_COLORS[tf]) for tf in legend_tfs]
    ax.legend(handles, legend_tfs, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "top10000_tf_proportions_stacked_v2_evaluated_only.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "top10000_tf_proportions_stacked_v2_evaluated_only.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_top10000_evaluated_only_pies(data: pd.DataFrame) -> None:
    top_k = 10000
    sub = data[data["top_k"] == top_k].copy()
    method_order = evaluated_signal_share(data, top_k=top_k).sort_values(ascending=False).index.tolist()
    eval_tfs = all_requested_tfs()
    order = eval_tfs + ["Other"]

    sub["plot_TF"] = np.where(sub["is_evaluated_for_method"], sub["TF"], "Other")
    plot = (
        sub.groupby(["method", "method_label", "plot_TF"], as_index=False)
        .agg(count=("count", "sum"), edge_share=("edge_share", "sum"))
        .rename(columns={"plot_TF": "TF"})
    )
    pivot = (
        plot.pivot_table(index="method", columns="TF", values="edge_share", aggfunc="sum", fill_value=0)
        .reindex(index=method_order, columns=order, fill_value=0)
    )

    fig, axes = plt.subplots(2, 3, figsize=(10.6, 6.9))
    axes = axes.ravel()
    for ax, method in zip(axes, method_order):
        vals = pivot.loc[method]
        vals = vals[vals.gt(0)]
        labels = [tf if vals.loc[tf] >= (0.12 if tf == "Other" else 0.035) else "" for tf in vals.index]
        colors = [TF_COLORS.get(tf, TF_COLORS["Other"]) if tf != "Other" else TF_COLORS["Other"] for tf in vals.index]
        wedges, texts = ax.pie(
            vals.to_numpy(float),
            labels=labels,
            colors=colors,
            startangle=90,
            counterclock=False,
            labeldistance=0.58,
            wedgeprops={"linewidth": 0.45, "edgecolor": "white"},
            textprops={"fontsize": 7.2, "color": "white", "ha": "center", "va": "center"},
        )
        for text, label in zip(texts, labels):
            if label == "":
                text.set_visible(False)
            elif label == "Other":
                text.set_color("#555555")
        ax.set_title(METHOD_LABELS[method], fontsize=10, pad=4)
        ax.set_aspect("equal")

    for ax in axes[len(method_order):]:
        ax.set_visible(False)

    legend_tfs = [tf for tf in order if tf != "Other" and pivot[tf].sum() > 0] + ["Other"]
    handles = [plt.Rectangle((0, 0), 1, 1, color=TF_COLORS[tf]) for tf in legend_tfs]
    fig.legend(handles, legend_tfs, frameon=False, loc="center left", bbox_to_anchor=(0.91, 0.52), fontsize=8)
    fig.suptitle(f"Evaluated/Signaling Regulator Proportions in Top {top_k:,} Predictions", y=0.98, fontsize=13)
    fig.tight_layout(rect=(0, 0, 0.90, 0.95))
    fig.savefig(FIG_DIR / "top10000_tf_proportions_pies_v2_evaluated_only.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "top10000_tf_proportions_pies_v2_evaluated_only.pdf", bbox_inches="tight")
    plt.close(fig)


def running_es(hit_ranks: np.ndarray, n_total: int) -> tuple[float, int]:
    hit_ranks = np.sort(hit_ranks.astype(np.int64))
    n_hits = len(hit_ranks)
    n_misses = n_total - n_hits
    if n_hits == 0 or n_misses <= 0:
        return np.nan, -1
    hit_index = np.arange(1, n_hits + 1, dtype=float)
    misses_before = hit_ranks - hit_index
    running_at_hits = hit_index / n_hits - misses_before / n_misses
    peak_idx = int(np.argmax(running_at_hits))
    return float(running_at_hits[peak_idx]), int(hit_ranks[peak_idx])


def running_curve(hit_ranks: np.ndarray, n_total: int, points: int = 1000) -> pd.DataFrame:
    hit_ranks = np.sort(hit_ranks.astype(np.int64))
    n_hits = len(hit_ranks)
    n_misses = n_total - n_hits
    x = np.unique(np.rint(np.linspace(1, n_total, min(points, n_total))).astype(np.int64))
    if n_hits == 0 or n_misses <= 0:
        return pd.DataFrame({"rank": x, "rank_fraction": x / n_total, "running_es": np.nan})
    hits_seen = np.searchsorted(hit_ranks, x, side="right")
    misses_seen = x - hits_seen
    es = hits_seen / n_hits - misses_seen / n_misses
    return pd.DataFrame({"rank": x, "rank_fraction": x / n_total, "running_es": es})


def random_es_values(n_total: int, n_hits: int, rng: np.random.Generator, n_perm: int = 200) -> np.ndarray:
    if n_hits == 0 or n_hits >= n_total:
        return np.array([np.nan])
    values = np.empty(n_perm, dtype=float)
    hit_index = np.arange(1, n_hits + 1, dtype=float)
    n_misses = n_total - n_hits
    for i in range(n_perm):
        ranks = np.sort(rng.choice(n_total, size=n_hits, replace=False) + 1)
        misses_before = ranks - hit_index
        values[i] = np.max(hit_index / n_hits - misses_before / n_misses)
    return values


def compute_rank_enrichment(include_curves: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(20260814)
    totals = {method: 0 for method in METHOD_ORDER}
    hit_ranks: dict[tuple[str, str], list[int]] = {
        (method, tf): [] for method in METHOD_ORDER for tf in enrichment_tfs(method)
    }
    tfs_union = set(all_requested_tfs())
    usecols = ["method", "TF", "rank"]
    for chunk in ranked_edge_chunks(usecols):
        chunk = chunk[chunk["method"].isin(METHOD_ORDER)]
        if chunk.empty:
            continue
        for method, n in chunk.groupby("method", sort=False).size().items():
            totals[method] += int(n)
        sub = chunk[chunk["TF"].isin(tfs_union)]
        for (method, tf), ranks in sub.groupby(["method", "TF"], sort=False)["rank"]:
            if (method, tf) in hit_ranks:
                hit_ranks[(method, tf)].extend(ranks.astype(int).tolist())

    rows = []
    curve_rows = []
    for method in METHOD_ORDER:
        n_total = totals[method]
        for tf in enrichment_tfs(method):
            ranks = np.array(hit_ranks[(method, tf)], dtype=np.int64)
            es, peak_rank = running_es(ranks, n_total)
            null_values = random_es_values(n_total, len(ranks), rng)
            null_mean = float(np.nanmean(null_values))
            nes = es / null_mean if null_mean and not np.isnan(null_mean) else np.nan
            p_value = float((np.sum(null_values >= es) + 1) / (np.sum(~np.isnan(null_values)) + 1)) if np.any(~np.isnan(null_values)) else np.nan
            rows.append({
                "method": method,
                "method_label": METHOD_LABELS[method],
                "TF": tf,
                "n_edges": n_total,
                "n_tf_edges": len(ranks),
                "edge_fraction": len(ranks) / n_total if n_total else np.nan,
                "ES": es,
                "null_mean_ES": null_mean,
                "NES": nes,
                "empirical_p_value": p_value,
                "leading_edge_rank": peak_rank,
                "leading_edge_rank_fraction": peak_rank / n_total if peak_rank > 0 else np.nan,
            })
            if include_curves:
                curve = running_curve(ranks, n_total)
                curve["method"] = method
                curve["TF"] = tf
                curve_rows.append(curve)
    curves = pd.concat(curve_rows, ignore_index=True) if curve_rows else pd.DataFrame()
    return pd.DataFrame(rows), curves


def plot_enrichment_heatmap(metrics: pd.DataFrame, value: str, stem: str, title: str, cmap: str) -> None:
    tf_order = all_requested_tfs()
    pivot = (
        metrics.pivot_table(index="TF", columns="method", values=value, aggfunc="first")
        .reindex(index=tf_order, columns=METHOD_ORDER)
    )
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    im = ax.imshow(pivot.to_numpy(float), aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(METHOD_ORDER)))
    ax.set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER])
    ax.set_yticks(np.arange(len(tf_order)))
    ax.set_yticklabels(tf_order)
    for i, tf in enumerate(tf_order):
        for j, method in enumerate(METHOD_ORDER):
            val = pivot.loc[tf, method]
            ax.text(j, i, "" if pd.isna(val) else f"{val:.2f}", ha="center", va="center", fontsize=7.5)
    ax.set_title(title)
    ax.spines[:].set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cbar.set_label(value)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_running_curves(curves: pd.DataFrame, metrics: pd.DataFrame) -> None:
    n_cols = 4
    n_rows = int(np.ceil(len(METHOD_ORDER) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.8 * n_cols, 4.0 * n_rows), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, method in zip(axes, METHOD_ORDER):
        sub = curves[curves["method"] == method]
        metric_sub = metrics[metrics["method"] == method].set_index("TF")
        for tf in enrichment_tfs(method):
            curve = sub[sub["TF"] == tf]
            if curve.empty:
                continue
            color = TF_COLORS.get(tf, "#777777")
            ax.plot(curve["rank_fraction"], curve["running_es"], color=color, label=tf, linewidth=1.45)
            peak_x = metric_sub.loc[tf, "leading_edge_rank_fraction"]
            peak_y = metric_sub.loc[tf, "ES"]
            if pd.notna(peak_x) and pd.notna(peak_y):
                ax.scatter([peak_x], [peak_y], color=color, s=11, zorder=3)
        ax.axhline(0, color="#777777", linewidth=0.7)
        ax.set_title(METHOD_LABELS[method].replace("\n", " "))
        ax.grid(color="#e8e8e8", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for ax in axes[len(METHOD_ORDER):]:
        ax.set_visible(False)
    for ax in axes[-n_cols:]:
        ax.set_xlabel("Rank fraction")
    for ax in axes[::n_cols]:
        ax.set_ylabel("Running ES")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper center", ncol=min(len(labels), 10), bbox_to_anchor=(0.5, 1.01), fontsize=8)
    fig.suptitle("Running TF Enrichment Across Full Ranked Edge Lists (v2 evaluated TFs)", y=1.06, fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "tf_rank_enrichment_running_curves_v2.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "tf_rank_enrichment_running_curves_v2.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    data = load_topk_counts()
    data.to_csv(TABLE_DIR / "regulator_prioritization_v2_fixed_topk_tf_proportions_full.csv", index=False)
    keep_tfs = selected_tfs(data)
    plot_stacked_topks(data, keep_tfs)
    plot_top10000_evaluated_stacked(data, keep_tfs)
    plot_top10000_evaluated_only_stacked(data)

    metrics, curves = compute_rank_enrichment()
    metrics.to_csv(TABLE_DIR / "regulator_prioritization_v2_tf_rank_enrichment.csv", index=False)
    curves.to_csv(TABLE_DIR / "regulator_prioritization_v2_tf_rank_enrichment_curves.csv", index=False)
    plot_enrichment_heatmap(metrics, "NES", "tf_rank_enrichment_nes_heatmap_v2", "Rank-Based TF Enrichment NES (v2 evaluated TFs)", "viridis")
    plot_enrichment_heatmap(metrics, "leading_edge_rank_fraction", "tf_rank_enrichment_leading_edge_peak_ranks_v2", "Leading-edge Peak Rank Fraction (v2 evaluated TFs)", "magma_r")
    plot_running_curves(curves, metrics)
    print(f"Wrote v2 regulator prioritization figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
