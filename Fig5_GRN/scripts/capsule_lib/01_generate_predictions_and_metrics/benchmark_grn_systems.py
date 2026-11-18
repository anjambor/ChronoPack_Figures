# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

import glob
import math
import pickle
import re
from dataclasses import dataclass
from pathlib import Path

import fig5_paths

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


BASE_DIR = fig5_paths.BENCH_ROOT
OUT_DIR = BASE_DIR / "benchmark_outputs"
TABLE_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"

CHRONO_ANALYSES_DIR = (fig5_paths.EXTERNAL_DIR / "projects/chrono_seq/analyses")
K562_ANALYSIS_DIR = CHRONO_ANALYSES_DIR / "k562-tnfa.hg38.rep"
GRN_INFERENCE_DIR = (fig5_paths.EXTERNAL_DIR / "projects/grn_inference")

PROTEIN_CODING_GENES_PATH = (fig5_paths.EXTERNAL_DIR / "ref/genome/hg38/annotations/protein_coding_genes.csv")
TRRUST_PATH = (fig5_paths.EXTERNAL_DIR / "ref/grns/ttrust/trrust_rawdata.human.tsv")
CHIP_ENCODE_GT_PATH = (fig5_paths.EXTERNAL_DIR / "downloads/k562-chip_encode/results/tf_gene.1kbp_200bp.stringent.tsv")
BEELINE_RESULTS_DIR = fig5_paths.BEELINE_RESULTS_DIR  # Fig5_GRN_v3: r9 P < 0.05 TFs + 2000 run (v2: 10% filter; v1: metacells)
PRIOR_RESULTS_DIR = fig5_paths.PRIOR_RESULTS_DIR  # Fig5_GRN_v3: motif-prior networks (FIG5_METHOD_SET=with_prior)
K562_TRADESEQ_PATH = K562_ANALYSIS_DIR / "intmd/sc-K562-1e4_to_tradeSeq/asso_df.1trajectory.lfdr.csv"
K562_OMNIPATH_PATH = K562_ANALYSIS_DIR / "intmd/omnipath_tf_ranks/omnipath_ranks.TNFa.csv"
THP1_TRADESEQ_PATH = CHRONO_ANALYSES_DIR / "thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq/asso_df.1trajectory.lfdr.csv"
THP1_OMNIPATH_PATH = K562_ANALYSIS_DIR / "intmd/omnipath_tf_ranks/omnipath_ranks.TLR4.csv"
EN2SYM_PATH = GRN_INFERENCE_DIR / "causalbench/data_cache/en2sym.pkl"
REPLOGLIE_DEG_DIR = (fig5_paths.EXTERNAL_DIR / "downloads/perturb_db/seurat_perturbdb_workflow/results/deg_tables")
YAO2023_DEGS_PATH = GRN_INFERENCE_DIR / "causalbench/notebooks/Yao2023/TFs_vs_non_targeting_DE_all.csv"

CORE_TFS = ["RELA", "NFKB1", "JUN", "ATF2", "ELK1", "JUNB", "EGR1"]
DISPLAY_TFS = ["RELA", "NFKB1", "JUN"]

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
METHOD_ORDER = fig5_paths.insert_prior(METHOD_ORDER)  # Fig5_GRN_v3: + PYSCENIC, CELLORACLE, INFERELATOR (with_prior)

METHOD_COLORS = {
    "ChronoPack_K562": "#D99000",
    "GENIE3": "#3D8B37",
    "SINCERITIES": "#69A761",
    "GRNBOOST2": "#008F77",
    "PPCOR": "#6A9A9A",
    "Replogle2022_K562": "#CC79A7",
    "ChronoPack_THP1": "#D55E00",
    "Yao2023_THP1": "#3366AA",
    **fig5_paths.PRIOR_COLORS,
}

METHOD_LABELS = {
    "ChronoPack_K562": "ChronoPack | K562 TNF",
    "GENIE3": "GENIE3 | K562 TNF",
    "SINCERITIES": "SINCERITIES | K562 TNF",
    "GRNBOOST2": "GRNBOOST2 | K562 TNF",
    "PPCOR": "PPCOR | K562 TNF",
    "Replogle2022_K562": "Replogle2022 | K562 perturb",
    "ChronoPack_THP1": "ChronoPack | THP1 LPS",
    "Yao2023_THP1": "Yao2023 | THP1 perturb",
    **{m: f"{label} | K562 TNF" for m, label in fig5_paths.PRIOR_LABELS.items()},
}

K562_NATIVE_METHODS = fig5_paths.insert_prior(["ChronoPack_K562", "GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"])
K562_PLUS_REPLOGLE_METHODS = K562_NATIVE_METHODS + ["Replogle2022_K562"]
THP1_METHODS = ["ChronoPack_THP1", "Yao2023_THP1"]


def strip_rows(methods: list[str], height: float) -> tuple[list[float], float, tuple[float, float]]:
    """Fig5_GRN_v3: strip-forest rows for any number of methods.

    Returns the row y positions (v2's 0.86 spacing), the figure height and the axes-fraction y of the K562 and THP1
    block labels. With v2's 8 methods these are v2's exact values (rows 0-6.02, the given height, 0.62 / 0.16); with
    more rows the height grows in proportion and the labels are centred on their blocks.
    """
    y = [round(0.86 * i, 2) for i in range(len(methods))]
    if len(methods) == 8:
        return y, height, (0.62, 0.16)
    top, bottom = -0.62, y[-1] + 0.55
    def frac(rows):
        mid = (y[rows[0]] + y[rows[-1]]) / 2
        return (bottom - mid) / (bottom - top)
    k562 = [i for i, m in enumerate(methods) if m not in THP1_METHODS]
    thp1 = [i for i, m in enumerate(methods) if m in THP1_METHODS]
    return y, height * len(methods) / 8, (frac(k562), frac(thp1))


@dataclass(frozen=True)
class MethodData:
    name: str
    edges: pd.DataFrame
    native_tfs: frozenset[str]
    native_targets: frozenset[str]
    context: str


def ensure_dirs() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_valid_targets() -> set[str]:
    genes = pd.read_csv(PROTEIN_CODING_GENES_PATH, header=None)[0].dropna().astype(str)
    rp_re = re.compile(r"^RP[SL][0-9]")
    # FIG5_UNIVERSE=shared: cut to the genes every non-Perturb-seq method ranks (fig5_paths.restrict_universe)
    return fig5_paths.restrict_universe({
        g for g in genes
        if not (g.startswith("MT-") or g.startswith("ENSG") or rp_re.match(g))
    })


def load_ground_truths(valid_targets: set[str]) -> dict[str, pd.DataFrame]:
    trrust = pd.read_csv(TRRUST_PATH, sep="\t", header=None)
    trrust.columns = ["TF", "target", "interaction_type", "ids"]
    trrust = trrust[["TF", "target"]].drop_duplicates()

    chip_mat = pd.read_csv(CHIP_ENCODE_GT_PATH, sep="\t", index_col=0)
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
        gt["TF"] = gt["TF"].astype(str)
        gt["target"] = gt["target"].astype(str)
        gt = gt[gt["TF"].isin(CORE_TFS) & gt["target"].isin(valid_targets)]
        out[name] = gt.drop_duplicates(["TF", "target"]).reset_index(drop=True)
    return out


def normalize_edges(df: pd.DataFrame, method: str, valid_targets: set[str]) -> pd.DataFrame:
    df = df[["TF", "target", "EdgeWeight"]].copy()
    df["TF"] = df["TF"].astype(str).str.strip()
    df["target"] = df["target"].astype(str).str.strip()
    df["EdgeWeight"] = pd.to_numeric(df["EdgeWeight"], errors="coerce")
    df = df.dropna(subset=["TF", "target", "EdgeWeight"])
    df = df[df["TF"] != df["target"]]
    df = df[df["target"].isin(valid_targets)]
    df = df[df["TF"].isin(CORE_TFS)]
    df = (
        df.sort_values("EdgeWeight", ascending=False, kind="mergesort")
          .drop_duplicates(["TF", "target"], keep="first")
          .reset_index(drop=True)
    )
    df["method"] = method
    return df


def load_beeline_method(name: str, valid_targets: set[str]) -> pd.DataFrame:
    if name == "PPCOR":
        path = BEELINE_RESULTS_DIR / name / "outFile.txt"
        rows = []
        for chunk in pd.read_csv(path, sep="\t", chunksize=1_000_000):
            chunk = chunk.rename(columns={"Gene1": "TF", "Gene2": "target"})
            chunk = chunk[chunk["TF"].isin(CORE_TFS)]
            chunk = chunk[chunk["target"].isin(valid_targets)]
            chunk = chunk[chunk["TF"] != chunk["target"]]
            if len(chunk):
                chunk["EdgeWeight"] = chunk["corVal"].abs()
                rows.append(chunk[["TF", "target", "EdgeWeight"]])
        df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["TF", "target", "EdgeWeight"])
    else:
        path = BEELINE_RESULTS_DIR / name / "rankedEdges.csv"
        rows = []
        for chunk in pd.read_csv(path, sep="\t", chunksize=1_000_000):
            chunk = chunk.rename(columns={"Gene1": "TF", "Gene2": "target"})
            chunk = chunk[chunk["TF"].isin(CORE_TFS)]
            chunk = chunk[chunk["target"].isin(valid_targets)]
            chunk = chunk[chunk["TF"] != chunk["target"]]
            if len(chunk):
                rows.append(chunk[["TF", "target", "EdgeWeight"]])
        df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["TF", "target", "EdgeWeight"])
    return normalize_edges(df, name, valid_targets)


def break_ties(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    score = out["EdgeWeight"].astype(float)
    dense = score.rank(method="dense", ascending=False)
    n_unique = int(dense.max()) if len(dense) else 0
    base = (n_unique - (dense - 1)).to_numpy(float)
    within = out.groupby("EdgeWeight", sort=False).cumcount().to_numpy(float)
    gsize = out.groupby("EdgeWeight", sort=False)["EdgeWeight"].transform("size").to_numpy(float)
    out["EdgeWeight"] = base + 1e-6 * ((gsize - 1 - within) / gsize)
    return out


def load_chronopack(
    name: str,
    asso_path: Path,
    omnipath_path: Path,
    valid_targets: set[str],
) -> pd.DataFrame:
    asso = pd.read_csv(asso_path)
    gene_col = "X" if "X" in asso.columns else asso.columns[0]
    asso = asso.rename(columns={gene_col: "target"})
    asso["target"] = asso["target"].astype(str)
    asso = asso[asso["target"].isin(valid_targets)].copy()
    asso["gene_score"] = 1 - pd.to_numeric(asso["lfdr"], errors="coerce")
    asso = asso.dropna(subset=["gene_score"])

    omni = pd.read_csv(omnipath_path)
    omni = omni[omni["TF"].isin(CORE_TFS)].copy()
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
    pred = normalize_edges(pred, name, valid_targets)
    return break_ties(pred)


def load_replogle(valid_targets: set[str]) -> pd.DataFrame:
    with open(EN2SYM_PATH, "rb") as f:
        en2sym = pickle.load(f)
    rows = []
    for fname in sorted(glob.glob(str(REPLOGLIE_DEG_DIR / "*"))):
        tf_ensg = Path(fname).name.split(".")[0]
        tf = en2sym.get(tf_ensg)
        if tf is None:
            continue
        tmp = pd.read_csv(fname, sep="\t")
        target_col = tmp.columns[0]
        tmp["TF"] = tf
        tmp["target"] = tmp[target_col].astype(str).map(en2sym)
        tmp = tmp.dropna(subset=["target"])
        tmp["EdgeWeight"] = -np.log10(pd.to_numeric(tmp["p_val"], errors="coerce").clip(lower=np.nextafter(0, 1)))
        rows.append(tmp[["TF", "target", "EdgeWeight"]])
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["TF", "target", "EdgeWeight"])
    return normalize_edges(df, "Replogle2022_K562", valid_targets)


def load_prior_method(name: str, valid_targets: set[str]) -> pd.DataFrame:
    """Fig5_GRN_v3: a motif-prior network (../Run_PriorGRN rankedEdges.csv, BEELINE format), loaded like BEELINE."""
    df = pd.read_csv(PRIOR_RESULTS_DIR / name / "rankedEdges.csv", sep="\t")
    df = df.rename(columns={"Gene1": "TF", "Gene2": "target"})
    return normalize_edges(df[["TF", "target", "EdgeWeight"]], name, valid_targets)


def load_yao(valid_targets: set[str]) -> pd.DataFrame:
    df = pd.read_csv(YAO2023_DEGS_PATH)
    df = df.rename(columns={"Gene name": "target"})
    df["EdgeWeight"] = -np.log10(pd.to_numeric(df["P-value"], errors="coerce").clip(lower=np.nextafter(0, 1)))
    return normalize_edges(df[["TF", "target", "EdgeWeight"]], "Yao2023_THP1", valid_targets)


def load_methods(valid_targets: set[str]) -> dict[str, MethodData]:
    methods: dict[str, MethodData] = {}
    for name in ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]:
        edges = load_beeline_method(name, valid_targets)
        methods[name] = MethodData(
            name=name,
            edges=edges,
            native_tfs=frozenset(edges["TF"].unique()),
            native_targets=frozenset(edges["target"].unique()),
            context="K562 TNF scRNA",
        )

    chrono_k562 = load_chronopack("ChronoPack_K562", K562_TRADESEQ_PATH, K562_OMNIPATH_PATH, valid_targets)
    methods["ChronoPack_K562"] = MethodData(
        "ChronoPack_K562", chrono_k562,
        frozenset(chrono_k562["TF"].unique()),
        frozenset(chrono_k562["target"].unique()),
        "K562 TNF time-course",
    )

    chrono_thp1 = load_chronopack("ChronoPack_THP1", THP1_TRADESEQ_PATH, THP1_OMNIPATH_PATH, valid_targets)
    methods["ChronoPack_THP1"] = MethodData(
        "ChronoPack_THP1", chrono_thp1,
        frozenset(chrono_thp1["TF"].unique()),
        frozenset(chrono_thp1["target"].unique()),
        "THP1 LPS time-course",
    )

    replogle = load_replogle(valid_targets)
    methods["Replogle2022_K562"] = MethodData(
        "Replogle2022_K562", replogle,
        frozenset(replogle["TF"].unique()),
        frozenset(replogle["target"].unique()),
        "K562 perturb-seq",
    )

    yao = load_yao(valid_targets)
    methods["Yao2023_THP1"] = MethodData(
        "Yao2023_THP1", yao,
        frozenset(yao["TF"].unique()),
        frozenset(yao["target"].unique()),
        "THP1 perturbation",
    )

    if fig5_paths.WITH_PRIOR:  # Fig5_GRN_v3
        for name in fig5_paths.PRIOR_METHODS:
            edges = load_prior_method(name, valid_targets)
            methods[name] = MethodData(name, edges, frozenset(edges["TF"].unique()), frozenset(edges["target"].unique()),
                                       "K562 TNF scRNA + motif prior")
    return {m: methods[m] for m in METHOD_ORDER if m in methods}


def ap_random_completion(hit_idx: np.ndarray, n_pred: int, n_universe: int, k_pos: int) -> float:
    hit_ranks = hit_idx + 1
    t = len(hit_ranks)
    prefix = float(np.sum(np.arange(1, t + 1) / hit_ranks)) if t else 0.0
    m = int(n_universe - n_pred)
    p2 = int(k_pos - t)
    tail = 0.0
    if m > 0 and p2 > 0:
        if m == 1:
            tail = (p2 / m) * (t + 1) / (n_pred + 1)
        else:
            j = np.arange(1, m + 1, dtype=float)
            p = p2 / m
            num = p * (t + 1 + (j - 1) * (p2 - 1) / (m - 1))
            tail = float(np.sum(num / (n_pred + j)))
    return (prefix + tail) / float(k_pos) if k_pos else np.nan


def evaluate_one(
    edges: pd.DataFrame,
    gt: pd.DataFrame,
    tf_set: set[str],
    target_universe: set[str],
    curve_points: int = 500,
) -> tuple[dict, pd.DataFrame]:
    target_universe = set(target_universe)
    tf_set = set(tf_set)
    gt_eval = gt[gt["TF"].isin(tf_set) & gt["target"].isin(target_universe)].drop_duplicates(["TF", "target"])
    positives = set(zip(gt_eval["TF"], gt_eval["target"]))
    k_pos = len(positives)
    n_universe = len(tf_set) * len(target_universe)
    rho = k_pos / n_universe if n_universe else np.nan

    pred = edges[edges["TF"].isin(tf_set) & edges["target"].isin(target_universe)].copy()
    pred = pred.sort_values("EdgeWeight", ascending=False, kind="mergesort")
    pred = pred.drop_duplicates(["TF", "target"], keep="first").reset_index(drop=True)
    n_pred = len(pred)
    hit_mask = pd.MultiIndex.from_frame(pred[["TF", "target"]]).isin(positives)
    hit_idx = np.flatnonzero(hit_mask)
    n_hits_pred = len(hit_idx)
    positive_coverage = n_hits_pred / k_pos if k_pos else np.nan

    if k_pos == 0 or n_universe == 0:
        row = {
            "k_pos": k_pos,
            "N": n_universe,
            "rho": rho,
            "n_pred": n_pred,
            "TP_emitted": n_hits_pred,
            "positive_coverage": positive_coverage,
            "EPR": np.nan,
            "AP": np.nan,
            "AUPRC_ratio": np.nan,
            "k_EPR": np.nan,
            "TP_at_k": np.nan,
        }
        return row, pd.DataFrame(columns=["k", "TP", "precision", "recall", "is_expectation"])

    k_eval = max(1, min(k_pos, n_universe))
    extra = [1, k_eval, n_universe]
    if n_pred:
        grid_pred = np.unique(np.linspace(1, min(n_pred, n_universe), min(curve_points, min(n_pred, n_universe)), dtype=int))
        extra.extend(grid_pred.tolist())
    if n_pred < n_universe:
        grid_tail = np.unique(np.linspace(max(n_pred + 1, 1), n_universe, curve_points, dtype=int))
        extra.extend(grid_tail.tolist())
    extra.extend((hit_idx + 1).tolist())
    k_all = np.unique(np.array([k for k in extra if 1 <= k <= n_universe], dtype=int))

    # Fig5_GRN_v3 fix (2026-10-22): hit_idx holds 0-based positions, so the hits within the top k are those with
    # position < k (side="left"). The capsule / v1 / v2 code used side="right", which also counted a hit at rank k + 1
    # (EPR, TP_at_k and PR-curve TP / precision; AP and AUPRCr are computed separately and were not affected).
    tp_exact = np.searchsorted(hit_idx, np.minimum(k_all, n_pred), side="left").astype(float)
    is_exp = k_all > n_pred
    tp_all = tp_exact.copy()
    if n_pred < n_universe:
        denom = n_universe - n_pred
        rho_tail = (k_pos - n_hits_pred) / denom if denom > 0 else 0.0
        tp_all[is_exp] = n_hits_pred + (k_all[is_exp] - n_pred) * rho_tail

    precision = tp_all / k_all
    recall = tp_all / k_pos
    curve = pd.DataFrame({
        "k": k_all,
        "TP": tp_all,
        "precision": precision,
        "recall": recall,
        "is_expectation": is_exp,
    })

    i_star = np.flatnonzero(k_all == k_eval)
    i_star = int(i_star[0]) if len(i_star) else int(np.argmin(np.abs(k_all - k_eval)))
    ap = ap_random_completion(hit_idx, n_pred, n_universe, k_pos)
    row = {
        "k_pos": k_pos,
        "N": n_universe,
        "rho": rho,
        "n_pred": n_pred,
        "TP_emitted": n_hits_pred,
        "positive_coverage": positive_coverage,
        "EPR": precision[i_star] / rho if rho > 0 else np.nan,
        "AP": ap,
        "AUPRC_ratio": ap / rho if rho > 0 else np.nan,
        "k_EPR": int(k_all[i_star]),
        "TP_at_k": float(tp_all[i_star]),
    }
    return row, curve


def universe_for_strategy(
    strategy: str,
    methods: dict[str, MethodData],
    method_subset: list[str],
    broad_targets: set[str],
) -> tuple[set[str], set[str]]:
    if strategy == "broad":
        return set(CORE_TFS), set(broad_targets)
    if strategy == "common":
        tf_sets = [methods[m].native_tfs & set(CORE_TFS) for m in method_subset]
        target_sets = [methods[m].native_targets & broad_targets for m in method_subset]
        return set.intersection(*map(set, tf_sets)), set.intersection(*map(set, target_sets))
    if strategy == "union":
        tf_sets = [methods[m].native_tfs & set(CORE_TFS) for m in method_subset]
        target_sets = [methods[m].native_targets & broad_targets for m in method_subset]
        return set.union(*map(set, tf_sets)), set.union(*map(set, target_sets))
    raise ValueError(strategy)


def build_evaluations(
    methods: dict[str, MethodData],
    gts: dict[str, pd.DataFrame],
    broad_targets: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scenarios = [
        ("system_broad_all", "broad", METHOD_ORDER),
        ("common_all", "common", METHOD_ORDER),
        ("common_k562_native", "common", K562_NATIVE_METHODS),
        ("common_k562_plus_replogle", "common", K562_PLUS_REPLOGLE_METHODS),
        ("common_thp1", "common", THP1_METHODS),
    ]
    metric_rows = []
    curve_rows = []
    coverage_rows = []

    for scenario, strategy, subset in scenarios:
        subset = [m for m in subset if m in methods]
        tf_set, target_set = universe_for_strategy(strategy, methods, subset, broad_targets)
        for bench_name, gt in gts.items():
            gt_eval = gt[gt["TF"].isin(tf_set) & gt["target"].isin(target_set)]
            coverage_rows.append({
                "scenario": scenario,
                "benchmark": bench_name,
                "strategy": strategy,
                "method_subset": ",".join(subset),
                "n_tfs": len(tf_set),
                "n_targets": len(target_set),
                "gt_edges_retained": len(gt_eval),
                "gt_edges_total": len(gt),
                "gt_retention": len(gt_eval) / len(gt) if len(gt) else np.nan,
                "tfs": ",".join(sorted(tf_set)),
            })
            for method in subset:
                row, curve = evaluate_one(methods[method].edges, gt, tf_set, target_set)
                row.update({
                    "scenario": scenario,
                    "benchmark": bench_name,
                    "strategy": strategy,
                    "scope": "network",
                    "TF": "ALL",
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "method_context": methods[method].context,
                    "n_tfs_eval": len(tf_set),
                    "n_targets_eval": len(target_set),
                })
                metric_rows.append(row)
                if scenario in {"system_broad_all", "common_all", "common_k562_native", "common_thp1"}:
                    c = curve.copy()
                    c["scenario"] = scenario
                    c["benchmark"] = bench_name
                    c["scope"] = "network"
                    c["TF"] = "ALL"
                    c["method"] = method
                    curve_rows.append(c)

            for tf in CORE_TFS:
                if tf not in tf_set:
                    continue
                for method in subset:
                    row, curve = evaluate_one(methods[method].edges, gt, {tf}, target_set)
                    row.update({
                        "scenario": scenario,
                        "benchmark": bench_name,
                        "strategy": strategy,
                        "scope": "per_tf",
                        "TF": tf,
                        "method": method,
                        "method_label": METHOD_LABELS[method],
                        "method_context": methods[method].context,
                        "n_tfs_eval": 1,
                        "n_targets_eval": len(target_set),
                    })
                    metric_rows.append(row)
                    if scenario in {"system_broad_all", "common_all", "common_k562_native", "common_thp1"} and tf in DISPLAY_TFS:
                        c = curve.copy()
                        c["scenario"] = scenario
                        c["benchmark"] = bench_name
                        c["scope"] = "per_tf"
                        c["TF"] = tf
                        c["method"] = method
                        curve_rows.append(c)

    metrics = pd.DataFrame(metric_rows)
    curves = pd.concat(curve_rows, ignore_index=True) if curve_rows else pd.DataFrame()
    coverage = pd.DataFrame(coverage_rows)
    return metrics, curves, coverage


def add_summary_tables(metrics: pd.DataFrame) -> pd.DataFrame:
    per_tf = metrics[(metrics["scope"] == "per_tf") & metrics["k_pos"].gt(0)].copy()
    rows = []
    for keys, sub in per_tf.groupby(["scenario", "benchmark", "method", "method_label"], dropna=False):
        rows.append({
            "scenario": keys[0],
            "benchmark": keys[1],
            "method": keys[2],
            "method_label": keys[3],
            "n_TFs": sub["TF"].nunique(),
            "median_EPR": sub["EPR"].median(),
            "median_AUPRC_ratio": sub["AUPRC_ratio"].median(),
            "mean_EPR": sub["EPR"].mean(),
            "mean_AUPRC_ratio": sub["AUPRC_ratio"].mean(),
            "median_positive_coverage": sub["positive_coverage"].median(),
            "mean_positive_coverage": sub["positive_coverage"].mean(),
        })
    return pd.DataFrame(rows)


def plot_summary(summary: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    for scenario in ["system_broad_all", "common_all", "common_k562_native", "common_thp1"]:
        sub = summary[summary["scenario"] == scenario].copy()
        if sub.empty:
            continue
        sub["method"] = pd.Categorical(sub["method"], METHOD_ORDER, ordered=True)
        sub = sub.sort_values("method")
        fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
        for ax, metric, title in zip(axes, ["median_EPR", "median_AUPRC_ratio"], ["Median EPR", "Median AUPRC ratio"]):
            sns.barplot(
                data=sub,
                y="method_label",
                x=metric,
                hue="benchmark",
                ax=ax,
                palette={"TRRUST": "#777777", "K562_ChIP": "#4C78A8"},
            )
            ax.axvline(1.0, color="black", lw=1, ls="--", alpha=0.6)
            ax.set_title(title)
            ax.set_xlabel(metric)
            ax.set_ylabel("")
            if ax.legend_ and metric == "median_AUPRC_ratio":
                ax.legend_.remove()
        axes[0].legend(frameon=False, loc="lower right")
        fig.suptitle(f"{scenario}: median per-TF performance")
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"summary_{scenario}.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIG_DIR / f"summary_{scenario}.pdf", bbox_inches="tight")
        plt.close(fig)


def plot_coverage(coverage: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    sub = coverage[coverage["scenario"].isin(["system_broad_all", "common_all", "common_k562_native", "common_thp1"])].copy()
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.barplot(data=sub, x="scenario", y="gt_retention", hue="benchmark", ax=ax)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Ground-truth positives retained")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ground_truth_retention_by_scenario.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "ground_truth_retention_by_scenario.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_pr_curves(curves: pd.DataFrame, scenario: str, scope: str = "network") -> None:
    sns.set_theme(style="white", font_scale=0.9)
    sub = curves[(curves["scenario"] == scenario) & (curves["scope"] == scope)].copy()
    if sub.empty:
        return
    tfs = ["ALL"] if scope == "network" else DISPLAY_TFS
    fig, axes = plt.subplots(len(tfs), 2, figsize=(9, 3.3 * len(tfs)), squeeze=False)
    for i, tf in enumerate(tfs):
        for j, bench in enumerate(["TRRUST", "K562_ChIP"]):
            ax = axes[i, j]
            panel = sub[(sub["TF"] == tf) & (sub["benchmark"] == bench)]
            for method in METHOD_ORDER:
                d = panel[panel["method"] == method]
                if d.empty:
                    continue
                ax.plot(
                    d["recall"],
                    d["precision"],
                    drawstyle="steps-post",
                    lw=1.4,
                    color=METHOD_COLORS.get(method),
                    label=METHOD_LABELS.get(method, method),
                )
            ax.set_title(f"{bench} - {tf}")
            ax.set_xlabel("Recall")
            ax.set_ylabel("Precision")
            ax.set_xlim(0, 1.0)
            ymax = panel["precision"].replace([np.inf, -np.inf], np.nan).dropna().quantile(0.995) if len(panel) else 1
            ax.set_ylim(0, min(1.0, max(0.05, ymax * 1.1)))
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.suptitle(f"Precision-recall curves: {scenario} ({scope})")
    fig.tight_layout(rect=[0, 0, 0.84, 0.97])
    fig.savefig(FIG_DIR / f"pr_curves_{scenario}_{scope}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"pr_curves_{scenario}_{scope}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_per_tf_heatmap(metrics: pd.DataFrame, scenario: str) -> None:
    sns.set_theme(style="white", font_scale=0.85)
    sub = metrics[
        (metrics["scenario"] == scenario)
        & (metrics["scope"] == "per_tf")
        & (metrics["TF"].isin(CORE_TFS))
        & (metrics["k_pos"].gt(0))
    ].copy()
    if sub.empty:
        return
    for bench in ["TRRUST", "K562_ChIP"]:
        b = sub[sub["benchmark"] == bench].copy()
        if b.empty:
            continue
        b["method"] = pd.Categorical(b["method"], METHOD_ORDER, ordered=True)
        for metric in ["EPR", "AUPRC_ratio", "positive_coverage"]:
            mat = b.pivot_table(index="method", columns="TF", values=metric, aggfunc="first")
            mat = mat.reindex(index=[m for m in METHOD_ORDER if m in mat.index], columns=CORE_TFS)
            fig, ax = plt.subplots(figsize=(8, max(3, 0.35 * len(mat))))
            sns.heatmap(mat, annot=True, fmt=".2g", cmap="viridis", ax=ax, cbar_kws={"label": metric})
            ax.set_title(f"{scenario} {bench}: per-TF {metric}")
            ax.set_xlabel("")
            ax.set_ylabel("")
            fig.tight_layout()
            fig.savefig(FIG_DIR / f"heatmap_{scenario}_{bench}_{metric}.png", dpi=300, bbox_inches="tight")
            fig.savefig(FIG_DIR / f"heatmap_{scenario}_{bench}_{metric}.pdf", bbox_inches="tight")
            plt.close(fig)


def main() -> None:
    ensure_dirs()
    valid_targets = load_valid_targets()
    gts = load_ground_truths(valid_targets)
    methods = load_methods(valid_targets)

    method_summary = pd.DataFrame([
        {
            "method": m,
            "method_label": METHOD_LABELS[m],
            "context": data.context,
            "native_tfs": ",".join(sorted(data.native_tfs)),
            "n_native_tfs": len(data.native_tfs),
            "n_native_targets": len(data.native_targets),
            "n_ranked_edges": len(data.edges),
        }
        for m, data in methods.items()
    ])
    method_summary.to_csv(TABLE_DIR / "method_native_coverage.csv", index=False)

    metrics, curves, coverage = build_evaluations(methods, gts, valid_targets)
    summary = add_summary_tables(metrics)

    metrics.to_csv(TABLE_DIR / "metrics_long.csv", index=False)
    curves.to_csv(TABLE_DIR / "pr_curves_long.csv", index=False)
    coverage.to_csv(TABLE_DIR / "benchmark_universe_coverage.csv", index=False)
    summary.to_csv(TABLE_DIR / "summary_by_method.csv", index=False)

    plot_summary(summary)
    plot_coverage(coverage)
    for scenario in ["system_broad_all", "common_all", "common_k562_native", "common_thp1"]:
        plot_pr_curves(curves, scenario, scope="network")
        plot_pr_curves(curves, scenario, scope="per_tf")
        plot_per_tf_heatmap(metrics, scenario)

    print("Wrote:")
    for p in sorted(TABLE_DIR.glob("*")):
        print(p)
    for p in sorted(FIG_DIR.glob("*.png")):
        print(p)


if __name__ == "__main__":
    main()
