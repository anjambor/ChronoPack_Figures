# Copied from grn_inference_v2/reproduction/scripts; path roots redirected to Fig5_GRN via fig5_paths.
from __future__ import annotations

import glob
import pickle
from pathlib import Path

import fig5_paths

import numpy as np
import pandas as pd


BASE_DIR = fig5_paths.BENCH_ROOT
TABLE_DIR = BASE_DIR / "benchmark_outputs" / "tables"

CHRONO_ANALYSES_DIR = (fig5_paths.EXTERNAL_DIR / "projects/chrono_seq/analyses")
K562_ANALYSIS_DIR = CHRONO_ANALYSES_DIR / "k562-tnfa.hg38.rep"
GRN_INFERENCE_DIR = (fig5_paths.EXTERNAL_DIR / "projects/grn_inference")

PROTEIN_CODING_GENES_PATH = (fig5_paths.EXTERNAL_DIR / "ref/genome/hg38/annotations/protein_coding_genes.csv")
TRRUST_PATH = (fig5_paths.EXTERNAL_DIR / "ref/grns/ttrust/trrust_rawdata.human.tsv")
CHIP_ENCODE_GT_PATH = (fig5_paths.EXTERNAL_DIR / "downloads/k562-chip_encode/results/tf_gene.1kbp_200bp.stringent.tsv")
PRIOR_RESULTS_DIR = fig5_paths.PRIOR_RESULTS_DIR  # Fig5_GRN_v3: motif-prior networks (FIG5_METHOD_SET=with_prior)
BEELINE_RESULTS_DIR = fig5_paths.BEELINE_RESULTS_DIR  # Fig5_GRN_v3: r9 P < 0.05 TFs + 2000 run (was K562-3635genes-231metacells)
K562_TRADESEQ_PATH = K562_ANALYSIS_DIR / "intmd/sc-K562-1e4_to_tradeSeq/asso_df.1trajectory.lfdr.csv"
K562_OMNIPATH_PATH = K562_ANALYSIS_DIR / "intmd/omnipath_tf_ranks/omnipath_ranks.TNFa.csv"
THP1_TRADESEQ_PATH = CHRONO_ANALYSES_DIR / "thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq/asso_df.1trajectory.lfdr.csv"
THP1_OMNIPATH_PATH = K562_ANALYSIS_DIR / "intmd/omnipath_tf_ranks/omnipath_ranks.TLR4.csv"
EN2SYM_PATH = GRN_INFERENCE_DIR / "causalbench/data_cache/en2sym.pkl"
REPLOGLE_DEG_DIR = (fig5_paths.EXTERNAL_DIR / "downloads/perturb_db/seurat_perturbdb_workflow/results/deg_tables")
YAO2023_DEGS_PATH = GRN_INFERENCE_DIR / "causalbench/notebooks/Yao2023/TFs_vs_non_targeting_DE_all.csv"

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

CUTOFFS = [1000, 5000, 10000]


def load_valid_targets() -> set[str]:
    genes = pd.read_csv(PROTEIN_CODING_GENES_PATH, header=None)[0].dropna().astype(str)
    # FIG5_UNIVERSE=shared: cut to the genes every non-Perturb-seq method ranks (fig5_paths.restrict_universe)
    return fig5_paths.restrict_universe({g for g in genes if not (g.startswith("MT-") or g.startswith("ENSG"))})


def load_reference_tfs(valid_targets: set[str]) -> set[str]:
    trrust = pd.read_csv(TRRUST_PATH, sep="\t", header=None)
    trrust_tfs = set(trrust[trrust[1].astype(str).isin(valid_targets)][0].dropna().astype(str).str.strip())

    chip_mat = pd.read_csv(CHIP_ENCODE_GT_PATH, sep="\t", index_col=0)
    chip_tfs = set(chip_mat.columns.astype(str).str.strip())
    return trrust_tfs | chip_tfs


def normalize_edges(df: pd.DataFrame, method: str, valid_targets: set[str]) -> pd.DataFrame:
    df = df[["TF", "target", "EdgeWeight"]].copy()
    df["TF"] = df["TF"].astype(str).str.strip()
    df["target"] = df["target"].astype(str).str.strip()
    df["EdgeWeight"] = pd.to_numeric(df["EdgeWeight"], errors="coerce")
    df = df.dropna(subset=["TF", "target", "EdgeWeight"])
    df = df[df["TF"] != df["target"]]
    df = df[df["target"].isin(valid_targets)]
    df = (
        df.sort_values("EdgeWeight", ascending=False, kind="mergesort")
        .drop_duplicates(["TF", "target"], keep="first")
        .reset_index(drop=True)
    )
    df["method"] = method
    return df


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


def load_beeline_method(name: str, valid_targets: set[str]) -> pd.DataFrame:
    if name == "PPCOR":
        path = BEELINE_RESULTS_DIR / name / "outFile.txt"
        score_col = "corVal"
    else:
        path = BEELINE_RESULTS_DIR / name / "rankedEdges.csv"
        score_col = "EdgeWeight"

    rows = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=1_000_000):
        chunk = chunk.rename(columns={"Gene1": "TF", "Gene2": "target"})
        chunk["TF"] = chunk["TF"].astype(str).str.strip()
        chunk["target"] = chunk["target"].astype(str).str.strip()
        chunk = chunk[chunk["target"].isin(valid_targets)]
        chunk = chunk[chunk["TF"] != chunk["target"]]
        if len(chunk):
            chunk["EdgeWeight"] = pd.to_numeric(chunk[score_col], errors="coerce")
            if name == "PPCOR":
                chunk["EdgeWeight"] = chunk["EdgeWeight"].abs()
            rows.append(chunk[["TF", "target", "EdgeWeight"]])
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["TF", "target", "EdgeWeight"])
    return normalize_edges(df, name, valid_targets)


def load_chronopack_method(
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
    omni["TF"] = omni["TF"].astype(str).str.strip()
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
    return break_ties(normalize_edges(pred, name, valid_targets))


def load_replogle(valid_targets: set[str]) -> pd.DataFrame:
    with open(EN2SYM_PATH, "rb") as f:
        en2sym = pickle.load(f)
    rows = []
    for fname in sorted(glob.glob(str(REPLOGLE_DEG_DIR / "*"))):
        tf_ensg = Path(fname).name.split(".")[0]
        tf = en2sym.get(tf_ensg)
        if tf is None:
            continue
        tmp = pd.read_csv(fname, sep="\t")
        target_col = tmp.columns[0]
        tmp["TF"] = tf
        tmp["target"] = tmp[target_col].astype(str).map(en2sym)
        tmp = tmp.dropna(subset=["target"])
        p_val = pd.to_numeric(tmp["p_val"], errors="coerce").clip(lower=np.nextafter(0, 1))
        tmp["EdgeWeight"] = -np.log10(p_val)
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
    p_val = pd.to_numeric(df["P-value"], errors="coerce").clip(lower=np.nextafter(0, 1))
    df["EdgeWeight"] = -np.log10(p_val)
    return normalize_edges(df[["TF", "target", "EdgeWeight"]], "Yao2023_THP1", valid_targets)


def load_methods(valid_targets: set[str]) -> dict[str, pd.DataFrame]:
    methods = {name: load_beeline_method(name, valid_targets) for name in ["GENIE3", "SINCERITIES", "GRNBOOST2", "PPCOR"]}
    methods["ChronoPack_K562"] = load_chronopack_method(
        "ChronoPack_K562", K562_TRADESEQ_PATH, K562_OMNIPATH_PATH, valid_targets
    )
    methods["ChronoPack_THP1"] = load_chronopack_method(
        "ChronoPack_THP1", THP1_TRADESEQ_PATH, THP1_OMNIPATH_PATH, valid_targets
    )
    methods["Replogle2022_K562"] = load_replogle(valid_targets)
    methods["Yao2023_THP1"] = load_yao(valid_targets)
    if fig5_paths.WITH_PRIOR:
        methods.update({name: load_prior_method(name, valid_targets) for name in fig5_paths.PRIOR_METHODS})
    return {name: methods[name] for name in METHOD_ORDER if name in methods}


def count_top_tfs(
    methods: dict[str, pd.DataFrame],
    tf_filter: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    long_rows = []
    wide_rows = []
    max_k = max(CUTOFFS)

    for method in METHOD_ORDER:
        if method not in methods:
            continue
        pred = methods[method]
        if tf_filter is not None:
            pred = pred[pred["TF"].isin(tf_filter)]
        pred = pred.sort_values("EdgeWeight", ascending=False, kind="mergesort").reset_index(drop=True)
        top = pred.head(max_k)
        tf_union = sorted(top["TF"].unique())
        counts_by_k = {}
        for cutoff in CUTOFFS:
            effective_k = min(cutoff, len(pred))
            counts = pred.head(effective_k)["TF"].value_counts()
            counts_by_k[cutoff] = counts
            for tf, count in counts.sort_index().items():
                long_rows.append(
                    {
                        "method": method,
                        "cutoff": cutoff,
                        "effective_k": effective_k,
                        "TF": tf,
                        "count": int(count),
                        "fraction": float(count / effective_k) if effective_k else np.nan,
                    }
                )

        for tf in tf_union:
            row = {
                "method": method,
                "TF": tf,
                "n_ranked_edges": len(pred),
                "n_tfs_in_top10000": top["TF"].nunique(),
            }
            for cutoff in CUTOFFS:
                effective_k = min(cutoff, len(pred))
                count = int(counts_by_k[cutoff].get(tf, 0))
                row[f"top_{cutoff}_effective_k"] = effective_k
                row[f"top_{cutoff}_count"] = count
                row[f"top_{cutoff}_fraction"] = float(count / effective_k) if effective_k else np.nan
            wide_rows.append(row)

    return pd.DataFrame(long_rows), pd.DataFrame(wide_rows)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    valid_targets = load_valid_targets()
    reference_tfs = load_reference_tfs(valid_targets)
    methods = load_methods(valid_targets)

    long, wide = count_top_tfs(methods)
    ref_long, ref_wide = count_top_tfs(methods, tf_filter=reference_tfs)

    long_out = TABLE_DIR / "top_predicted_tf_counts_by_method_long.csv"
    wide_out = TABLE_DIR / "top_predicted_tf_counts_by_method_wide.csv"
    ref_long_out = TABLE_DIR / "top_reference_tf_counts_by_method_long.csv"
    ref_wide_out = TABLE_DIR / "top_reference_tf_counts_by_method_wide.csv"
    long.to_csv(long_out, index=False)
    wide.to_csv(wide_out, index=False)
    ref_long.to_csv(ref_long_out, index=False)
    ref_wide.to_csv(ref_wide_out, index=False)
    print(f"Wrote {long_out}")
    print(f"Wrote {wide_out}")
    print(f"Wrote {ref_long_out}")
    print(f"Wrote {ref_wide_out}")


if __name__ == "__main__":
    main()
