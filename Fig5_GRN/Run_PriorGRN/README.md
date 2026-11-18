# Run_PriorGRN: motif-prior GRN methods on the r9 single cells

Started 2026-10-21. Adds three popular GRN methods that combine scRNA-seq with a TF-binding-motif prior to the
Run_BEELINE comparison (GENIE3, GRNBOOST2, PPCOR, SINCERITIES), on exactly the same cells and genes.

## Why these methods (and not SCENIC+ or Pando)

There is no ATAC-seq (single-cell or bulk) for the K562 + TNF-alpha or THP1 + LPS experiments. SCENIC+ and Pando
need per-cell accessibility as model input, not only as a prior: SCENIC+ derives candidate regions from scATAC topics
and links regions to genes by regressing expression on accessibility across metacells; Pando models each target as a
sum of TF x region-accessibility terms, keeping regions whose accessibility correlates with the target. A promoter
base GRN cannot fill that role (with constant accessibility both reduce to methods already included: SCENIC+ to
pySCENIC, Pando to a CellOracle-like regression). The user chose (2026-10-21) to run the three methods below, which use
a motif prior as designed:

| Method | Version | Prior | Network model | EdgeWeight |
|---|---|---|---|---|
| PYSCENIC | pySCENIC 0.12.1 (docker aertslab/pyscenic:0.12.1) | cisTarget hg38 refseq_r80 v10 gene-based motif rankings, 500 bp up / 100 bp down and 10 kb up / 10 kb down of the TSS; motif2TF hgnc v10nr_clust | GRNBoost2 (pySCENIC's TF list as regulators) -> ctx module pruning (`--mask_dropouts`, activating modules) | GRNBoost2 importance of the regulon edges |
| CELLORACLE | CellOracle 0.18.0 (docker kenjikamimoto126/celloracle_ubuntu:0.18.0) | CellOracle's built-in hg38 promoter base GRN (`hg38_gimmemotifsv5_fpr2`; 39,315 promoter windows of 1.1 kb, 1,094 TFs, 27,314 genes) | KNN imputation + bagging ridge per target (alpha 10, 20 bags), one GRN unit = all cells | `coef_abs` (CellOracle's ranking) |
| INFERELATOR | Inferelator 0.6.4 (pip, `.conda_env`) | the same CellOracle base GRN as a 0/1 gene x TF prior | TF activity (TFA) + BBSR, 10 bootstraps | `combined_confidences` |

Notes:
- CellOracle's `load_human_promoter_base_GRN()` defaults to **hg19**; the hg38 file is used here (copied from the
  local CellOracle repo, `input/`).
- The base GRN is permissive (e.g. EGR1 has a motif in 68% of promoter windows, FOS 46%); all 12 Fig5 TFs are in it.
- No method sees a Fig5 ground truth: Inferelator runs with `use_no_gold_standard` (its gold standard is only used
  for its own precision/recall report).
- pySCENIC predicts only the cisTarget-supported regulon edges, so its networks are much smaller than the others;
  CellOracle ranks every base-GRN edge it fits; Inferelator ranks every edge with confidence > 0.

## Data

Datasets `{K562,THP1}-r9-sc-cpfilt-p05-TFs{1000,2000}`, i.e. the Run_BEELINE runs with no expression filter and the
BEELINE paper's selection at an uncorrected GAM P < 0.05 (significant TFs + the N most-varying significant genes),
Fig5 TFs forced. Genes and cells are read from Run_BEELINE's `ExpressionData.csv` and checked against the r9 AnnData.
Expression is r9's log1p layer (= log1p of its cpm layer) for all three methods, the matrix the BEELINE runs use
(CellOracle is given the cpm layer as "raw" counts and applies log1p itself; the script checks the result equals
the log1p layer).

## Steps

| Step | Script | What |
|---|---|---|
| setup | `scripts/slurm/setup.sbatch` | pulls both containers into `resources/containers/`, downloads the cisTarget files (sha1-checked) into `resources/cistarget/`, creates the project-local `.conda_env` (inferelator) |
| 01 | `scripts/01_prepare_inputs.py` (`.conda_env`) | per dataset: `output/inputs/<dataset>/{expr.h5ad, expr_log1p.tsv, tfs_pyscenic.txt, prior.tsv, tfs_inferelator.txt, summary.json}` |
| PYSCENIC | `scripts/slurm/method.sbatch` -> `arboreto_with_multiprocessing.py`, `pyscenic ctx`, `scripts/02_pyscenic_edges.py` | GRNBoost2 without dask (avoids the shutdown hangs seen in Run_BEELINE), seed 777 |
| CELLORACLE | `scripts/03_celloracle.py` (container) | tutorial settings (PCA elbow <= 50 dims, k = 2.5% of cells) |
| INFERELATOR | `scripts/04_inferelator.py` (`.conda_env`) | |

Outputs: `output/networks/<dataset>/<METHOD>/rankedEdges.csv` (BEELINE format: tab-separated `Gene1 Gene2 EdgeWeight`,
sorted by weight), `summary.json`, `time_seconds.txt`, plus the method's full output (`adjacencies.tsv` +
`regulons.csv`, `links_all.csv`, `network_all.tsv`; large ones git-ignored).

```bash
cd Run_PriorGRN
./run_prior_grn.sh setup        # once
./run_prior_grn.sh all          # step 01, then 3 methods x 4 datasets (12 jobs)
```

Scoring with the Fig5 metrics is not wired in yet (Run_BEELINE steps 05-07 are pinned to K562 and the r9_TFs dataset
names).

## Results

All 12 runs completed on 2026-10-22 (jobs 943020-943031). Not yet scored with the Fig5 metrics.

| Dataset | Method | Edges | TFs | Targets | Runtime | Fig5 TFs with no edges |
|---|---|---|---|---|---|---|
| K562-r9-sc-cpfilt-p05-TFs1000 | PYSCENIC | 5,630 | 86 | 1,197 | 6 min | RELA, ATF2, CREB1 |
| K562-r9-sc-cpfilt-p05-TFs1000 | CELLORACLE | 48,530 | 137 | 1,277 | 2 min | - |
| K562-r9-sc-cpfilt-p05-TFs1000 | INFERELATOR | 16,126 | 137 | 1,310 | 60 min | - |
| K562-r9-sc-cpfilt-p05-TFs2000 | PYSCENIC | 7,555 | 101 | 1,930 | 10 min | RELA, ATF2, ATF1 |
| K562-r9-sc-cpfilt-p05-TFs2000 | CELLORACLE | 86,198 | 148 | 2,214 | 4 min | - |
| K562-r9-sc-cpfilt-p05-TFs2000 | INFERELATOR | 23,803 | 148 | 2,258 | 96 min | - |
| THP1-r9-sc-cpfilt-p05-TFs1000 | PYSCENIC | 6,807 | 78 | 1,040 | 2 min | RELA, ELK1 |
| THP1-r9-sc-cpfilt-p05-TFs1000 | CELLORACLE | 37,951 | 108 | 1,149 | 7 min | ELK1 |
| THP1-r9-sc-cpfilt-p05-TFs1000 | INFERELATOR | 8,886 | 107 | 1,180 | 38 min | ELK1 |
| THP1-r9-sc-cpfilt-p05-TFs2000 | PYSCENIC | 7,182 | 76 | 1,280 | 2 min | ATF2, ELK1 |
| THP1-r9-sc-cpfilt-p05-TFs2000 | CELLORACLE | 47,992 | 110 | 1,460 | 8 min | ELK1 |
| THP1-r9-sc-cpfilt-p05-TFs2000 | INFERELATOR | 9,676 | 110 | 1,452 | 40 min | ELK1 |

ELK1 has no edges in THP1 because it is not among the r9 THP1 genes. pySCENIC drops TFs whose module has no enriched motif (e.g. no RELA regulon on THP1 TFs + 1000), so its networks are small and miss some Fig5 TFs.
