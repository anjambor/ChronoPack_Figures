# Run_BEELINE: BEELINE comparators on the r9 single-cell K562 + TNF-alpha and THP1 + LPS data

Reruns the four BEELINE methods that Fig5_GRN compares against (GENIE3, GRNBOOST2, PPCOR, SINCERITIES) on the
same r9 cells ChronoPack fits, following the BEELINE paper's scRNA-seq procedure (Pratapa et al. 2020,
*Nat Methods*, Methods). No ChronoPack or tradeSeq result is used to build the inputs.

There are two contexts, each with two gene sets, **TFs + 1000** and **TFs + 2000**:
- **K562:** the r9 K562 + TNF-alpha cells. Datasets `K562-r9-sc-TFs{1000,2000}`.
- **THP1:** the THP1 cells of the r9 THP1 + K562 + LPS experiment, added 2026-10-20. Datasets
  `THP1-r9-sc-TFs{1000,2000}`.

How the predictions are used:
- `../Fig5_GRN` is unchanged; its comparators come from the original K562 metacell run (see "Original run" below).
- `../Fig5_GRN_v2` uses `K562-r9-sc-TFs2000`.
- The THP1 predictions are not used by any figure yet. Fig5's only THP1 comparator is Yao 2023.

## Design (K562)

| | Original Fig5 run (Jan 2026) | Run_BEELINE (r9) |
|---|---|---|
| Cells | 231 k-means metacells (20 per time point, pre-r9 processing) | the 6,857 r9 prefiltered K562 + TNF-alpha single cells, time tags 1-12 |
| Expression | metacell counts, normalize 1e4 + log1p | r9 `log1p` layer (= log1p CPM; the dataset's own normalization, as the paper prescribes) |
| Gene filter | 5,000 seurat HVGs (batch = trajectory cluster), time-tag-bias / MT / RP / ENSG removed | BEELINE: expressed (count > 0) in >= 10% of cells: 3,611 of 11,866 r9 genes |
| Gene set | 3,635 + 8 TFs = 3,643 | BEELINE `generateExpInputs.py -p 0 -t -n N`: all 366 expressed BEELINE TFs + the N most-varying non-TF genes, + 5 forced Fig5 TFs: **1,371** (N = 1000) and **2,371** (N = 2000) |
| Pseudotime | time tag | time tag (only SINCERITIES uses it) |
| Algorithms | BEELINE defaults; PPCOR pearson, pVal 0.01; SINCERITIES nBins 10 | same (SINCERITIES' quantile binning of the 12 time tags into 10 bins merges tags 1+2 and 6+7) |

**Deviations from the paper's gene selection** (agreed 2026-10-19):

1. **No time-test filter.** The paper keeps genes whose GAM-along-pseudotime P is < 0.01 after Bonferroni
   correction, then takes significant TFs + the N most-varying significant genes. Step 02 computes this test
   (`gam(z ~ lo(t))` from the R `gam` package, as in the Slingshot workflow; P = parametric-ANOVA `Pr(>F)` of
   `lo(t)`), but on these sparse single cells only 277 of 3,611 genes pass (58 TFs), too few for TFs + 1000.
   Looser cutoffs do not reach it either (uncorrected P < 0.05: 1,472). `generateExpInputs.py` is therefore run
   with its p-value filter off (`-p 0`), which keeps every TF and ranks the other genes by variance only. The P
   values stay in `output/prep/<context>/GeneOrdering.csv` (columns `VGAMpValue`, `NparPValue`, `OverallPValue`) for
   reference. THP1 is the same: 160 of 3,936 genes pass.
2. **Forced Fig5 TFs.** RELA, NFKB1, JUN, ATF2, ELK1, JUNB, EGR1, CREB1, ATF1, FOS, REL and IRF3 are always in
   the input. Seven pass the 10% filter and are BEELINE TFs already; RELA (6.9% of cells), JUN (9.1%), CREB1
   (6.1%), ELK1 (3.6%) and FOS (1.5%) do not and are appended (`output/prep/K562/forced_tfs.csv`).

## THP1 context (added 2026-10-20)

Same procedure, code and settings as K562; only the cells differ (`BEELINE_CONTEXT=THP1`).

| | THP1 |
|---|---|
| Cells | the 787 r9 prefiltered THP1 cells of the THP1 + K562 + LPS experiment (`adata_thp1_prefiltered.fixed.h5ad`); the same 787 cells as ChronoPack's r9 THP1 tradeSeq fit. Time tags 1-12: 97, 96, 74, 65, 68, 60, 80, 66, 57, 46, 38 and 40 cells |
| Expression | r9 `log1p` layer (= log1p CPM) |
| Gene filter | expressed in >= 10% of cells: 3,936 of 10,117 r9 genes |
| Gene set | all 477 expressed BEELINE TFs + the N most-varying non-TF genes + 4 forced Fig5 TFs: **1,481** (N = 1000) and **2,481** (N = 2000) |
| Forced TFs | RELA (4.4% of cells), EGR1 (1.7%), ATF1 (5.8%) and FOS (6.7%) are below the filter and appended. The other Fig5 TFs pass it, except ELK1, which is not among the r9 THP1 genes. ELK1 is therefore absent (`output/prep/THP1/forced_tfs.csv`), so the THP1 runs have no ELK1 predictions. |
| refNetwork | the K562 ChIP network, as for K562 (BEELINE's own evaluation only; Fig5 also scores THP1 methods against it) |
| PPCOR | **not applicable for THP1** (user decision, 2026-10-22). Every THP1 dataset has more genes (1,211-3,004) than cells (787), so ppcor inverts the singular Spearman correlation matrix with a pseudo-inverse and its p-values have n - p < 0 degrees of freedom (all NA). BEELINE's parser keeps edges with p <= 0.01 (ranked by \|r\|) and p > 0.01 (weight 0), so NA edges are dropped: every THP1 PPCOR `rankedEdges.csv` holds only the self-edges (one per gene). Not re-ranked by \|r\| without the filter: the pseudo-inverse values are not partial correlations, and K562's PPCOR ranks only the p <= 0.01 tier (1.2% of edges) by \|r\|. |

## ChronoPack expression filter (added 2026-10-21)

The same procedure with ChronoPack's expression filter instead of BEELINE's 10% filter
(`BEELINE_GENE_FILTER=chronopack`). Prep files go to `output/prep/<context>-cpfilt/`; the datasets are
`<context>-r9-sc-cpfilt-TFs{1000,2000}`.

- **The filter:** none on top of the r9 prefiltered genes. Those genes are exactly what ChronoPack fits with tradeSeq
  and scores in its GRN inference. Its filters are applied upstream:
  - detected in at least 3 cells (`01_PreprocessDatasets`);
  - at least 10 total counts;
  - MT, RP and ENSG genes and time-tag-bias genes removed (`03_PrepareForTradeSeq`).
- **Everything else is unchanged:** `-p 0`, all BEELINE TFs plus the N most-varying non-TF genes, the forced Fig5
  TFs, and the algorithm settings.

| | K562 | THP1 |
|---|---|---|
| Genes eligible | all 11,866 r9 genes (10% filter: 3,611) | all 10,117 (10% filter: 3,936) |
| BEELINE TFs | 1,032 (10% filter: 366 + 5 forced) | 1,004 (10% filter: 477 + 4 forced) |
| Gene sets | **2,032** (N = 1000), **3,032** (N = 2000) | **2,004**, **3,004** |
| Forced Fig5 TFs appended | none (all are BEELINE TFs in the r9 genes) | none (ELK1 is still absent from the r9 THP1 genes) |

**What changes compared with the 10% filter:** only the TFs.
- **Non-TF genes are the same:** the N most-varying non-TF genes are identical in both modes, for both contexts and
  both N. The 10% filter never removed a gene that would have made the variance cut-off; the least-expressed selected
  non-TF gene is in 13–25% of cells.
- **TFs are added:** the extra genes are 661 (K562) and 523 (THP1) TFs expressed in under 10% of cells. These enter
  as both regulators and targets.
- **So the sparse targets stay out:** the non-TF targets that BEELINE does not cover remain uncovered. This variant
  tests the TF part of the filter only.

GENIE3 now requests 90G and excludes `bigcompute`, where `THP1-r9-sc-TFs2000` GENIE3 hung (2026-10-20).

## Paper selection with an uncorrected P < 0.05 cut-off (added 2026-10-21)

Datasets `<context>-r9-sc-cpfilt-p05-TFs{1000,2000}`, built with `BEELINE_GENE_FILTER=chronopack BEELINE_PVAL=0.05`. They
use the no-filter prep files; only step 03 is rerun.

The selection is BEELINE's own, through `generateExpInputs.py -p 0.05 -t -n N` with no `-c`:
1. Drop genes with GAM P ≥ 0.05.
2. Include every significant TF.
3. Add the N significant non-TF genes with the highest variance.
4. Append the forced Fig5 TFs.

The only difference from the paper is the threshold: uncorrected P < 0.05 instead of Bonferroni P < 0.01. At the
paper's threshold, 303 K562 and 172 THP1 genes pass, too few for TFs + N.

| | Significant TFs | Non-TF genes | Forced TFs added | Genes |
|---|---|---|---|---|
| K562 TFs + 1000 | 324 | 1,000 | 7 | 1,331 |
| K562 TFs + 2000 | 324 | 2,000 | 7 | 2,331 |
| THP1 TFs + 1000 | 204 | 1,000 | 7 | 1,211 |
| THP1 TFs + 2000 | 204 | 1,358 (all that pass) | 7 | **1,569** |

THP1 TFs + 2000 is short of 2,000 non-TF genes because only 1,358 pass the cut-off. The user accepted this on
2026-10-21.

**Why these runs were added:** the variance-only selection (`-p 0`) leaves out many of ChronoPack's top targets,
with or without the 10% filter. These runs cover more of them. Coverage of ChronoPack's top-k targets
(`Fig5_GRN/output/chronopack/<ctx>/grn.target_scores.csv` order):

| Top k | K562 TFs+1000 | K562 TFs+2000 | THP1 TFs+1000 | THP1 TFs+2000 |
|---|---|---|---|---|
| 100 | 61% → 83% | 76% → 99% | 79% → 99% | 92% → 99% |
| 500 | 52% → 71% | 69% → 90% | 54% → 86% | 73% → 89% |
| 1000 | 43% → 60% | 61% → 82% | 41% → 70% | 59% → 76% |

(10% filter runs → P < 0.05 runs.) An emulation of `generateExpInputs.py` was checked to reproduce the existing
`-p 0` gene sets exactly. The check is step 09, `scripts/09_chronopack_target_coverage.py`, and its table is
`output/chronopack_target_coverage.csv`.

The `-cpfilt` (no filter, `-p 0`) runs are kept for comparison. They differ from the 10% runs only by about 600
added low-expression TFs.

## Run

From this directory (submits Slurm jobs):

```bash
./run_beeline.sh                     # per context (K562, THP1): prep (steps 01-03), then its 8 algorithm jobs
CONTEXTS=THP1 ./run_beeline.sh       # one context
SKIP_PREP=1 ./run_beeline.sh         # algorithm jobs only
ALGOS="GRNBOOST2" DATASETS="K562-r9-sc-TFs1000" SKIP_PREP=1 CONTEXTS=K562 ./run_beeline.sh
BEELINE_GENE_FILTER=chronopack ./run_beeline.sh   # ChronoPack's expression filter (submitted 2026-10-21)
BEELINE_GENE_FILTER=chronopack BEELINE_PVAL=0.05 PREP_STEPS=03 ./run_beeline.sh   # + paper selection, P < 0.05
```

`BEELINE_CONTEXT` (K562 by default) selects the context for steps 01-03. Their files go to `output/prep/<context>/`.
Steps 05 and 06 are K562 only.

Steps 05, 06 and 10 score with Fig5's code through `scripts/fig5_code.py`: since 2026-11-16 `..` (Fig5_GRN_v4)'s copy
(before: Fig5_GRN_v3's, the same code), which has the EPR top-k fix (hits at ranks 1..k; Fig5's original code also
counted rank k + 1). **2026-11-16 (user request: no dependencies on `../old/`):** the metacell round of Fig5 v1 and
the checks against v1's tables were dropped from steps 05, 06, 10 and 11, step 07 (TRRUST targets lost vs the metacell
input) was removed, step 09 reads ChronoPack's target scores from `..` (Fig5_GRN_v4) (identical to v1's) and step 11
takes Fig5's 5b plotting code from `..` (Fig5_GRN_v4). Rerun 2026-11-16: every BEELINE row is unchanged; only
ChronoPack_THP1 rows moved (v4 rankings: CREB1 no longer scored), e.g. THP1 TRRUST P < 0.05 TFs + 2000 median EPR
4.77 -> 3.13, AUPRCr 4.07 -> 4.26.

| Step | Script | What |
|---|---|---|
| 01 | `scripts/01_prepare_expression.py` (chronopack env) | r9 cells, `log1p` layer, 10% expression filter -> `output/prep/<context>/ExpressionData.full.csv`, `PseudoTime.csv`, `gene_filter.csv`, `forced_tfs.csv`, `prep_summary.json` |
| 02 | `scripts/02_gene_ordering.R` (`.conda_env`, R 4.3 + gam 1.22-6) | BEELINE gene ordering: per-gene variance and GAM P values -> `output/prep/<context>/GeneOrdering.csv` |
| 03 | `scripts/03_generate_inputs.py` | BEELINE's `generateExpInputs.py` (unchanged) for N = 1000 / 2000, forced TFs appended, refNetwork from the K562 ChIP ground truth, one config per algorithm -> `beeline/BEELINE-data/{inputs,config-files}/scRNA-Seq/<context>-r9-sc-TFs{1000,2000}/`, `output/inputs_summary.csv` (all contexts), `output/gene_selection_*.csv` |
| 04 | `scripts/slurm/beeline_algo.sbatch` (BEELINE env) | `BLRunner.py` for one (dataset, algorithm) -> `beeline/outputs/scRNA-Seq/<dataset>/<ALGO>/` (`rankedEdges.csv`; PPCOR also `outFile.txt`; `time.txt` = runtime and memory) |
| 05 | `scripts/05_fig5_per_tf_metrics.py` (`sbatch scripts/slurm/fig5_metrics.sbatch`) | Fig5 5b per-TF EPR / AUPRCr of the r9 predictions, scored with Fig5's code -> `output/fig5_metrics/` |
| 06 | `scripts/06_gene_universe_sensitivity.py` (`sbatch scripts/slurm/gene_universe.sbatch`) | Fig5's 04e matched-gene sensitivity with each run's BEELINE predictions -> `output/gene_universe/` |
| 09 | `scripts/09_chronopack_target_coverage.py` (chronopack env) | share of ChronoPack's top-k targets in each TFs + N gene set (current runs and the P-cut-off selections; also checks its `generateExpInputs.py` emulation) -> `output/chronopack_target_coverage.csv` |
| 10 | `scripts/10_score_all_rounds.py` (`sbatch scripts/slurm/rounds.sbatch`, ~40 min) | Fig5 5b metrics + 04e matched-gene analysis for every round (10% filter, no filter, P < 0.05) x TFs + 1000 / 2000 x K562 / THP1; THP1 BEELINE scored with the THP1 TF set; THP1 PPCOR not scored. Checks: K562 10% rows reproduce step 06; ChronoPack / Replogle / Yao rows equal Fig5's table -> `output/rounds/{per_tf_metrics,median_by_universe,runs}.csv` |
| 11 | `scripts/11_plot_p05_figures.py` (`sbatch scripts/slurm/p05_figures.sbatch`) | 5b-style strip-forests for the P < 0.05 round (full and matched universes) and round-comparison medians -> `output/rounds/figures/` |

Steps 01-03 run in `scripts/slurm/prep.sbatch` (K562 ~3 min, THP1 ~1 min). Resources per algorithm are set in `run_beeline.sh`
(GENIE3 40 CPUs, GRNBOOST2 24, PPCOR / SINCERITIES 4).

Results: `REPORT.md` (its metacell columns are from the runs before 2026-11-16, kept as a record).

## BEELINE code (`beeline/`)

`beeline/Beeline/` is a copy (2026-10-19) of `~/projects/grn_inference/BEELINE/Beeline` (BLRunner.py, BLRun/,
BLEval/, generateExpInputs.py, BLEvaluator.py), the local fork that runs the BEELINE containers with Singularity
instead of Docker. The runners treat `beeline/` as the BEELINE root (bound to `/data` in the containers), so all
inputs and outputs stay in this folder. The Singularity images are used in place from
`~/projects/grn_inference/BEELINE/singularity_images/` (`arboreto.sif` for GENIE3 / GRNBOOST2, `ppcor.sif`,
`sincerities.sif`). BEELINE's own evaluation (`BLEvaluator.py`) is not run; Fig5 evaluates the predictions.

## Inputs (`input/`)

- `adata_prefiltered.1e4-variant.fixed.h5ad`: symlink to the r9 K562 + TNF-alpha AnnData (6,857 cells x 11,866 genes).
- `adata_thp1_prefiltered.fixed.h5ad`: symlink to the r9 THP1 AnnData of the THP1 + K562 + LPS experiment
  (`thp1-k562-lps.hg38.rep/output/`, 787 cells x 10,117 genes).
- `human-tfs.csv`: BEELINE-Networks human TF list (1,563 TFs; BEELINE paper), copied.
- `tf_gene.1kbp_200bp.stringent.tsv`: K562 ENCODE ChIP ground truth, copied from `../Fig5_GRN/input/capsule/`; used
  only to write BEELINE's `refNetwork.csv`.

## Original run (for reference)

Fig5's current BEELINE predictions (`Fig5_GRN/input/capsule/external/projects/grn_inference/BEELINE/outputs/scRNA-Seq/K562-3635genes-231metacells/`)
are byte-identical to `~/projects/grn_inference/BEELINE/outputs/...` from Slurm job 806370 (2026-01-19,
`run-k562-pipeline.slurm`). Its input was exported by
`chrono_seq/analyses/k562-tnfa.hg38.rep/notebooks/not_used/old/GRN_Inference.v1.ipynb` to `intmd/old/to_beeline/`,
from the k-means metacells of `not_used/old/Preprocess_GetTrajectories_MapFates.ipynb` and the gene filter of
`not_used/old/not_used/PrepareForTradeSeq.KMeansMetacellsOriginal.ipynb`.

## External dependencies

- `~/miniconda3/envs/chronopack`: steps 01, 03.
- `.conda_env/` (created here: conda-forge r-base 4.3, r-gam, r-data.table): step 02.
- `~/miniconda3/envs/BEELINE` (pandas 0.23, PyYAML): `generateExpInputs.py` and `BLRunner.py`.
- `singularity` and the BEELINE images in `~/projects/grn_inference/BEELINE/singularity_images/`.
