# Fig5_GRN_v4: Fig5_GRN_v3 with ChronoPack's current TF rankings

> **Folder renamed 2026-11-17:** `Fig5_GRN_v4/` -> `Fig5_GRN/` (user request). Below, `Fig5_GRN_v4` means this folder; the
> earlier `Fig5_GRN` it was compared with (and other superseded versions) were deleted on 2026-11-16 and are in git
> history (commit 00ed526, `old/`). Script and output names keep their version tags.

Created 2026-10-22 (user request). This is `../Fig5_GRN_v3` with one change: **ChronoPack's TF rankings come from the
current ChronoPack** (`cp.tl.grn.rank_regulators` with its defaults since commit e0c89fb, 2026-10-22), so Fig5 is
consistent with the package. Everything else (BEELINE and motif-prior comparators, Replogle, Yao, references,
scoring, plotting, composite) is v3's.

- **Step 01** now computes the TF rankings (TNF for K562, TLR4 for THP1) from the OmniPath snapshot the ChronoPack
  analyses use (`input/chronopack/omnipath_snapshot_2025-10-21`, symlink to
  `ChronoPack/analyses/resources/`; no network access) instead of reading v3's saved 2026-09-30 r9 rankings.
  `output/chronopack/provenance.json` records the ChronoPack version / commit (e0c89fb, source unmodified), the
  snapshot and the `rank_regulators` defaults. `LIVE_OMNIPATH=1` uses the current OmniPath instead.
- **What the new defaults change** (ChronoPack CHANGELOG): conflicting OmniPath signs take the consensus sign, tied
  signs pass signal without direction, and cAMP / Ca2+ / DAG second-messenger edges are added. On the same snapshot
  with the old options, `rank_regulators` reproduces v3's rankings exactly (step 03), so every difference comes from
  the new defaults: Spearman over 501 TFs 0.97 (K562) / 0.94 (THP1), same top 10 in a slightly different order,
  TFs with positive priority 205 -> 195 (K562) and 217 -> 211 (THP1).
- **Unchanged:** target scores, and the target order of every TF in both networks (ChronoPack's edge weight is the TF's
  priority share times the target score). So every per-TF metric and PR curve of 5b-5d is v3's, except:
- **ChronoPack_THP1 CREB1 is no longer scored.** Its TLR4 signal is net inhibitory under the consensus signs
  (priority_act 0; rank 55 -> 214), so `infer` drops it, and as for any method's unpredicted TF (e.g. pySCENIC RELA)
  it is left out of that method's scores and plots (user decision 2026-10-22). CREB1 is still scored for every other
  method, including ChronoPack_K562.
- `input/`: `chronopack/<context>/tf_ranking.csv` symlinks removed (rankings are step 01 outputs); v3's `v2/` and
  `chronopack_notebook_reference/` not copied (step 03 compares with v3). `GENE_UNIVERSE_REPORT.md` stays in v3 (its
  analyses are per-TF on K562 and unaffected).

## v4 vs v3 (step 03, `output/comparison_to_v3.txt`, 2026-10-22: all checks passed)

| | v3 | v4 |
|---|---|---|
| 5b TRRUST ChronoPack_THP1 median EPR / AUPRCr (5b TFs) | 4.77 / 4.07 (9 TFs) | 3.13 / 4.26 (8 TFs) |
| 5b K562 ChIP ChronoPack_THP1 median EPR / AUPRCr | 1.58 / 1.50 (6 TFs) | 1.56 / 1.48 (5 TFs) |
| 5c Chrono \| THP1 macro-PR curves | with CREB1 | without CREB1 |
| 5e ChronoPack_K562 NES (RELA / NFKB1 / JUN / CREB1 / ATF2 / ATF1) | 156 / 145 / 136 / 116 / 124 / 51.7 | 154 / 145 / 134 / 107 / 124 / 50.1 |
| 5e fold over best comparator (same order) | 22.0 / 8.0 / 14.6 / 20.1 / 8.8 / 1.8 | 21.6 / 7.9 / 14.2 / 18.6 / 8.9 / 1.8 |
| 5e top-10,000 edge share, ChronoPack_K562 RELA / NFKB1 / JUN | 23.7% / 20.5% / 5.8% | 21.6% / 18.0% / 7.8% |

- Every other per-TF metric and PR curve is identical to v3 (120 groups expression, 159 with_prior); with_prior
  differs from v3 in the same way.
- 5e: the other methods' edge counts, ES and leading-edge ranks are identical; their NES move by Monte Carlo noise
  (max 0.23, median 0.03; with_prior max 0.30) because the permutation nulls of all methods are drawn from one random
  stream in method order with ChronoPack_K562 first (as described for the method sets below).
- 5d (RELA, NFKB1, JUN) is unchanged.
- Composite: `output/composite/Fig5_composite_v4.{png,pdf}` (TRRUST strip: 80 of 104 dots labelled).

Run: `sbatch scripts/slurm/fig5_pipeline.sbatch` (ends with `FIG5_V4_PIPELINE_OK`). Step 03 (the checks against
v3, `03_compare_to_v3.py`) was removed on 2026-11-16 when v3 moved to `../old/` (no dependencies on `old/`); its
last report, `output/comparison_to_v3.txt` (all checks passed), is kept as the record above. Side analyses (`scripts/slurm/tf_balanced_global.sbatch`) were not rerun for v4; step
04 reads the step 01 ranking.

## Comparator runs inside this folder (2026-11-16)

`Run_BEELINE/` (GENIE3, GRNBOOST2, PPCOR, SINCERITIES on the r9 single cells) and `Run_PriorGRN/` (pySCENIC,
CellOracle, Inferelator) moved into this folder (user request). `input/beeline_r9/` and `input/prior_grn/` are relative
symlinks into them; their project-local conda envs (`Run_*/.conda_env`) were recloned at the new paths with
`conda create --clone` (R gam 1.22.6; Python inferelator/anndata/pyarrow; `pip freeze` identical). Run them from their
own folders as before (e.g. `cd Run_BEELINE && sbatch scripts/slurm/rounds.sbatch`).

## with_prior v2 figures: shared gene universe, no Perturb-seq (2026-11-11, user request)

`output/with_prior/benchmark_outputs/figures/5{b,c,d,e}_v2/` (tables in `tables/shared_universe_v2/`) are the
with_prior panels with two changes, set by `FIG5_UNIVERSE=shared` (`fig5_paths.py`):

1. **Gene universe** = the genes every non-Perturb-seq method can rank: the BEELINE / motif-prior expression input
   (`Run_BEELINE` `K562-r9-sc-cpfilt-p05-TFs2000`, 2,331 genes incl. its 331 TFs as targets; identical for the
   motif-prior methods) that ChronoPack_K562 (10,202 targets) and ChronoPack_THP1 (9,481) also rank, within Fig5's
   protein-coding filter: **2,077 genes** (71 of the 2,331 are not protein-coding; 183 more are not in ChronoPack_THP1's
   fit). `shared_universe_genes.txt` lists them. Predictions, reference positives (TRRUST, K562 ChIP) and the random
   tail of `evaluate_one` are all cut to this universe (both `load_valid_targets` copies call
   `fig5_paths.restrict_universe`).
2. **Replogle2022_K562 and Yao2023_THP1 removed** from every panel (`fig5_paths.insert_prior` drops them); the THP1
   block keeps ChronoPack_THP1 only. 5e never included them.

The full-universe default is unchanged: 19,338 protein-coding targets, a choice of the evaluation code, not of the
Perturb-seq data. Run: `sbatch scripts/slurm/fig5_shared_universe.sbatch` (02, 02b, 02c in the build tree
`output/with_prior/shared_universe_build/`, then `02d_install_shared_universe_v2.py`: checks no Perturb-seq rows,
all targets in the universe, N = 2,077 for every scored per-TF row; moves outputs; ends `SHARED_UNIVERSE_V2_OK`).
`run_fig5_v4.sh` runs it after the full with_prior set (whose step 02 clears `output/with_prior/benchmark_outputs/`).
Axis limits are wider for these figures (`scripts/shared_universe_axes.py`); the K562 ChIP zoom is
`_xmax0p5_ymin0p3_ymax0p65` instead of the two y 0.2 zooms.

5b medians, full -> shared universe (AUPRCr; EPR):

| | TRRUST AUPRCr | TRRUST EPR | K562 ChIP AUPRCr | K562 ChIP EPR |
|---|---|---|---|---|
| ChronoPack_K562 | 2.27 -> 4.44 | 4.22 -> 1.66 | 1.71 -> 1.02 | 1.79 -> 1.04 |
| GENIE3 | 1.30 -> 1.52 | 0.97 -> 0.00 | 1.33 -> 0.96 | 1.48 -> 0.93 |
| SINCERITIES | 1.23 -> 1.01 | 0.45 -> 0.00 | 1.34 -> 1.02 | 1.48 -> 1.00 |
| GRNBOOST2 | 1.29 -> 1.73 | 1.47 -> 1.47 | 1.07 -> 1.02 | 1.04 -> 1.00 |
| PPCOR | 1.21 -> 1.15 | 0.00 -> 0.00 | 1.34 -> 1.01 | 1.48 -> 1.00 |
| pySCENIC | 1.51 -> 2.71 | 3.27 -> 9.40 | 1.03 -> 1.00 | 1.01 -> 1.00 |
| CellOracle | 1.80 -> 2.98 | 2.71 -> 1.66 | 1.21 -> 1.09 | 1.15 -> 1.03 |
| Inferelator | 1.12 -> 1.31 | 1.41 -> 0.00 | 1.06 -> 1.01 | 1.02 -> 1.00 |
| ChronoPack_THP1 | 4.26 -> 7.70 | 3.13 -> 3.32 | 1.48 -> 0.99 | 1.56 -> 0.99 |

Caveats: TRRUST positives per TF become small (ChronoPack_K562 rows: RELA 46, NFKB1 47, JUN 25, CREB1 13, ATF2 7,
ATF1 3), so EPR (top-k with k = positives) is coarse and median EPRs swing (e.g. 0 for most comparators); AUPRCr is
the steadier metric. K562 ChIP becomes random (~1) for every method, as in the 04e matched-universe analysis: the
full-universe ChIP gaps are gene coverage. pySCENIC's TRRUST median covers 3 TFs (no RELA / ATF2 / ATF1 regulons).

## Perturb-seq as a functional reference (2026-11-11, user request)

`output/with_prior/benchmark_outputs/figures/perturbseq_reference_v2/` (tables in `tables/perturbseq_reference_v2/`,
step `scripts/02e_perturbseq_reference.py`, `sbatch scripts/slurm/fig5_perturbseq_reference.sbatch`; also run by
`run_fig5_v4.sh` and `fig5_shared_universe.sbatch`). Replogle2022 (K562, 11 Fig5 TFs) and Yao2023 (THP1 LPS, 12 TFs)
are used as **references** (TF perturbation vs non-targeting DE, the TF's own gene excluded), not as comparator
methods; the methods are the with_prior networks without them, scored from the step 02 ranked-edge tables. Both gene
universes (user choice): the shared 2,077 genes and Fig5's 19,338 protein-coding genes, each cut to the genes measured
in the reference.

**Why rank-based first (user choice):** with a binary cut-off these references are nearly empty. **Correction
(2026-11-11):** the first version used the tables' adjusted P (Seurat `p_val_adj`), which is Bonferroni, and called it
FDR; binary positives are now Benjamini-Hochberg q < 0.05 per TF (`q_bh`; Bonferroni kept as `p_bonf`). BH positives
over all measured genes: Replogle 0-6 per TF except EGR1 87; Yao 8-167 (RELA 167, STAT1 162, NFKB1 79, FOS 60, JUNB
56). In the shared universe: Replogle EGR1 14, others 0-2; Yao RELA 45, NFKB1 17, STAT1 13, others 0-4.

- `perturbseq_rank_response_AUROC_{shared,full}`: per TF, the method's top 100 targets vs the other universe genes,
  AUROC on the DE response score -log10 P (direction-agnostic; 0.5 = random). `*` = GSEA permutation FDR < 0.05.
- `perturbseq_rank_ES_z_{shared,full}`: weighted GSEA ES of the same sets as a z-score against 1,000 random sets of
  the same size (0 = random). NES (random sets give |NES| ~ 1) and permutation P are in `perturbseq_rank_metrics.csv`
  (also K = 50 and 200).
- `perturbseq_fdr05_AUPRCr_{shared,full}`: BH FDR < 0.05 DEGs as positives, Fig5's evaluate_one AUPRCr (EPR in the
  table), TFs with >= 5 positives.
- Curves (`--curves-only` redraws them from the metric tables):
  `perturbseq_topk_response_curves_{shared,full}`: mean response percentile of each method's top k targets
  (k = 1-500, log axis; 0.5 = random), averaged over TFs; dashed where some TFs have fewer than k predictions.
  `perturbseq_gsea_curves_{Replogle_K562,Yao_THP1}_{shared,full}`: GSEA running enrichment of each method's top-100
  set, one panel per TF, over the 5-95% band and median of 200 random 100-gene sets. Weighted GSEA on -log10 P gives
  random sets clearly positive curves (most genes have near-zero weight), so read the methods against the band, not 0.
  Curve tables: `perturbseq_topk_response_curves.csv`, `perturbseq_gsea_running_curves.csv` (400-point grid + the
  extremes; the extreme matches the table ES), `perturbseq_gsea_random_set_band.csv`, `perturbseq_gsea_hits.csv`
  (the hit rug under each GSEA panel: one row per method, where its top-100 targets fall).
- scMINER Fig. 4-style panels (user request after reviewing scMINER / SCING / CausalBench):
  `perturbseq_roc_global_{shared,full}` (4c: all TF-gene pairs of a method's TFs pooled on edge weight, BH FDR < 0.05
  DEGs as labels), `perturbseq_roc_per_tf_{shared,full}` (4g: TFs with >= 5 DEGs; AUROC per panel; checked against
  Mann-Whitney U), `perturbseq_summary_{shared,full}` (4d: median over TFs of response AUROC and -log10 GSEA P, TF
  bootstrap 95% CI, TF dots, paired two-sided Wilcoxon vs the in-context ChronoPack). Tables `perturbseq_roc_auc.csv`,
  `perturbseq_roc_curves.csv` (501-point FPR grid), `perturbseq_summary_stats.csv`. Not done: scMINER's signed
  (positive / negative target, CRISPRa vs CRISPRi) GSEA, because ChronoPack's signs are per TF (`tf_ranking.csv`),
  not per edge, and the BEELINE methods are unsigned.
  Shared-universe AUROC (Yao, BH labels): global ChronoPack_K562 / THP1 0.84 / 0.84 vs comparators 0.47-0.60; per TF
  RELA 0.83 / 0.88 (best comparator CellOracle 0.67), NFKB1 0.74 / 0.74 (GENIE3 0.64), STAT1 0.63 / 0.59 (only
  ChronoPack predicts STAT1); Replogle EGR1 (no ChronoPack prediction) GENIE3 0.68. The global ROC pools TFs
  on edge weight, and ChronoPack's weight includes the TF priority, so it also rewards ranking RELA / NFKB1 (most of
  the DEGs) first. Full-universe ROC: DEGs outside a method's ranked genes fall in the tail after all ranked genes,
  so ChronoPack (ranks ~10k genes) gets AUROC < 0.5 for several TFs (e.g. THP1 ATF1 0.16, REL 0.27) from gene
  coverage, not ranking; use the shared universe for ROC comparisons. Replogle global ROC in the shared universe rests
  on 20 positive pairs (14 of them EGR1).
- Less strict / threshold-free per-TF ROC (user request):
  `perturbseq_roc_soft_per_tf_{shared,full}`: soft-label ROC, each gene weighted by 1 - local FDR from a
  beta-uniform mixture fit to the TF's P values (Pounds & Morris 2003; `soft`, `lfdr`, `pi0` per gene / TF); TPR / FPR
  = cumulative responder / non-responder mass; AUROC checked against a weighted Mann-Whitney U. TFs with responder mass
  < 2 (no detectable response: Replogle FOS, RELA; Yao ATF1, CREB1, ELK1, IRF3, JUN, STAT3) are omitted. The BUM
  spreads responder mass over hundreds of weakly shifted genes (e.g. Yao FOS 932, RELA 663 vs 41 / 148 BH hits), so
  soft AUROCs sit close to 0.5; only RELA and NFKB1 (Yao) separate methods (shared: ChronoPack 0.63-0.66 / 0.61-0.62 vs
  <= 0.58). `perturbseq_roc_threshold_sweep_{shared,full}`: AUROC with the top N responders by P as positives
  (N = 10-500); grey dashed line = BH count. ChronoPack's lead on Yao RELA / NFKB1 holds from N = 10 to 200 and fades to
  random by N = 500; Replogle stays at random for every N. Tables `perturbseq_roc_soft_auc.csv`,
  `perturbseq_roc_soft_curves.csv`, `perturbseq_roc_threshold_sweep.csv`.

Results (median over TFs):

| | Replogle K562 response AUROC (full / shared) | Yao THP1 response AUROC (full / shared) | Yao BH FDR AUPRCr (full / shared) |
|---|---|---|---|
| ChronoPack_K562 | 0.502 / 0.502 | 0.509 / 0.498 | 0.94 / 3.64 |
| GENIE3 | 0.478 / 0.482 | 0.529 / 0.517 | 1.00 / 2.56 |
| SINCERITIES | 0.492 / 0.489 | 0.509 / 0.500 | 1.02 / 0.94 |
| GRNBOOST2 | 0.492 / 0.492 | 0.522 / 0.518 | 1.27 / 2.05 |
| PPCOR | 0.496 / 0.497 | 0.514 / 0.500 | 1.03 / 1.02 |
| pySCENIC | 0.485 / 0.498 | 0.517 / 0.500 | 1.30 / 3.80 |
| CellOracle | 0.508 / 0.502 | 0.508 / 0.512 | 1.13 / 1.82 |
| Inferelator | 0.489 / 0.493 | 0.525 / 0.524 | 1.32 / 1.81 |
| ChronoPack_THP1 | 0.504 / 0.507 | 0.537 / 0.527 | 0.96 / 12.2 |

- Replogle K562: every method is at random (no significant GSEA cell); these knockdowns barely change expression in
  unstimulated K562, so this reference cannot separate methods.
- Yao THP1: weak signal, concentrated in RELA, NFKB1, FOS and JUNB. ChronoPack_THP1 (the only in-context network) has
  the highest median AUROC and the strongest single cells (RELA 0.78, NFKB1 0.69 full universe). Shared-universe BH
  AUPRCr medians cover RELA, NFKB1 and STAT1 (STAT1 only for ChronoPack). Response AUROC / GSEA columns are unaffected
  by the BH correction (no threshold).
- ChronoPack ranks targets in the same order for every TF, so its high RELA / NFKB1 scores reflect its target score
  ranking NF-kB-responsive genes first, not TF-specific target choice; its other TFs are near random.
- Responses include indirect effects (functional influence, not direct binding); K562 networks vs Yao and THP1 vs
  Replogle are cross-context comparisons.

---

**The rest of this README is v3's**, kept for the code and method description. Its numbers are v3's; where v4
differs, the table above applies. File names: `run_fig5_v4.sh`, Slurm jobs `fig5v4_*`, composite
`Fig5_composite_v4`; v3's step 03 (`03_compare_to_v2.py`) is replaced by `03_compare_to_v3.py`.

## v3: GRN benchmark with P < 0.05 BEELINE comparators and motif-prior methods

Created 2026-10-22. This is `../Fig5_GRN_v2` (panels 5b-5e) with two changes:

1. **BEELINE predictions** (GENIE3, SINCERITIES, GRNBOOST2, PPCOR) come from the r9 single-cell **P < 0.05** run
   of `Run_BEELINE` (`K562-r9-sc-cpfilt-p05-TFs2000`): no expression filter, the BEELINE paper's gene selection
   with an uncorrected GAM P < 0.05 (significant TFs + the 2,000 most-varying significant genes), Fig5 TFs forced;
   2,331 genes. v2 used the 10%-expression-filter run (`K562-r9-sc-TFs2000`, 2,371 genes); v1 used metacells.
2. **Two method sets.** Every cross-method panel is drawn twice:
   - `expression` (default): ChronoPack, the BEELINE methods, Replogle, Yao, i.e. v2's methods.
     Outputs in `output/benchmark_outputs/`.
   - `with_prior`: the same plus the K562 motif-prior methods of `Run_PriorGRN` on the same dataset:
     pySCENIC (cisTarget hg38 motif rankings), CellOracle and Inferelator (CellOracle's hg38 promoter base GRN).
     Outputs in `output/with_prior/benchmark_outputs/`.

Comparators stay K562-only, as in v1/v2 (user decision 2026-10-22): the THP1 block keeps ChronoPack_THP1 and Yao2023.
Everything else (ChronoPack r9 predictions, Replogle, Yao, ground truths, scoring and plotting code) is v2's, except:

3. **EPR top-k fix** (2026-10-22, user request; see below).
4. **No sEPR figures** (2026-10-22, user request): the sEPR strip-forests and sEPR violins (v4 `generate_6b`) are no
   longer made. 5b keeps the EPR / AUPRCr strip-forests (v2-style and label-placed v3) and the EPR@500 ChIP variant.
5. **Figure set changes** (2026-10-22, user request): the 5c micro-PR figures are no longer drawn (their table,
   `system_broad_filtered_publication_v4_micro_pr_matched_metrics.csv`, is still written); the 5e composite is split
   into its two panels, `5e/regulator_prioritization_v2_nes_heatmap_signaling_only_tf_rows` (rank-enrichment heatmap +
   fold over the best comparator) and `5e/regulator_prioritization_v2_top10000_composition_stacked_signaling_only`
   (share of the top 10,000 edges by TF).

v1 and v2 now live in `../old/Fig5_GRN/` and `../old/Fig5_GRN_v2/` (moved 2026-10-22); they keep the original code.

## EPR top-k fix

`benchmark_grn_systems.evaluate_one` counts the true targets in the top k of a TF's ranking (k = the TF's number of
reference targets) with `np.searchsorted` on 0-based hit positions. The capsule, v1 and v2 used `side="right"`, which
also counts a hit sitting at rank k + 1; v3 uses `side="left"` (hits at ranks 1..k). This changes TP_at_k / EPR and
the TP / precision / recall points of the PR curves; AP and AUPRCr are computed separately and do not change.

Checks: step 02c recounts the top k of every per-TF row with plain slicing and requires equality (161 rows in
with_prior, max |delta| 5e-13); the old rule would change 8 of them. In the 5b panels only two values move, both
TRRUST NFKB1: ChronoPack_K562 6.36 -> 6.15 and PPCOR 2.12 -> 1.91. The others are 5b-excluded or THP1 ChIP rows
(ChronoPack_K562 FOS, Replogle ATF2, ChronoPack_THP1 JUN / CREB1 / ATF1 / ATF2, each by 0.01 or less). No 5b
median changes. Step 03 accepts v2 differences only of this kind (one hit fewer, never on PR expectation rows): 7 of
62 non-BEELINE per-TF rows and 41,143 of 146,010 PR-curve points. `Run_BEELINE` imports v1's code, so its EPRs keep
the old rule (1 BEELINE row differs).

Why the P < 0.05 selection: `Run_BEELINE/README.md` and `REPORT.md` ("All rounds"). The motif-prior methods
(and why SCENIC+ and Pando are not included: no ATAC-seq for these experiments): `Run_PriorGRN/README.md`.

## Run

From this directory, on a compute node (~10 min):

```bash
sbatch scripts/slurm/fig5_pipeline.sbatch     # 01 once; 02, 02b, 02c per method set; 03 checks
```

A successful run ends with `FIG5_V3_PIPELINE_OK`; logs in `output/logs/`. A single method set can be rebuilt with
`FIG5_METHOD_SET=with_prior python scripts/02_reproduce_figures.py` (chronopack env). Side analyses
(`STEPS="04a 04 04b 04c 04d 04e" sbatch scripts/slurm/tf_balanced_global.sbatch`) run on the expression set only
(user decision); of them only 04e has been rerun for v3 (both method sets, 2026-10-22).

| Step | Script | What |
|---|---|---|
| 01 | `01_run_chronopack_grn.py` | ChronoPack GRN inference on r9 data (identical to v1/v2), shared by both method sets |
| 02 | `02_reproduce_figures.py` | benchmark tables and panels 5b-5e for `FIG5_METHOD_SET` |
| 02b | `02b_plot_expr_curves_1lineage.py` | 5d expression curves (no comparator methods) |
| 02c | `02c_k562_chip_epr_at_500.py` | 5b K562 ChIP strip-forests with EPR at k = 500 |
| 03 | `03_compare_to_v2.py` | checks -> `output/comparison_to_v2.txt` (below) |
| 04e | `04e_gene_universe_sensitivity.py` (`STEPS=04e [FIG5_METHOD_SET=with_prior] sbatch scripts/slurm/tf_balanced_global.sbatch`) | matched-gene universes, coverage null and oracle -> `output/[with_prior/]gene_universe/` (below) |

No notebook was made for v3 (v2's `Fig5_v2.ipynb` was not copied).

## Code changes from v2

- **`scripts/fig5_paths.py`:** `BEELINE_RESULTS_DIR = input/beeline_r9/K562-r9-sc-cpfilt-p05-TFs2000`;
  `PRIOR_RESULTS_DIR = input/prior_grn/K562-r9-sc-cpfilt-p05-TFs2000`; `FIG5_METHOD_SET` -> `METHOD_SET`,
  `WITH_PRIOR`; `BENCH_ROOT = output/with_prior` for with_prior; `PRIOR_LABELS`, `PRIOR_COLORS` (purple / lavender /
  brown); `insert_prior(methods, after)` puts the three prior methods after the last BEELINE method of a list when
  with_prior is on (a no-op for expression); `V2_DIR` replaces `V1_DIR`.
- **Method lists** (`capsule_lib/`): every list that v2 hard-codes now goes through `insert_prior`, and the label /
  colour dicts gain the prior methods: `benchmark_grn_systems` (`METHOD_ORDER`, `K562_NATIVE_METHODS`, colours,
  labels), `count_top_predicted_tfs`, `compare_regulator_prioritization`, `plot_regulator_prioritization_v2`, the
  publication-figure modules (`ALL_METHOD_SHORT`, `VIOLIN_METHOD_LABELS`, PR labels), `generate_publication_figures_v4`
  (`K562_METHODS_ALL`), micro-PR and 5d PR-curve modules, and 5e (`PANEL_B_COMPARATOR_METHODS`, NES
  `K562_METHODS` / `COMPARATOR_METHODS`).
- **PR panels** (5c macro/micro, 5d per-TF PR curves and heatmaps) show GENIE3 as the only BEELINE method, as in
  v1/v2; with_prior adds the three prior methods after GENIE3.
- **Loaders:** `load_prior_method` (rankedEdges.csv, normalised like BEELINE) in `count_top_predicted_tfs` and
  `benchmark_grn_systems`; step 02's ChronoPack loader adds the prior methods for with_prior.
- **Layouts that assumed 8 methods:** strip-forest row positions, figure height and the K562 / THP1 block labels
  come from `benchmark_grn_systems.strip_rows` (v2's exact values for 8 rows); the violins, the 5d heatmaps (height and
  absolute margins) and the K562 / THP1 row groups scale with the number of rows; the 5e "BASELINE METHODS" header
  spans however many comparators there are. With the expression set every layout value is v2's.
- **Step 03:** `03_compare_to_v1.py` replaced by `03_compare_to_v2.py`.
- **EPR fix:** `side="left"` in `evaluate_one`; 02c also writes
  `tables/publication_figures/system_broad_filtered_publication_v2_per_tf_epr_check.csv` (independent top-k recount);
  `generate_publication_figures_v4_filtered_tfs.main` no longer calls `generate_6b` (sEPR figures).
- **04e:** writes to `output/[with_prior/]gene_universe/` (per method set).
- **Figure set:** `generate_publication_figures_v4_filtered_tfs.main` builds the micro-PR table but no longer plots it;
  `plot_regulator_prioritization_v2_composite_stacked.main` saves the left and right panels as separate figures
  (`OUT_STEM_LEFT`, `OUT_STEM_RIGHT`; segment labels need a share >= 0.025 instead of 0.018 so short segments stay
  unlabelled).
- **Renames:** `run_fig5_v3.sh`, Slurm job names `fig5v3_*`; `nb_helpers.py` and `execute_notebook.sbatch` removed.

## Inputs (`input/`)

- `beeline_r9/K562-r9-sc-cpfilt-p05-TFs2000`: symlink into `Run_BEELINE/beeline/outputs/scRNA-Seq/`.
- `prior_grn/K562-r9-sc-cpfilt-p05-TFs2000`: symlink into `Run_PriorGRN/output/networks/` (`PYSCENIC`,
  `CELLORACLE`, `INFERELATOR` `rankedEdges.csv`).
- `capsule/`, `chronopack/`, `chronopack_notebook_reference/`: as in v2.
- `v2/`: v2's step 02 tables (`tables/*.csv` except the 360 MB ranked-edge table, and
  `tables/publication_figures/*.csv`), copied from `../Fig5_GRN_v2/output/` on 2026-10-22, for step 03.

## Checks (step 03, 2026-10-22: all passed)

- ChronoPack K562 and THP1 outputs are identical to the r9 GRNInference notebook's.
- Expression set vs v2: every non-BEELINE per-TF metric and PR curve is identical apart from the EPR fix's
  one-hit-fewer top-k counts (micro-PR and top-k shares identical).
- BEELINE per-TF rows equal `Run_BEELINE/output/rounds/per_tf_metrics.csv` (r9_p05_TFs2000, full universe; 44 of
  44 rows) in n_pred, TP_emitted and AUPRCr (EPR: old rule there, 1 row differs).
- 02c: every per-TF TP_at_k equals an independent plain top-k recount.
- with_prior vs expression: every per-TF metric and PR curve of the expression-set methods is identical; only the
  three prior methods are added.
- Not checked, expected: ChronoPack_THP1's rank-enrichment NES differs between the method sets by up to 14 (the
  shared-RNG effect described in v1's README: permutation nulls are drawn from one stream in method order, and the
  prior methods come before ChronoPack_THP1). 5e shows only the K562 NES, which is identical.

## Results

Medians over the 5b TFs (TRRUST: RELA, NFKB1, JUN, CREB1, ATF1, ATF2; K562 ChIP: the same without NFKB1; THP1
methods add REL and IRF3 and keep ELK1). A method is scored only on TFs it predicts for.

| Method | TRRUST EPR | TRRUST AUPRCr | ChIP EPR | ChIP AUPRCr | AUPRCr in v2 (TRRUST / ChIP) |
|---|---|---|---|---|---|
| ChronoPack_K562 | 4.22 | 2.27 | 1.79 | 1.71 | 2.27 / 1.71 |
| GENIE3 | 0.97 | 1.30 | 1.48 | 1.33 | 1.11 / 1.42 |
| SINCERITIES | 0.45 | 1.23 | 1.48 | 1.34 | 1.07 / 1.45 |
| GRNBOOST2 | 1.47 | 1.29 | 1.04 | 1.07 | 1.06 / 1.08 |
| PPCOR | 0.00 | 1.21 | 1.48 | 1.34 | 1.06 / 1.44 |
| pySCENIC (3 / 2 TFs) | 3.27 | 1.51 | 1.01 | 1.03 | - |
| CellOracle | 2.71 | 1.80 | 1.15 | 1.21 | - |
| Inferelator | 1.41 | 1.12 | 1.02 | 1.06 | - |
| Replogle2022_K562 | 0.42 | 0.90 | 1.87 | 1.70 | unchanged |
| ChronoPack_THP1 | 4.77 | 4.07 | 1.58 | 1.50 | unchanged |
| Yao2023_THP1 | 0.00 | 1.70 | 1.39 | 1.34 | unchanged |

- **BEELINE, P < 0.05 vs 10% filter:** TRRUST AUPRCr rises (1.06-1.11 -> 1.21-1.30) and K562 ChIP falls
  (1.42-1.45 -> 1.33-1.34), as in `Run_BEELINE/REPORT.md`. On RELA ChIP EPR, GENIE3 (2.04), SINCERITIES (2.21)
  and PPCOR (2.01) are now below ChronoPack (2.25); in v2 SINCERITIES and GENIE3 were above it.
- **Motif-prior methods:** CellOracle is the strongest comparator on TRRUST (median AUPRCr 1.80) but stays below
  ChronoPack_K562 (2.27); on ChIP all three are near random (1.03-1.21). pySCENIC has no RELA, ATF1 or ATF2 regulon
  on this dataset (cisTarget pruning), so its medians cover only 3 TRRUST / 2 ChIP TFs and are not comparable with
  the others'; the missing TFs are blank in the 5d heatmaps and 5e.
- **5c micro-PR (EPR / AUPRCr):** TRRUST: ChronoPack_K562 5.95 / 3.79, GENIE3 4.50 / 1.57 (v2: 2.90 / 1.37),
  pySCENIC 3.78 / 2.05, CellOracle 4.65 / 1.85, Inferelator 1.87 / 1.05. K562 ChIP: ChronoPack_K562 1.13 / 1.24,
  GENIE3 1.37 / 1.32, pySCENIC 1.01 / 1.03, CellOracle 1.13 / 1.19, Inferelator 1.03 / 1.05.
- **5e regulator prioritization:** the best comparator per TF is a BEELINE method in both sets; ChronoPack_K562's
  fold over it is RELA 22.0x, CREB1 20.1x, JUN 14.6x, ATF2 8.8x, NFKB1 8.0x, ATF1 1.8x. The prior methods' NES are
  at most 5.0 (Inferelator, JUN). Share of the top 10,000 edges from the signalling TFs: pySCENIC 4.5%, CellOracle
  3.9%, Inferelator 3.4% (ChronoPack 50.0% K562, 39.6% THP1).

## Does the target-gene universe drive the differences? (04e, 2026-10-22)

Full write-up of this exploratory analysis (gene sets, enrichment, random completion, 04e / 04f / 04g, JUN):
`GENE_UNIVERSE_REPORT.md`.

ChronoPack ranks ~10,200 targets per TF; the BEELINE methods rank the 2,260 genes of their input (GRNBOOST2 and the
motif-prior methods fewer). 5b scores everyone on all 19,338 protein-coding genes, with unranked genes filled at
random. 04e separates gene coverage from ranking (`output/with_prior/gene_universe/`, all 11 methods):
- **coverage null**: the score of the method's own ranked genes in random order (what knowing only *which* genes a
  method ranks is worth);
- **shared_K562**: only the 2,260 BEELINE-input genes (ChronoPack ranks all of them), so every method is scored on
  the same genes;
- **full_cp_hvg**: full universe, but ChronoPack's edges cut to the BEELINE genes.

Median AUPRCr over the 5b TFs (EPR in brackets):

| Method | TRRUST full | TRRUST coverage null | TRRUST shared | ChIP full | ChIP coverage null | ChIP shared |
|---|---|---|---|---|---|---|
| ChronoPack_K562 | 2.27 (4.22) | 0.99 | 4.57 | 1.71 (1.79) | 1.60 | 1.03 |
| ChronoPack_K562, BEELINE genes only (full_cp_hvg) | 2.30 (4.33) | | | 1.34 (1.48) | | |
| GENIE3 | 1.30 | 1.27 | 1.55 | 1.33 | 1.34 | 0.96 |
| SINCERITIES | 1.23 | 1.30 | 1.05 | 1.34 | 1.33 | 1.04 |
| GRNBOOST2 | 1.29 | 1.20 | 1.89 | 1.07 | 1.06 | 1.02 |
| PPCOR | 1.21 | 1.33 | 1.09 | 1.34 | 1.34 | 1.01 |
| pySCENIC (3 / 2 TFs) | 1.51 | 1.43 | 2.28 | 1.03 | 1.03 | 1.00 |
| CellOracle | 1.80 | 1.42 | 3.36 | 1.21 | 1.19 | 1.10 |
| Inferelator | 1.12 | 1.15 | 1.28 | 1.06 | 1.06 | 1.02 |

- **TRRUST: the gene universe does not explain ChronoPack's lead.** ChronoPack's coverage null is 0.99, so its
  2.27 comes from ranking, and cutting it to the BEELINE genes leaves it unchanged (2.30). On the shared genes its
  lead grows (4.57 vs at most 1.89 for BEELINE and 3.36 for CellOracle). The BEELINE methods' full-universe scores
  equal their coverage nulls (1.21-1.30 vs 1.20-1.33): their TRRUST signal there is which genes they include, not their
  order. CellOracle is the one comparator whose ranking adds clearly (1.80 vs 1.42).
- **K562 ChIP: the gene universe explains essentially all of ChronoPack's lead.** ChronoPack's coverage null is
  already 1.60 of its 1.71: ChIP-bound genes are enriched among the genes it ranks. Cut to the BEELINE genes it
  scores 1.34, the same as GENIE3 / SINCERITIES / PPCOR (1.33-1.34), and on the shared genes every method is at
  random (0.96-1.10). Replogle (1.70 vs null 1.73) and Yao (1.34 vs 1.33) behave the same way.
- **EPR on the shared genes** rests on very few positives per TF (TRRUST: 3-56 within 2,260 genes; CREB1 / ATF1 / ATF2
  mostly 0 hits), so its medians (e.g. ChronoPack 1.26, pySCENIC 7.53 over 3 TFs) are noisy; AUPRCr is the
  reliable comparison there. Per-TF values: `per_tf_metrics_by_universe.csv`.

## Two-tier completion: comparators given ChronoPack's candidate genes (04f, 2026-10-22)

`scripts/04f_two_tier_completion.py` (`STEPS=04f FIG5_METHOD_SET=with_prior sbatch scripts/slurm/tf_balanced_global.sbatch`)
completes each K562 comparator's ranking in two tiers: its own ranked targets in its own order, then the rest of
ChronoPack_K562's ~10,200 candidate genes in random order (20 draws), then the usual random tail over the remaining
universe. Every comparator thus gets ChronoPack's gene-set advantage. Medians over the 5b TFs
(`output/with_prior/two_tier_completion/`):

| Method | TRRUST AUPRCr 5b -> two-tier | K562 ChIP AUPRCr 5b -> two-tier | K562 ChIP EPR 5b -> two-tier |
|---|---|---|---|
| ChronoPack_K562 | 2.27 | 1.71 | 1.79 |
| GENIE3 | 1.30 -> 1.24 | 1.33 -> 1.70 | 1.48 -> 1.74 |
| SINCERITIES | 1.23 -> 1.16 | 1.34 -> 1.71 | 1.48 -> 1.75 |
| GRNBOOST2 | 1.29 -> 1.26 | 1.07 -> 1.63 | 1.04 -> 1.66 |
| PPCOR | 1.21 -> 1.16 | 1.34 -> 1.71 | 1.48 -> 1.74 |
| pySCENIC (3 / 2 TFs) | 1.51 -> 1.48 | 1.03 -> 1.61 | 1.01 -> 1.65 |
| CellOracle | 1.80 -> 1.70 | 1.21 -> 1.71 | 1.15 -> 1.71 |
| Inferelator | 1.12 -> 1.03 | 1.06 -> 1.61 | 1.02 -> 1.66 |

- **K562 ChIP:** ChronoPack's lead disappears: GENIE3, SINCERITIES, PPCOR and CellOracle reach 1.70-1.71 median
  AUPRCr (ChronoPack 1.71). Per TF only RELA keeps a margin (ChronoPack 1.91 vs 1.68-1.82).
- **TRRUST:** the comparators get slightly worse (ChronoPack's extra genes are depleted of TRRUST targets, 0.74x),
  so ChronoPack's lead (2.27 vs at most 1.70) is unaffected.

## Composite figure (step 05, 2026-10-22)

`scripts/05_composite_figure.py` (`sbatch scripts/slurm/composite.sbatch`; also the last step of `run_fig5_v3.sh`)
rebuilds the user's metacell-era composite `ref/FigureComposite.v6.png` (5,000 x 4,111 px, 300 dpi) from this folder's
expression method set: `output/composite/Fig5_composite_v3.{png,pdf}`.
- **Layout:** each subfigure is placed in the bounding box its counterpart occupies in the reference (`REF_BOXES`,
  measured from the reference's non-white pixels; `output/composite/panel_boxes.csv` lists reference vs placed boxes,
  all within a few px). Panel letters at the reference positions and height. Panel (a), the method schematic, is
  left blank (no assets here). `overlay_with_reference.png` blends the two for checking.
- **Subfigures** are redrawn from this folder's tables in the reference's style, not copied from the 5b-5e outputs
  (which are unchanged); each is rendered at 600 dpi with its figure size adjusted to the reference box's aspect ratio
  (`output/composite/panels/`):
  - (b) 5b strip-forests: y-axis lines, bracketed K562 / THP-1 groups, bold 24-pt titles and bold axis labels, larger
    method / TF / tick labels, Median legend in the ChIP panel only. TF labels in the reference's hand-placed style
    (`scripts/strip_label_layout.py`): each label right beside its dot (no leader lines), dots nudged vertically within
    their row to make room, tied values (e.g. EPR = 0) as compact stacks labelled at their ends, and a label left out
    when no clear, unambiguous spot exists (TRRUST 82 of 106 dots labelled, K562 ChIP 81 of 82). The layout is
    optimised from v3's values (seeded annealing, so reproducible); `panels/b_{trrust,chip}_labels.csv` lists which
    TFs are labelled. `COMPOSITE_STRIP_LABELS=all` restores v3's every-TF placement;
  - (c) mean macro-PR: TRRUST with the reference legend (Chrono | THP1, Chrono | K562, GENIE3 | K562, Yao | THP1,
    Replogle | K562), K562 ChIP zoom (x 0-0.5, y 0.2-0.6) with the title inside the axes; thicker lines;
  - (d) 5e rank-enrichment heatmap without the fold column, columns labelled at a 40-degree angle, colour bar beside
    the heatmap with ticks and label on its right;
  - (e) 5d TRRUST heatmap with "THP-1"; per-TF TRRUST PR curves without legend (y up to 1.0 for RELA / NFKB1, 0.1 for
    JUN); the ChronoPack expression curves (02b's PNGs).
- All text is Liberation Sans (Arial's metrics, registered from the env's R package fontLiberation as in 02b), black.
- Panels are cropped to their non-white pixels before placement (the reference boxes were measured that way), so
  each panel's content fills its reference box. TRRUST strip plot: EPR axis from -0.75 (5b: -0.5), so the EPR = 0
  dots sit as far from the y axis as in the reference; both strip plots use the reference's title gap and tick
  length / label gap (axes frames within ~10 px of the reference's).
- **PNG** = raster canvas (5,000 x 4,111 px, 300 dpi). **PDF** = vector: each panel's own PDF (transparent
  background, TrueType text), cropped like the PNG and placed in the same box with pypdf; no raster images
  (panel d's heatmap is drawn with pcolormesh and its colour bar unrasterized). pypdf is pure Python, installed in
  the gitignored `Fig5_GRN_v3/.pylib` (`python -m pip install --no-deps --target .pylib pypdf`).
- Known difference: label choices follow rules, not the reference's exact picks (the values changed).

## Figures

Under each `benchmark_outputs/figures/` (expression and with_prior):
- **5b:** EPR / AUPRCr strip-forests for TRRUST and K562 ChIP, plus the K562 ChIP EPR@500 variant. Each comes in two
  versions with identical data and statistics: `_v2` (labels beside each dot, may overlap; PDF + PNG +
  `_summary.csv` / `_clipped_points.csv`) and `_v3` (overlapping dots spread within their row, labels placed clear of
  each other with leader lines; PDF only). Both have the same size (v3's TRRUST was drawn taller until 2026-10-22,
  now at v2's height with refined label placement; a few labels in the densest rows still touch). No sEPR strip-forests or violins.
- **5c:** mean macro-PR curves (TRRUST, K562 ChIP, and the K562 ChIP zooms). No micro-PR figures.
- **5d:** per-TF TRRUST PR curves (RELA, NFKB1, JUN), the TRRUST heatmaps and the ChronoPack expression curves.
- **5e:** the two former composite panels as separate figures (rank-enrichment heatmap with fold over the best
  comparator; top-10,000 composition bars).

The with_prior strip-forests have 11 rows (pySCENIC, CellOracle, Inferelator after PPCOR); the PR panels and 5d
heatmaps have 8.
