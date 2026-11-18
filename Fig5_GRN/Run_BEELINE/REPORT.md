# BEELINE comparators: metacell (current Fig5) vs r9 single-cell reruns: report (2026-10-20)

K562 + TNF-alpha GRN benchmark, panel 5b metrics (per-TF EPR and AUPRC ratio, TRRUST and K562 ENCODE ChIP).
Three sets of BEELINE predictions (GENIE3, GRNBOOST2, PPCOR, SINCERITIES) are compared, each scored with Fig5's own
code against the same ChronoPack_K562 predictions:

- **metacell**: the predictions Fig5 currently uses (BEELINE run of 2026-01-19).
- **r9 TFs + 1000** and **r9 TFs + 2000**: the reruns in this folder on the r9 single cells.

Fig5_GRN itself is unchanged. Scripts and run instructions: `README.md`.

**Update 2026-11-16.** Run_BEELINE no longer depends on `../old/` (user request): the metacell round and the checks
against Fig5 v1 were dropped and the steps now score with `..` (Fig5_GRN_v4)'s code. The **metacell** columns and the
"metacell reproduces Fig5 v1" checks below are the earlier results, kept as a record; they are no longer regenerated
(step 07's TRRUST-targets-lost table was removed). All r9 BEELINE numbers are unchanged by the rerun. ChronoPack_THP1
values in `output/rounds/` follow v4's TF rankings (CREB1 no longer scored), e.g. THP1 TRRUST P < 0.05 TFs + 2000
median EPR 4.77 -> 3.13 and AUPRCr 4.07 -> 4.26; this report's ChronoPack_K562 numbers are unchanged.

**EPR top-k fix (2026-10-22).** All EPRs in this report use the corrected top-k count (hits at ranks 1..k). Fig5's
original `evaluate_one`, which steps 05, 06 and 10 imported until 2026-10-22, also counted a hit at rank k + 1. The
steps now import Fig5_GRN_v3's fixed copy (`scripts/fig5_code.py`; see `../Fig5_GRN_v3/README.md`) and were rerun.
AUPRCr is unaffected. Values that changed in this report: TRRUST EPR of ChronoPack NFKB1 (6.36 -> 6.15) and of
metacell SINCERITIES JUN (1.79 -> 0.89), the metacell SINCERITIES TRRUST median EPR (0.85 -> 0.45) and three
matched-gene ChIP median EPRs (by 0.01). The metacell rows still reproduce Fig5 v1 (`old/Fig5_GRN`) apart from these
one-hit-lower top-k counts (2 of 60 per-TF rows, 18 of 209 matched-gene rows).

All eight BEELINE runs are complete. GENIE3 on r9 TFs + 1000 needed a second attempt: the first hung in dask
2 minutes in.

## Summary

1. **Under Fig5's scoring, the r9 reruns do worse on TRRUST and better on K562 ChIP than the metacell run.**
   - **TRRUST:** for RELA, NFKB1 and JUN, AUPRCr falls in 22 of the 24 method × input × TF comparisons. The
     exceptions are GRNBOOST2 on JUN (1.41 → 1.51 and 1.41). RELA drops to about random (0.97–1.08).
     ChronoPack's TRRUST lead therefore widens.
   - **K562 ChIP:** for GENIE3, PPCOR and SINCERITIES on RELA and JUN, AUPRCr rises from 1.11–1.14 to 1.34–1.65. On RELA ChIP
     EPR, SINCERITIES (2.47) and GENIE3 (2.39) on TFs + 2000 pass ChronoPack (2.25). ChronoPack still has the
     highest ChIP AUPRCr.
2. **Both shifts are gene coverage, not ranking quality.** Fig5 scores every method on all 19,338 protein-coding
   genes and treats unranked genes as random guesses.
   - **TRRUST:** BEELINE's 10%-of-cells expression filter removes most TRRUST targets of RELA, NFKB1 and JUN.
     These are sparse inflammatory genes: the median target is detected in 2.7–4.1% of cells. The r9 inputs hold 27–39 of
     about 300 RELA/NFKB1 targets, versus about 90 in the metacell input.
   - **ChIP:** the r9 high-variance genes are rich in ChIP-bound genes (18% of TFs + 2000 genes are RELA ChIP
     targets, versus 10% of the metacell genes).
3. **Matched-gene sensitivity analysis (Fig5's 04e, repeated on r9):**
   - **TRRUST:** ChronoPack's lead survives and grows. Median AUPRCr for ChronoPack versus the best BEELINE method is
     7.18 vs 1.89 on TFs + 1000 and 6.16 vs 1.57 on TFs + 2000; the metacell comparison was 3.40 vs 1.49. JUN is the
     exception: there, GENIE3 ties ChronoPack (3.49 and 3.30 vs 3.44 and 3.24).
   - **ChIP:** every method, ChronoPack included, is at random (AUPRCr 0.98–1.03), so all full-scoring ChIP
     differences on r9 are coverage.
4. **The r9 runs follow the BEELINE paper more closely than the metacell run.** They use single cells, BEELINE's
   10% filter and BEELINE's own TFs + N script and TF list. The metacell run used pooled cells, Seurat HVGs and a
   hand-picked TF list.

## Methods

### BEELINE inputs

| | metacell (current Fig5) | r9 TFs + 1000 | r9 TFs + 2000 |
|---|---|---|---|
| Cells | 231 k-means metacells (k = 20 per time point on 6 PCs, counts summed; pre-r9 processing) | 6,857 r9 prefiltered single cells | same |
| Expression | metacell counts, normalized to 1e4, log1p | r9 `log1p` layer (log1p CPM) | same |
| Gene filter | 5,000 Seurat HVGs (batch = trajectory cluster), minus time-tag-biased, MT, RP and ENSG genes | BEELINE: expressed in ≥ 10% of cells (3,611 of 11,866 genes) | same |
| Gene selection | the 3,635 filtered HVGs + 8 hand-picked TFs | BEELINE's `generateExpInputs.py -p 0 -t -n 1000`: all 366 expressed BEELINE TFs + 1,000 highest-variance non-TF genes | same, `-n 2000` |
| Forced TFs | ATF1, ATF2, CEBPB, CREB1, ELK1, ELK4, JUND, RELA | the 12 Fig5 TFs; 5 were added (RELA, JUN, ELK1, CREB1, FOS) | same |
| Genes in input | 3,643 | 1,371 | 2,371 |
| Protein-coding targets that can be scored | 2,918 | 1,363 | 2,338 |
| Pseudotime (SINCERITIES only) | time tag | time tag | time tag |
| Algorithm settings | BEELINE defaults; PPCOR pearson, pVal 0.01; SINCERITIES nBins 10 | same | same |
| Run | Slurm job 806370, 2026-01-19 | this folder, 2026-10-19/20 | same |

Runtimes (wall clock):

| Method | metacell | r9 TFs + 1000 | r9 TFs + 2000 |
|---|---|---|---|
| GENIE3 | 1 h 50 min | 2 h 53 min | 5 h 47 min |
| GRNBOOST2 | 11 min | 46 min | 1 h 28 min |
| PPCOR | 3.5 min | 28 s | 10 min |
| SINCERITIES | 35 min | 4 min | 1 h 2 min |

**How each run relates to the BEELINE paper** (Pratapa et al. 2020, scRNA-seq procedure):

- **Paper:** single cells; drop genes expressed in < 10% of cells; keep genes whose GAM-along-pseudotime P is < 0.01
  after Bonferroni; take significant TFs + the 500 or 1,000 most-varying significant genes.
- **r9 runs:** follow the paper except in two places, both agreed on 2026-10-19:
  - **No time-test filter (`-p 0`).** Only 277 genes pass it on these sparse single cells, too few for TFs + 1000.
    The test is still computed and stored (`output/prep/K562/GeneOrdering.csv`). It is applied only after the 10%
    filter, so it could not restore any of the genes discussed below.
  - **The Fig5 TFs are forced in.**
- **Metacell run:** departs from the paper in cell input (pooled cells), gene filter (no detection cutoff), gene
  ranking (Seurat HVGs rather than variance) and TF set (8 TFs rather than the BEELINE TF list). It also had no
  time test.

### Evaluation (Fig5 panel 5b)

The scoring is Fig5's own code, imported read-only. Only the BEELINE results directory changes (`scripts/05_fig5_per_tf_metrics.py`).

- **Edges:** BEELINE edges are loaded with Fig5's `count_top_predicted_tfs.load_beeline_method` (PPCOR scored by
  |corVal|, the other methods by EdgeWeight), restricted to reference TFs.
- **Ground truth:** TRRUST and K562 ENCODE ChIP (`tf_gene.1kbp_200bp.stringent.tsv`), targets restricted to 19,338
  protein-coding genes (MT, RP and ENSG removed).
- **Per-TF metrics:** `evaluate_one` over all 19,338 targets. A target the method does not rank is placed in a
  random tail (scored at the expected random hit rate).
  - **EPR:** precision at k = number of reference targets, divided by the random rate; capped at 25.
  - **AUPRCr:** average precision divided by the random rate.
- **TFs:** the v2 per-TF table includes FOS and ELK1. The medians below use the 5b strip-forest TF sets, which drop
  them for K562:
  - TRRUST: RELA, NFKB1, JUN, CREB1, ATF1, ATF2;
  - ChIP: RELA, JUN, CREB1, ATF1, ATF2 (there is no NFKB1 ChIP experiment).
- **Check:** rescoring the metacell outputs reproduces all 60 BEELINE rows of Fig5's per-TF table (maximum
  difference 4e-16).

### Sensitivity analysis (Fig5 04e, repeated)

`scripts/06_gene_universe_sensitivity.py` imports `Fig5_GRN/scripts/04e_gene_universe_sensitivity.py` unchanged and
swaps in each run's BEELINE predictions. It compares two target universes:

- **full:** the Fig5 universe above.
- **shared:** only the targets that both that run's BEELINE input and ChronoPack_K562 rank. On this set ChronoPack,
  GENIE3, PPCOR and SINCERITIES each rank every gene. GRNBOOST2 prunes edges, so it still ranks only part of it.
  This removes gene coverage and leaves only ranking order.
  - **metacell:** 2,478 genes (as in 04e).
  - **r9:** ChronoPack ranks every r9 BEELINE gene, so the shared set is the whole r9 input (1,363 and 2,338 genes).

Changes from 04e:
- only K562 methods are scored, so 04e's shared_THP1 universe is dropped;
- the oracle and coverage-null columns are computed as in 04e but not discussed here.

**Check:** the metacell run reproduces Fig5's 04e table (209 rows).

## Results under Fig5 scoring

### TRRUST, per TF

EPR:

| TF | Method | metacell | r9 TFs + 1000 | r9 TFs + 2000 |
|---|---|---|---|---|
| RELA | ChronoPack | 6.06 | 6.06 | 6.06 |
| | GENIE3 | 1.51 | 1.95 | 1.30 |
| | GRNBOOST2 | 3.46 | 0.46 | 0.33 |
| | PPCOR | 2.60 | 1.51 | 1.73 |
| | SINCERITIES | 2.81 | 0.87 | 0.65 |
| NFKB1 | ChronoPack | 6.15 | 6.15 | 6.15 |
| | GENIE3 | 6.36 | 3.18 | 3.18 |
| | GRNBOOST2 | 4.88 | 1.72 | 1.91 |
| | PPCOR | 2.54 | 1.91 | 2.12 |
| | SINCERITIES | 1.70 | 2.76 | 1.70 |
| JUN | ChronoPack | 6.26 | 6.26 | 6.26 |
| | GENIE3 | 6.26 | 5.37 | 5.37 |
| | GRNBOOST2 | 4.47 | 3.58 | 4.47 |
| | PPCOR | 1.79 | 2.68 | 0.89 |
| | SINCERITIES | 0.89 | 0.89 | 0.89 |

AUPRCr:

| TF | Method | metacell | r9 TFs + 1000 | r9 TFs + 2000 |
|---|---|---|---|---|
| RELA | ChronoPack | 4.04 | 4.04 | 4.04 |
| | GENIE3 | 1.40 | 1.13 | 1.08 |
| | GRNBOOST2 | 1.42 | 0.98 | 0.97 |
| | PPCOR | 1.67 | 1.08 | 1.04 |
| | SINCERITIES | 1.54 | 1.07 | 1.02 |
| NFKB1 | ChronoPack | 4.09 | 4.09 | 4.09 |
| | GENIE3 | 3.40 | 1.72 | 2.17 |
| | GRNBOOST2 | 2.84 | 1.26 | 1.55 |
| | PPCOR | 1.54 | 1.31 | 1.16 |
| | SINCERITIES | 1.44 | 1.16 | 1.08 |
| JUN | ChronoPack | 1.84 | 1.84 | 1.84 |
| | GENIE3 | 1.93 | 1.71 | 1.50 |
| | GRNBOOST2 | 1.41 | 1.51 | 1.41 |
| | PPCOR | 1.46 | 1.25 | 1.07 |
| | SINCERITIES | 1.53 | 1.16 | 1.08 |

- **EPR at random:** about 4.6 hits for RELA, 4.7 for NFKB1 and 1.1 for JUN. A single hit therefore moves JUN's EPR
  by about 0.9, so JUN EPR differences of 1–2 hits are noise.
- **JUN AUPRCr:** in the current Fig5, GENIE3 is above ChronoPack (1.93 vs 1.84). On r9 it is below it (1.71 on
  TFs + 1000, 1.50 on TFs + 2000).

### K562 ChIP, per TF

| TF | Method | EPR: metacell | EPR: r9 TFs + 1000 | EPR: r9 TFs + 2000 | AUPRCr: metacell | AUPRCr: r9 TFs + 1000 | AUPRCr: r9 TFs + 2000 |
|---|---|---|---|---|---|---|---|
| RELA | ChronoPack | 2.25 | 2.25 | 2.25 | 1.91 | 1.91 | 1.91 |
| | GENIE3 | 1.45 | 2.39 | 2.39 | 1.13 | 1.41 | 1.57 |
| | GRNBOOST2 | 1.24 | 1.24 | 1.22 | 1.10 | 1.12 | 1.10 |
| | PPCOR | 1.33 | 2.39 | 2.21 | 1.13 | 1.39 | 1.54 |
| | SINCERITIES | 1.22 | 2.38 | 2.47 | 1.11 | 1.45 | 1.65 |
| JUN | ChronoPack | 1.79 | 1.79 | 1.79 | 1.75 | 1.75 | 1.75 |
| | GENIE3 | 1.22 | 1.35 | 1.58 | 1.14 | 1.34 | 1.45 |
| | GRNBOOST2 | 1.05 | 1.05 | 1.05 | 1.06 | 1.09 | 1.08 |
| | PPCOR | 1.20 | 1.35 | 1.58 | 1.12 | 1.34 | 1.45 |
| | SINCERITIES | 1.20 | 1.34 | 1.58 | 1.14 | 1.34 | 1.46 |

**When ChIP EPR does not measure ranking:** EPR looks at the top k genes. When a method ranks fewer genes than k,
that top k is just every gene it ranked plus random filler, so EPR depends only on which genes it ranked, not on
their order. This applies to:

- **JUN ChIP (k = 4,123):** every BEELINE run in every input, since none ranks more than about 2,900 genes. That is
  why GENIE3, PPCOR and SINCERITIES all score exactly 1.58 on TFs + 2000.
- **RELA ChIP (k = 1,435):** GRNBOOST2 everywhere, and GENIE3, PPCOR and SINCERITIES on TFs + 1000.

ChronoPack ranks 10,201 genes, so its top k is always a real cutoff.

### Medians over the 5b TFs

| Benchmark | Method | EPR: metacell | EPR: r9 TFs + 1000 | EPR: r9 TFs + 2000 | AUPRCr: metacell | AUPRCr: r9 TFs + 1000 | AUPRCr: r9 TFs + 2000 |
|---|---|---|---|---|---|---|---|
| TRRUST | ChronoPack | 4.22 | 4.22 | 4.22 | 2.27 | 2.27 | 2.27 |
| | GENIE3 | 0.76 | 2.17 | 0.65 | 1.77 | 1.26 | 1.11 |
| | GRNBOOST2 | 2.92 | 1.09 | 0.16 | 1.42 | 1.12 | 1.06 |
| | PPCOR | 0.89 | 1.71 | 0.45 | 1.55 | 1.19 | 1.06 |
| | SINCERITIES | 0.45 | 0.43 | 0.32 | 1.53 | 1.14 | 1.07 |
| K562 ChIP | ChronoPack | 1.79 | 1.79 | 1.79 | 1.71 | 1.71 | 1.71 |
| | GENIE3 | 1.22 | 1.35 | 1.58 | 1.13 | 1.30 | 1.42 |
| | GRNBOOST2 | 1.05 | 1.05 | 1.05 | 1.05 | 1.09 | 1.08 |
| | PPCOR | 1.20 | 1.35 | 1.58 | 1.12 | 1.32 | 1.44 |
| | SINCERITIES | 1.19 | 1.34 | 1.58 | 1.11 | 1.31 | 1.45 |

## Why the results move: gene coverage

Reference targets inside each BEELINE input gene set (of all targets among the 19,338 genes):

| Benchmark | TF | Reference targets | In metacell input (2,918 genes) | In r9 TFs + 1000 (1,363) | In r9 TFs + 2000 (2,338) | Ranked by ChronoPack (10,201) |
|---|---|---|---|---|---|---|
| TRRUST | RELA | 299 | 90 | 27 | 37 | 148 |
| TRRUST | NFKB1 | 302 | 89 | 29 | 39 | 144 |
| TRRUST | JUN | 147 | 48 | 18 | 22 | 79 |
| K562 ChIP | RELA | 1,435 | 286 | 251 | 418 | 1,258 |
| K562 ChIP | JUN | 4,123 | 816 | 650 | 1,072 | 3,599 |

- **TRRUST, fewer targets in the gene set.** The r9 inputs have 2–3× fewer TRRUST targets than the metacell input,
  although they are 47–80% its size. Among RELA's 301 raw TRRUST targets:
  - 148 are in the r9 data at all;
  - 47 pass BEELINE's 10%-of-cells filter (median fraction of cells expressing them: 3.6%);
  - 37 survive the TFs + 2000 variance cut.

  NFKB1 (145 → 49 → 39) and JUN (80 → 25 → 22) follow the same pattern. Summing counts into metacells made these
  sparse genes detectable, so the metacell input kept them. The paper's procedure removes them.
- **ChIP, more targets in the gene set.** The share of RELA ChIP targets among the input genes is 7.4% across all
  genes, 9.8% in the metacell input and 17.9% in r9 TFs + 2000. Ranked at random, the top 1,435 of the TFs + 2000
  genes would hold about 256 RELA targets (EPR ≈ 2.4), which is what GENIE3, PPCOR and SINCERITIES score.

**Examples: TRRUST targets lost on r9 that ChronoPack ranks highly.** These were in the metacell input, are in
neither r9 input, and fall inside ChronoPack's top k (`output/fig5_metrics/trrust_targets_lost_vs_metacell.csv`).
ChronoPack gives every TF the same target order scaled by a TF weight, so its RELA and NFKB1 ranks agree to within
one position.

| Gene | ChronoPack rank (RELA / NFKB1) | Metacell GENIE3 rank for NFKB1 | Share of r9 cells expressing it | Why it is missing from r9 |
|---|---|---|---|---|
| BIRC3 | 9 | 17 | 10.3% | passes the 10% filter; variance rank 2,102 (cutoff 2,000) |
| TNFRSF9 | 16 | 1 | 6.7% | fails the 10% filter |
| IER3 | 26 | 252 | 11.7% | passes the filter; variance rank 2,694 |
| SERPINE1 | 37–38 | 111 | 4.5% | fails the filter |
| CXCL2 | 48–49 | 501 | 6.6% | fails the filter |
| TAPBP | 76–77 | 22 | 13.8% | passes the filter; variance rank 2,452 |
| PLAU | 91–92 | 36 | 2.7% | fails the filter |
| TNFAIP3 | 98–99 | 235 | 3.5% | fails the filter |
| BCL3 | 105–106 | 2 | 2.2% | fails the filter |
| NFKBIZ | 150–151 | 154 | 3.2% | fails the filter |
| CCL2 | 177–178 | 432 | 1.4% | fails the filter |
| EBI3 | 282–283 | 8 | 1.1% | fails the filter |

- **Totals:** across all ranks, 68 RELA, 69 NFKB1 and 38 JUN TRRUST targets were in the metacell input but are
  missing from r9 TFs + 2000.
- **JUN:** RELB (ChronoPack rank 61) and PLAU (91) are lost the same way.
- **NFKB1 GENIE3:** the metacell GENIE3 put several of these at the top of its NFKB1 ranking (TNFRSF9 #1, BCL3 #2,
  EBI3 #8). Losing them explains most of its EPR drop from 6.36 to 3.18.

## Sensitivity analysis: matched genes

Shared-set sizes:
- metacell: 2,478 genes;
- r9 TFs + 1000: 1,363 genes;
- r9 TFs + 2000: 2,338 genes.

TRRUST targets inside the shared set:
- RELA: 83, 27 and 37 (metacell, TFs + 1000, TFs + 2000);
- NFKB1: 79, 29 and 39;
- JUN: 42, 18 and 22.

### TRRUST, per TF

| TF | Method | EPR: metacell | EPR: r9 TFs + 1000 | EPR: r9 TFs + 2000 | AUPRCr: metacell | AUPRCr: r9 TFs + 1000 | AUPRCr: r9 TFs + 2000 |
|---|---|---|---|---|---|---|---|
| RELA | ChronoPack | 8.27 | 18.70 | 22.20 | 6.78 | 13.42 | 16.48 |
| | GENIE3 | 0.36 | 1.87 | 3.42 | 0.91 | 1.50 | 1.54 |
| | GRNBOOST2 | 1.80 | 0.00 | 0.00 | 1.50 | 0.85 | 0.88 |
| | PPCOR | 1.44 | 0.00 | 0.00 | 1.38 | 1.14 | 1.24 |
| | SINCERITIES | 0.72 | 0.00 | 0.00 | 1.08 | 1.02 | 1.07 |
| NFKB1 | ChronoPack | 8.74 | 14.59 | 19.98 | 7.00 | 10.91 | 15.70 |
| | GENIE3 | 6.35 | 8.10 | 10.76 | 4.93 | 5.92 | 9.32 |
| | GRNBOOST2 | 5.56 | 3.24 | 3.07 | 4.08 | 2.60 | 4.68 |
| | PPCOR | 1.59 | 3.24 | 6.15 | 1.27 | 2.78 | 2.03 |
| | SINCERITIES | 1.19 | 0.00 | 1.54 | 1.01 | 1.63 | 1.43 |
| JUN | ChronoPack | 4.21 | 4.21 | 4.83 | 2.60 | 3.44 | 3.24 |
| | GENIE3 | 5.62 | 12.62 | 14.49 | 1.55 | 3.49 | 3.30 |
| | GRNBOOST2 | 4.21 | 12.62 | 9.66 | 1.48 | 2.97 | 2.73 |
| | PPCOR | 1.40 | 0.00 | 0.00 | 0.89 | 1.35 | 0.97 |
| | SINCERITIES | 0.00 | 0.00 | 0.00 | 1.00 | 0.93 | 1.05 |

The counts behind these EPRs are small. With 22 JUN targets among 2,338 genes, random ranking puts 0.2 of them in
the top 22. GENIE3's 14.49 is 3 hits and ChronoPack's 4.83 is 1 hit.

### Medians over the 5b TFs

| Benchmark | Method | AUPRCr: metacell | AUPRCr: r9 TFs + 1000 | AUPRCr: r9 TFs + 2000 | EPR: metacell | EPR: r9 TFs + 1000 | EPR: r9 TFs + 2000 |
|---|---|---|---|---|---|---|---|
| TRRUST | ChronoPack | 3.40 | 7.18 | 6.16 | 2.11 | 2.10 | 2.42 |
| | GENIE3 | 1.23 | 1.89 | 1.57 | 0.18 | 0.93 | 1.71 |
| | GRNBOOST2 | 1.49 | 1.76 | 1.22 | 0.90 | 0.00 | 0.00 |
| | PPCOR | 1.23 | 1.45 | 1.27 | 0.70 | 0.00 | 0.00 |
| | SINCERITIES | 1.01 | 0.98 | 1.07 | 0.00 | 0.00 | 0.00 |
| K562 ChIP | ChronoPack | 1.17 | 1.03 | 1.03 | 1.17 | 1.00 | 1.01 |
| | GENIE3 | 1.05 | 0.99 | 0.98 | 1.06 | 0.98 | 0.98 |
| | GRNBOOST2 | 1.08 | 1.00 | 1.00 | 1.08 | 1.00 | 1.00 |
| | PPCOR | 1.02 | 1.00 | 1.00 | 1.01 | 1.00 | 1.00 |
| | SINCERITIES | 1.06 | 1.01 | 1.01 | 1.05 | 1.00 | 0.99 |

On r9 matched genes, the per-TF ChIP values for RELA and JUN are 0.76–1.10 for every method.

## Interpretation for Fig5

- **TRRUST.** ChronoPack's lead over the BEELINE methods holds on matched genes for every input and grows on r9.
  The exception is JUN, where GENIE3 matches ChronoPack's AUPRCr on very few targets (TFs + 1000: 3.49 vs 3.44;
  TFs + 2000: 3.30 vs 3.24). Switching
  Fig5 to an r9 run would widen the visible TRRUST gap for RELA and NFKB1. Part of that widening is coverage:
  BEELINE's 10% filter excludes the sparse NF-kB targets that ChronoPack still ranks.
- **K562 ChIP.** On matched genes neither ChronoPack nor the BEELINE methods rank ChIP targets better than chance.
  Every full-scoring ChIP difference, in any input, is therefore which genes each method covers. Switching to r9
  would narrow ChronoPack's visible ChIP lead. On RELA ChIP EPR, SINCERITIES and GENIE3 would pass ChronoPack.
- **Which BEELINE input to use.** The r9 runs are the defensible choice by protocol: the same cells as ChronoPack,
  and BEELINE's published scRNA-seq procedure apart from the time test and the forced TFs. The metacell run is not
  a BEELINE-standard input.
  - **TFs + 1000 vs TFs + 2000:** the two give the same qualitative picture. TFs + 2000 covers more genes and is
    closer in size to the metacell input.
  - **Possible reviewer objection:** "the comparators were denied the relevant genes." Two answers: the matched-gene
    results above, and, if needed, a sensitivity run with a looser expression filter (e.g. 2%), labelled as a
    deviation from BEELINE.

## Caveats

- TRRUST target counts per TF are small (147–302 overall; 18–39 inside the r9 gene sets). EPR in particular moves by
  large steps per hit.
- ChronoPack values come from Fig5's current r9 outputs and are identical in all three columns. Only the BEELINE
  predictions change.
- The time-test filter of the BEELINE paper was not applied. It runs only on genes that already pass the 10% filter
  and would have shrunk the r9 gene sets further (277 genes pass), so it would not restore the lost targets.

## Files

| Step | Script | Output |
|---|---|---|
| 05 | `scripts/05_fig5_per_tf_metrics.py` (`sbatch scripts/slurm/fig5_metrics.sbatch`) | `output/fig5_metrics/per_tf_metrics.csv`, `per_tf_comparison.csv` |
| 06 | `scripts/06_gene_universe_sensitivity.py` (`sbatch scripts/slurm/gene_universe.sbatch`) | `output/gene_universe/per_tf_metrics_by_universe.csv`, `median_by_universe.csv`, `positive_coverage_per_tf.csv` |
| 07 | `scripts/07_trrust_targets_lost.py` (Slurm, chronopack env, ~1 min) | `output/fig5_metrics/trrust_targets_lost_vs_metacell.csv` |

Logs: `output/logs/beeline_fig5_metrics_*`, `beeline_gene_universe_*`, `beeline_lost_targets_*`.

# All rounds, K562 and THP1 (added 2026-10-22)

Step 10 scores every BEELINE round with Fig5's 5b and 04e code (fixed EPR, rerun 2026-10-22); step 11 draws the figures
(`output/rounds/figures/`). Rounds: **10% filter** (BEELINE filter, all TFs + top-N variance genes; the runs
above), **no filter** (ChronoPack's expression filter, same `-p 0` selection), **P < 0.05** (no filter, BEELINE's
selection with an uncorrected GAM P < 0.05: significant TFs + the N most-varying significant genes). THP1 BEELINE is
scored like ChronoPack_THP1 (THP1 TF set: + REL, IRF3; FOS dropped). THP1 PPCOR is not applicable (README).
K562 GRNBOOST2, no filter, TFs + 2000 (18 h run) was added on 2026-10-22.

Checks: the K562 metacell and 10% rows reproduce step 06 (627 rows); the ChronoPack, Replogle and Yao rows equal
Fig5's published per-TF table (322 rows).

## Medians over the 5b TFs, TFs + 2000 (AUPRCr; EPR in brackets)

| Context, benchmark | Universe | ChronoPack | BEELINE, 10% filter | BEELINE, no filter | BEELINE, P < 0.05 |
|---|---|---|---|---|---|
| K562 TRRUST | Fig5 | 2.27 (4.22) | 1.06-1.11 (0.16-0.65) | 1.20-1.52 (0.89-3.85) | 1.21-1.30 (0.00-1.47) |
| K562 TRRUST | matched | see note | 1.07-1.57 | 1.04-1.98 | 1.05-1.89 |
| K562 ChIP | Fig5 | 1.71 (1.79) | 1.08-1.45 (1.05-1.58) | 1.08-1.46 (1.05-1.67) | 1.07-1.34 (1.04-1.48) |
| K562 ChIP | matched | 1.03-1.08 | 0.98-1.01 | 0.95-1.03 | 0.96-1.04 |
| THP1 TRRUST | Fig5 | 4.07 (4.77) | 1.50-1.63 (0.00-2.93) | 1.53-1.69 (0.00-3.04) | 1.49-2.20 (0.00-0.97) |
| THP1 TRRUST | matched | see note | 1.47-1.85 | 1.03-1.89 | 1.21-1.71 |
| THP1 ChIP | Fig5 | 1.50 (1.58) | 1.04-1.30 (1.03-1.41) | 1.02-1.33 (1.01-1.49) | 1.02-1.12 (1.02-1.14) |
| THP1 ChIP | matched | 0.93-0.95 | 0.99-1.03 | 0.99-1.03 | 1.00-1.02 |

Note: on matched genes ChronoPack's own value depends on the round's gene set. K562 TRRUST: 6.16 (10%), 2.93 (no
filter), 4.57 (P < 0.05); THP1 TRRUST: 4.81, 4.10, 3.39. With TFs + 1000, P < 0.05: K562 TRRUST matched ChronoPack
3.17 vs best BEELINE 2.56 (GENIE3); Fig5 universe 2.27 vs 1.26-1.40.

## Findings

1. **P < 0.05 removes most of the 10%-filter TRRUST penalty for K562 BEELINE, but ChronoPack's lead holds.** K562
   TRRUST median AUPRCr rises from 1.06-1.11 (10% filter) to 1.21-1.30, still well below ChronoPack's 2.27 and below
   the metacell run (1.42-1.77). On matched genes ChronoPack leads every BEELINE method in every round and context
   (K562 4.57 vs <= 1.89; THP1 3.39 vs <= 1.71).
2. **The K562 ChIP inflation from the 10% gene set shrinks.** BEELINE ChIP AUPRCr falls from 1.42-1.45 to 1.33-1.34;
   RELA ChIP EPR: GENIE3 2.39 -> 2.04, SINCERITIES 2.47 -> 2.21, PPCOR 2.21 -> 2.01, all now below ChronoPack (2.25).
3. **ChIP on matched genes is still random for every method** (median AUPRCr 0.92-1.11 in every r9 round, both
   contexts). All full-universe
   ChIP differences remain gene coverage.
4. **THP1 BEELINE is close to K562 BEELINE on TRRUST (1.5-2.2) and well below ChronoPack_THP1 (4.07).** SINCERITIES
   reaches 2.20 on P < 0.05 TFs + 2000, from REL (5.09) and IRF3 (3.89). Median EPR for THP1 BEELINE is 0 on TRRUST in most rounds
   because five THP1 TFs have few targets in the gene set; AUPRCr is the more informative metric there.
5. **The no-filter round sits between.** Adding ~600 lowly expressed TFs (no other gene changes) raises K562 GENIE3
   TRRUST to 1.52 in the Fig5 universe. Many TRRUST targets are themselves TFs: the share of each TF's TRRUST
   targets inside the gene set rises to 0.17-0.29. ChIP is unchanged.

Figures: `p05_strip_forest_{TRRUST,K562_ChIP}_{full,matched}_TFs{1000,2000}.pdf` (5b layout; K562 and THP1 blocks)
and `rounds_median_{AUPRCr,EPR}_TFs{1000,2000}.pdf` (medians per round, Fig5 and matched universes).
