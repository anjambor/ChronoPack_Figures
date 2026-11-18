"""Step 09 (run with the chronopack env python; ~1 min). Coverage of ChronoPack's top-k targets by the BEELINE paper's gene selection with a looser GAM cut-off.

Paper procedure (generateExpInputs.py -p P -t -n N, no -c): keep genes with GAM P < P; all significant TFs; the N
significant non-TF genes with the highest variance (all of them if fewer than N); plus the forced Fig5 TFs
(Run_BEELINE deviation 2). Emulated here from output/prep/<ctx>[-cpfilt]/GeneOrdering.csv, with the selection
re-checked against generateExpInputs.py's own logic (sort by P, filter, drop TFs, sort by variance, take N).
Compared with the current runs (10% filter or none, -p 0: all TFs + top-N variance).
ChronoPack target order: Fig5_GRN/output/chronopack/<ctx>/grn.target_scores.csv (identical to Fig5 v1's, which
this step read until 2026-11-16), restricted to the context's r9 genes.
"""
from pathlib import Path
import pandas as pd

RB = Path(__file__).resolve().parents[1]  # Run_BEELINE (inside Fig5_GRN since 2026-11-16)
FIG5 = RB.parent  # Fig5_GRN
tfs_all = set(pd.read_csv(RB / "input/human-tfs.csv").TF)
FORCED = ["RELA", "NFKB1", "JUN", "ATF2", "ELK1", "JUNB", "EGR1", "CREB1", "ATF1", "FOS", "REL", "IRF3"]
KS = [50, 100, 250, 500, 1000, 2000]


def select(go: pd.DataFrame, pcut: float, n: int, genes_in_expr: set) -> tuple[set, dict]:
    """generateExpInputs.py with -p pcut (no Bonferroni) -t -n n, then forced TFs appended."""
    g = go[go.index.isin(genes_in_expr)].sort_values(go.columns[0])
    if pcut:
        g = g[g[go.columns[0]] < pcut]
    tf = {t for t in tfs_all if t in genes_in_expr and t in g.index}
    rest = g.drop(index=list(tf))
    nontf = set(rest.index) if n > len(rest) else set(rest.sort_values(go.columns[1], ascending=False).index[:n])
    sel = tf | nontf
    forced = [t for t in FORCED if t not in sel and t in genes_in_expr]
    return sel | set(forced), {"n_sig_TFs": len(tf), "n_nonTF": len(nontf), "n_forced_added": len(forced),
                               "n_genes": len(sel) + len(forced)}


rows, sizes = [], []
for ctx in ("K562", "THP1"):
    ts = pd.read_csv(FIG5 / f"output/chronopack/{ctx}/grn.target_scores.csv")
    gf = pd.read_csv(RB / f"output/prep/{ctx}-cpfilt/gene_filter.csv", index_col=0).frac_cells_expressed
    r9 = set(gf.index)
    cp = [x for x in ts.target if x in r9]
    cp_nt = [x for x in cp if x not in tfs_all]
    sets = {}
    for n in (1000, 2000):
        sets[f"current 10% filter, TFs+{n}"] = set(pd.read_csv(RB / f"output/gene_selection_{ctx}-r9-sc-TFs{n}.csv").gene)
        for filt, tag in (("no filter", "-cpfilt"), ("10% filter", "")):
            go = pd.read_csv(RB / f"output/prep/{ctx}{tag}/GeneOrdering.csv", index_col=0)
            expr = set(go.index) | {t for t in FORCED if t in r9}  # ExpressionData.full.csv rows
            for p in (0.05, 0.01):
                s, info = select(go, p, n, expr)
                name = f"{filt}, P<{p}, TFs+{n}"
                sets[name] = s
                sizes.append({"context": ctx, "selection": name, **info})
        # the no-filter -p 0 selection must reproduce the existing cpfilt run (checks the emulation)
        go = pd.read_csv(RB / f"output/prep/{ctx}-cpfilt/GeneOrdering.csv", index_col=0)
        s0, _ = select(go, 0, n, set(go.index))
        assert s0 == set(pd.read_csv(RB / f"output/gene_selection_{ctx}-r9-sc-cpfilt-TFs{n}.csv").gene), "emulation mismatch"
    for k in KS:
        r = {"context": ctx, "k": k}
        for name, s in sets.items():
            r[name] = sum(x in s for x in cp[:k]) / k
            r[name + " [non-TF]"] = sum(x in s for x in cp_nt[:k]) / k
        rows.append(r)
df = pd.DataFrame(rows)
out = RB / "output" / "chronopack_target_coverage.csv"
df.to_csv(out, index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 60); pd.set_option("display.max_colwidth", 60)
print("emulation check: -p 0 selection reproduces the existing cpfilt gene sets: OK\n")
print(pd.DataFrame(sizes).to_string(index=False), "\n")
for n in (1000, 2000):
    cols = [c for c in df.columns if c.endswith(f"TFs+{n}")]
    print(f"== TFs+{n}: share of ChronoPack's top-k targets in the gene set")
    print(df[["context", "k"] + cols].round(3).to_string(index=False))
    print(f"== TFs+{n}: share of ChronoPack's top-k non-TF targets")
    print(df[["context", "k"] + [c + " [non-TF]" for c in cols]].round(3).to_string(index=False), "\n")
print("COVERAGE_PCUT_OK")
