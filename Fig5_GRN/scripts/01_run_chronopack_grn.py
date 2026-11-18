#!/usr/bin/env python
"""ChronoPack GRN inference for K562 (TNF) and THP1 (TLR4) on ChronoPack r9 data.

Same calls as ChronoPack/analyses/grn/GRNInference.ipynb (commit e0c89fb): rank_regulators ->
score_targets -> infer, written per context to output/chronopack/<context>/.

Fig5_GRN_v4: the TF rankings are computed with the installed ChronoPack's rank_regulators defaults from the
OmniPath snapshot the analyses use (fig5_paths.OMNIPATH_SNAPSHOT; no network access), instead of v3's saved
2026-09-30 r9 rankings. --live-omnipath queries the current OmniPath instead of the snapshot.
Also writes output/chronopack/provenance.json (ChronoPack version and commit, snapshot, ranking options).
"""
from __future__ import annotations

import argparse
import inspect
import json
import subprocess
from pathlib import Path

import chronopack as cp

from fig5_paths import CHRONOPACK_INPUT_DIR, CHRONOPACK_OUTPUT_DIR, CONTEXTS, OMNIPATH_SNAPSHOT

GIT = "/new-stg/home/alex/miniconda3/envs/metacellR/bin/git"  # no git on PATH on the cluster

PROTEIN_CODING = CHRONOPACK_INPUT_DIR / "protein_coding_genes.csv"
HUMAN_TFS = CHRONOPACK_INPUT_DIR / "hs_tfs.txt"
LIGAND_RECEPTORS = CHRONOPACK_INPUT_DIR / "interaction_input_CellChatDB.csv"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--live-omnipath", action="store_true", help="Rank TFs from the current OmniPath instead of the snapshot.")
    parser.add_argument("--contexts", nargs="+", default=list(CONTEXTS), choices=list(CONTEXTS))
    args = parser.parse_args()

    for context in args.contexts:
        name, source, is_ligand = CONTEXTS[context]
        out_dir = CHRONOPACK_OUTPUT_DIR / context
        out_dir.mkdir(parents=True, exist_ok=True)
        tfs = cp.tl.grn.rank_regulators(
            source=source, is_ligand=is_ligand, human_tf_path=HUMAN_TFS, ligand_receptor_path=LIGAND_RECEPTORS,
            omnipath_snapshot=None if args.live_omnipath else OMNIPATH_SNAPSHOT,
        )
        targets = cp.tl.grn.score_targets(
            CHRONOPACK_INPUT_DIR / context / "asso_df.1trajectory.nknots4.csv", valid_targets=PROTEIN_CODING
        )
        grn = cp.tl.grn.infer(targets, tfs, valid_targets=PROTEIN_CODING, name=name)

        targets.write(out_dir, prefix="target_scores")
        tfs.to_csv(out_dir / "tf_ranking.csv", index=False)
        grn.write(out_dir, prefix="grn")
        print(f"{name}: {len(grn.edges)} edges, {grn.edges['regulator'].nunique()} regulators, {grn.edges['target'].nunique()} targets -> {out_dir}")
        print(targets.summary.to_string(index=False))
    write_provenance(args.live_omnipath)


def write_provenance(live: bool) -> None:
    repo = Path(cp.__file__).resolve().parents[2]
    try:
        commit = subprocess.run([GIT, "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run([GIT, "-C", str(repo), "status", "--porcelain", "--", "src"], capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    defaults = {k: v.default for k, v in inspect.signature(cp.tl.grn.rank_regulators).parameters.items()
                if v.default is not inspect.Parameter.empty}
    info = {"chronopack_version": getattr(cp, "__version__", None), "chronopack_repo": str(repo), "chronopack_commit": commit,
            "chronopack_src_modified": bool(dirty) if dirty is not None else None,
            "omnipath": "live" if live else str(OMNIPATH_SNAPSHOT.resolve()),
            "rank_regulators_defaults": {k: str(v) for k, v in defaults.items()}}
    (CHRONOPACK_OUTPUT_DIR / "provenance.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
