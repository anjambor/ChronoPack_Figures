"""Small helpers shared by the figure notebooks (one copy per Fig* folder).

- ``run(cmd, ...)`` runs one pipeline step as a subprocess from the figure folder, streaming its
  output into the notebook and raising if it fails.
- ``needs_run(outputs, rerun)`` decides whether a step runs: always when ``rerun`` is True,
  otherwise only if one of its outputs is missing.
- ``show_*`` display figures, tables and text reports inline.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd
from IPython.display import HTML, Image, Markdown, display

FIG_DIR = Path(__file__).resolve().parents[1]


def needs_run(outputs: Iterable[str | Path], rerun: bool) -> bool:
    missing = [str(p) for p in outputs if not (FIG_DIR / p).exists()]
    if rerun:
        return True
    if missing:
        print("Outputs missing, running this step:", *missing, sep="\n  ")
        return True
    print("Using existing outputs (set RERUN = True to regenerate).")
    return False


def run(cmd: Sequence[str | Path], env: Mapping[str, str] | None = None, prepend_path: str | Path | None = None) -> None:
    """Run ``cmd`` from the figure folder; ``env`` adds/overrides variables, ``prepend_path`` is put on PATH."""
    full_env = os.environ.copy()
    # Jupyter exports MPLBACKEND (inline backend) to the kernel; drop it so steps render as on the command line.
    full_env.pop("MPLBACKEND", None)
    if prepend_path is not None:
        full_env["PATH"] = f"{prepend_path}{os.pathsep}{full_env.get('PATH', '')}"
    full_env.update({k: str(v) for k, v in (env or {}).items()})
    cmd = [str(c) for c in cmd]
    print("$", " ".join(cmd), flush=True)
    proc = subprocess.Popen(cmd, cwd=FIG_DIR, env=full_env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(line)
    if proc.wait() != 0:
        raise RuntimeError(f"Step failed with exit code {proc.returncode}: {' '.join(cmd)}")


def show_image(path: str | Path, width: int = 1000, caption: str | None = None) -> None:
    path = FIG_DIR / path
    if caption:
        display(Markdown(f"**{caption}** — `{path.relative_to(FIG_DIR)}`"))
    display(Image(filename=str(path), width=width))


def show_table(path: str | Path, n: int | None = 15, **read_kw) -> pd.DataFrame:
    df = pd.read_csv(FIG_DIR / path, **read_kw)
    with pd.option_context("display.max_columns", 40, "display.width", 200):
        display(df if n is None else df.head(n))
    return df


def show_text(path: str | Path) -> None:
    text = (FIG_DIR / path).read_text()
    display(HTML(f"<pre style='font-size:11px'>{text}</pre>"))
