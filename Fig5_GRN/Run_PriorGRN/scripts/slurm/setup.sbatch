#!/usr/bin/env bash
# One-time setup, all inside Run_PriorGRN/ (no system or shared-env changes). Submit from Run_PriorGRN/:
#   sbatch scripts/slurm/setup.sbatch
# - resources/containers/pyscenic_0.12.1.sif      docker://aertslab/pyscenic:0.12.1 (pySCENIC + arboreto)
# - resources/containers/celloracle_0.18.0.sif    docker://kenjikamimoto126/celloracle_ubuntu:0.18.0
# - resources/cistarget/                          hg38 refseq_r80 mc_v10_clust gene-based rankings (500 bp up /
#                                                 100 bp down and 10 kb up / 10 kb down), motif2TF table, TF list
# - .conda_env/                                   project-local Python 3.10 env with inferelator 0.6.4
#SBATCH --job-name=priorgrn_setup
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=output/logs/%x_%j.out
#SBATCH --error=output/logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the Run_PriorGRN directory}"
[[ -f run_prior_grn.sh ]] || { echo "Submit from the Run_PriorGRN directory" >&2; exit 2; }
R=$PWD
mkdir -p resources/containers resources/cistarget

curl -sI -m 30 https://pypi.org >/dev/null || { echo "no internet on $(hostname)" >&2; exit 3; }

# cisTarget databases and annotations (pySCENIC); sha1 checked against the server's sums
B=https://resources.aertslab.org/cistarget
cd resources/cistarget
for f in hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
         hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather; do
  [[ -s $f ]] || curl -sS -L --retry 5 -o $f "$B/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/$f"
  [[ -s $f.sha1sum.txt ]] || curl -sS -L --retry 5 -o $f.sha1sum.txt "$B/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/$f.sha1sum.txt"
  sha1sum -c $f.sha1sum.txt
done
[[ -s motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl ]] || curl -sS -L --retry 5 -o motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl "$B/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
[[ -s hs_hgnc_tfs.txt ]] || curl -sS -L --retry 5 -o hs_hgnc_tfs.txt https://raw.githubusercontent.com/aertslab/pySCENIC/0.12.1/resources/hs_hgnc_tfs.txt
wc -l hs_hgnc_tfs.txt motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
cd "$R"

# Containers: built on the node's local /scratch (unpacking layers on the NFS home took > 9 h for one 0.3 GB image,
# job 942924), then copied into the project; cache and temp removed afterwards
L=/scratch/$USER/priorgrn_$SLURM_JOB_ID
export SINGULARITY_CACHEDIR=$L/cache SINGULARITY_TMPDIR=$L/tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"
trap 'rm -rf "$L"' EXIT
pull() {  # pull <sif name> <docker ref>
  [[ -s resources/containers/$1 ]] && return
  singularity pull "$L/$1" "$2"
  cp "$L/$1" "resources/containers/$1.part" && mv "resources/containers/$1.part" "resources/containers/$1"
}
pull pyscenic_0.12.1.sif docker://aertslab/pyscenic:0.12.1
pull celloracle_0.18.0.sif docker://kenjikamimoto126/celloracle_ubuntu:0.18.0
singularity exec --cleanenv resources/containers/pyscenic_0.12.1.sif bash -c 'pyscenic --help | head -3; command -v arboreto_with_multiprocessing.py'
singularity exec --cleanenv resources/containers/celloracle_0.18.0.sif python -c 'import celloracle, scanpy; print("celloracle", celloracle.__version__, "scanpy", scanpy.__version__)'

# Inferelator env (project-local prefix)
CONDA=/new-stg/home/alex/miniconda3/bin/conda
[[ -x .conda_env/bin/python ]] || $CONDA create -y -q -p "$R/.conda_env" python=3.10 pip
.conda_env/bin/pip install -q "inferelator==0.6.4" "anndata>=0.9" pyarrow
.conda_env/bin/python -c 'import inferelator, anndata, pyarrow; print("inferelator OK", anndata.__version__, pyarrow.__version__)'
echo PRIORGRN_SETUP_OK
