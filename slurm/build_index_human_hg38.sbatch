#!/bin/bash -l
# -----------------------------------------------------------------------------
# Example SLURM job: human GRCh38/hg38 (GENCODE v40) + RepBase TE index
#
# Edit the lines marked EDIT, then submit from the repository folder:
#   sbatch slurm/build_index_human_hg38.sbatch
# -----------------------------------------------------------------------------
#SBATCH --job-name=kallisto_TE_index_hg38
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --partition=batch                 # EDIT: your cluster's partition
##SBATCH --mail-user=you@institution.edu  # EDIT (optional): remove one "#" to enable
##SBATCH --mail-type=END,FAIL

set -euo pipefail

REPBASE=/path/to/RepBase28.04.fasta.tar.gz   # EDIT: RepBase archive you downloaded
OUTDIR=/path/to/references/hg38_GencodeV40_RepBase28.04_TE   # EDIT

conda activate kallisto_te                   # EDIT if your env has another name

bash scripts/index/build_index.sh \
    --species human \
    --repbase "$REPBASE" \
    --name hg38_GencodeV40_RepBase28.04 \
    --outdir "$OUTDIR" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
