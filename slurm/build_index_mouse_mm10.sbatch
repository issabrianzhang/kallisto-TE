#!/bin/bash -l
# -----------------------------------------------------------------------------
# Example SLURM job: mouse GRCm38/mm10 (GENCODE M25) + RepBase TE index
#
# Edit the lines marked EDIT, then submit from the repository folder:
#   sbatch slurm/build_index_mouse_mm10.sbatch
#
# Resources: kb ref on a mammalian genome needs roughly 32-64 GB RAM and a
# few hours. Partition names differ between clusters; ask your HPC team.
# -----------------------------------------------------------------------------
#SBATCH --job-name=kallisto_TE_index_mm10
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

REPBASE=/path/to/RepBase25.05.fasta.tar.gz   # EDIT: RepBase archive you downloaded
OUTDIR=/path/to/references/mm10_GencodeM25_RepBase25.05_TE   # EDIT

conda activate kallisto_te                   # EDIT if your env has another name
# module load ...                            # EDIT: if your cluster needs modules

bash scripts/index/build_index.sh \
    --species mouse \
    --repbase "$REPBASE" \
    --name mm10_GencodeM25_RepBase25.05 \
    --outdir "$OUTDIR" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
