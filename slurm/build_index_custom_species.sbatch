#!/bin/bash -l
# -----------------------------------------------------------------------------
# Example SLURM job: ANY species (here Arabidopsis thaliana TAIR10)
# Walk-through: docs/02_custom_species.md
#
# Edit the lines marked EDIT, then submit from the repository folder:
#   sbatch slurm/build_index_custom_species.sbatch
#
# Small genomes (Arabidopsis ~120 Mb) need far less: ~8-16 GB RAM, < 1 hour.
# -----------------------------------------------------------------------------
#SBATCH --job-name=kallisto_TE_index_custom
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --partition=batch                 # EDIT: your cluster's partition

set -euo pipefail

NAME=TAIR10_Ensembl63                                   # EDIT: label for output files
GENOME=/path/to/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz   # EDIT
GTF=/path/to/Arabidopsis_thaliana.TAIR10.63.gtf.gz              # EDIT
TE_FASTA=/path/to/Arabidopsis_thaliana_TE_consensus.fa         # EDIT (RepBase extract, Dfam, EDTA, ...)
OUTDIR=/path/to/references/${NAME}_TE                         # EDIT

conda activate kallisto_te

bash scripts/index/build_index.sh \
    --species custom \
    --name "$NAME" \
    --genome "$GENOME" \
    --gtf "$GTF" \
    --te-fasta "$TE_FASTA" \
    --outdir "$OUTDIR" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
