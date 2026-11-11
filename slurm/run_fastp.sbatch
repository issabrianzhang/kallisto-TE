#!/bin/bash -l
# -----------------------------------------------------------------------------
# Example SLURM job: trim and QC every bulk RNA-seq sample in one folder
# Walk-through: docs/03_qc_trimming.md
#
# Edit the lines marked EDIT, then submit from the repository folder:
#   mkdir -p logs
#   sbatch slurm/run_fastp.sbatch
#
# fastp is fast: roughly a few minutes per sample of ~30 million read pairs.
# -----------------------------------------------------------------------------
#SBATCH --job-name=fastp_trim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --partition=batch                 # EDIT: your cluster's partition

set -euo pipefail

FASTQ_DIR=/path/to/project/00_fastq       # EDIT: folder with the raw *.fq.gz files
OUT_DIR=/path/to/project/01_trimmed       # EDIT: where trimmed reads + reports go

conda activate kallisto_te                # EDIT if your env has another name

bash scripts/qc/run_fastp.sh \
    --input   "$FASTQ_DIR" \
    --outdir  "$OUT_DIR" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
