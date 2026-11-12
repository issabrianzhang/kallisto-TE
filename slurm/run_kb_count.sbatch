#!/bin/bash -l
# -----------------------------------------------------------------------------
# Example SLURM job: count genes + TEs per cell in droplet single-cell RNA-seq
# Walk-through: docs/05_singlecell_kb_count.md
#
# Edit the lines marked EDIT, then submit from the repository folder:
#   mkdir -p logs
#   sbatch slurm/run_kb_count.sbatch
#
# Use the RAW FASTQ files: single-cell reads are never trimmed.
# Resources: our own 10x v4 mouse runs used 48 CPUs, --mem=196G and
# --memory 120G. The index itself needs under 2 GB; most of the memory goes to
# sorting barcodes (--memory), so less usually works too. If the job runs out
# of memory, raise --mem and keep --memory well below it.
# -----------------------------------------------------------------------------
#SBATCH --job-name=kb_count
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --partition=batch                       # EDIT: your cluster's partition

set -euo pipefail

FASTQ_DIR=/path/to/project/00_fastq_10x          # EDIT: raw 10x FASTQs (R1/R2 of every lane)
INDEX=/path/to/references/mm10_TE/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx   # EDIT
TECH=10xv3                                       # EDIT: 10xv2, 10xv3, 10xv4, ... (kb --list)
OUT_DIR=/path/to/project/03_kb_count             # EDIT: where the results go

conda activate kallisto_te                       # EDIT if your env has another name

bash scripts/quant/run_kb_count.sh \
    --input   "$FASTQ_DIR" \
    --index   "$INDEX" \
    --tech    "$TECH" \
    --outdir  "$OUT_DIR" \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --memory  64G
