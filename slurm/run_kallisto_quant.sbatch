#!/bin/bash -l
# -----------------------------------------------------------------------------
# Example SLURM job: quantify genes + TEs in every bulk RNA-seq sample
# Walk-through: docs/04_bulk_kallisto_quant.md
#
# Edit the lines marked EDIT, then submit from the repository folder:
#   mkdir -p logs
#   sbatch slurm/run_kallisto_quant.sbatch
#
# Resources: our own mouse runs used 64G and 24-30 CPUs. Most of the memory is
# for loading the index, so it hardly depends on the number of reads.
# -----------------------------------------------------------------------------
#SBATCH --job-name=kallisto_quant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --partition=batch                       # EDIT: your cluster's partition

set -euo pipefail

FASTQ_DIR=/path/to/project/01_trimmed/trimmed   # EDIT: trimmed reads from step 2
INDEX=/path/to/references/mm10_TE/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx   # EDIT
OUT_DIR=/path/to/project/02_quant               # EDIT: where the results go

conda activate kallisto_te                      # EDIT if your env has another name

# --strand auto works out the library strandedness from the first sample.
# Single-end reads also use -l 200 -s 20 (fragment length mean / sd) unless you
# change them. Add -b 0 to skip bootstraps if you only use DESeq2 / edgeR.
bash scripts/quant/run_kallisto_quant.sh \
    --input   "$FASTQ_DIR" \
    --index   "$INDEX" \
    --outdir  "$OUT_DIR" \
    --strand  auto \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
