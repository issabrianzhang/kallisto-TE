#!/bin/bash -l
# =============================================================================
# OPTIONAL: RepeatMasker annotation of the genome with the SAME TE library
#
# Not needed to build or use the kallisto-TE index. kallisto-TE quantifies TE
# expression at the consensus (subfamily) level. Run this only if you also
# want locus-level TE coordinates (GFF / .out) that use exactly the same TE
# names as the index, e.g. to see which genomic copies of an up-regulated
# subfamily sit near differentially expressed genes, or for locus-level
# counting with featureCounts/TElocal.
#
# Expect several days for a mammalian genome with the -s (slow) setting.
#
# Usage (edit the EDIT lines first):
#   sbatch scripts/optional/repeatmasker_genome.sbatch
#
# Inputs:
#   TE_LIB  the cleaned library made by build_index.sh
#           (OUTDIR/te/<name>_TE_library.fa)
#   GENOME  the uncompressed genome FASTA
# =============================================================================
#SBATCH --job-name=RepeatMasker
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --partition=batch                 # EDIT

set -euo pipefail

TE_LIB=/path/to/references/mm10_TE/te/mm10_TE_library.fa       # EDIT
GENOME=/path/to/GRCm38.primary_assembly.genome.fa               # EDIT (uncompressed)
OUTDIR=/path/to/references/mm10_TE/repeatmasker                 # EDIT

module load RepeatMasker   # EDIT: or `conda install -c bioconda repeatmasker`
mkdir -p "$OUTDIR"

# Each rmblast job uses ~4 threads, so -pa is about CPUs / 4.
PA=$(( ${SLURM_CPUS_PER_TASK:-16} / 4 ))
[[ $PA -ge 1 ]] || PA=1

# -lib     custom TE consensus library (same names as the kallisto index)
# -s       slow / most sensitive search
# -a       write alignments (.align), needed for divergence (age) estimates
# -gff     also write a GFF of the TE loci
# -xm      write the cross_match style .out.xm file
# -engine  rmblast search engine
RepeatMasker \
    -lib "$TE_LIB" \
    -pa "$PA" -s -a -gff -xm \
    -engine rmblast \
    -dir "$OUTDIR" \
    "$GENOME"
