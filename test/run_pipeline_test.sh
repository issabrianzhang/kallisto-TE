#!/usr/bin/env bash
# =============================================================================
# run_pipeline_test.sh  --  2-minute check that the Snakemake pipeline works
#
# Builds the toy genome + TE library (test/run_toy_test.sh), simulates small
# bulk and 10x single-cell FASTQ files, and runs every pipeline step on them:
# build index -> fastp -> kallisto quant -> kb count.
#
# Usage (from the repository folder, with the conda env active):
#   bash test/run_pipeline_test.sh        # writes to test/pipeline_out/
# =============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$REPO/test/pipeline_out"
cd "$REPO"
command -v snakemake >/dev/null 2>&1 || { echo "snakemake not found: conda activate kallisto_te"; exit 1; }

echo "1/4 toy genome, gene annotation and TE library"
bash test/run_toy_test.sh > /dev/null
rm -rf "$OUT"; mkdir -p "$OUT"

echo "2/4 simulating reads"
python3 test/simulate_toy_reads.py test/toy_out/ref/index/cdna.fa "$OUT/fastq"

cat > "$OUT/config.yaml" <<EOF
outdir: $OUT/results
index:
  species: custom
  name: toy
  dfam: false
  repbase: $REPO/test/toy_out/RepBaseToy.fasta.tar.gz
  repbase_libs: "mousub.ref rodrep.ref"
  genome: $REPO/test/toy_out/genome.fa.gz
  gtf: $REPO/test/toy_out/genes.gtf.gz
bulk:
  fastq_dir: $OUT/fastq/bulk
  bootstraps: 0
singlecell:
  fastq_dir: $OUT/fastq/singlecell
  tech: 10xv3
  memory: 1G
resources:
  build_index: {threads: 2, mem_mb: 4000, runtime: 60}
  trim_bulk:   {threads: 2, mem_mb: 4000, runtime: 60}
  quant_bulk:  {threads: 2, mem_mb: 4000, runtime: 60}
  kb_count:    {threads: 2, mem_mb: 4000, runtime: 60}
EOF

echo "3/4 running the pipeline"
snakemake --cores 2 --configfile "$OUT/config.yaml" --directory "$OUT" \
    --snakefile "$REPO/workflow/Snakefile" > "$OUT/snakemake.log" 2>&1 \
    || { echo "FAILED - see $OUT/snakemake.log and $OUT/results/logs/"; tail -20 "$OUT/snakemake.log"; exit 1; }

echo "4/4 checking the results"
R="$OUT/results"
fail=0
check() { if eval "$2"; then echo "  ok   $1"; else echo "  FAIL $1"; fail=1; fi; }
check "index built"                    "test -s $R/reference/index/toy_genes_TE.idx"
check "4 bulk samples trimmed"         "test \$(tail -n +2 $R/01_trimmed/fastp_summary.tsv | wc -l) -eq 4"
check "4 bulk samples quantified (rf)" "test \$(grep -c $'\trf\t' $R/02_quant/quant_summary.tsv) -eq 4"
check "1 single-cell sample, 30 cells" "awk -F'\t' 'NR==2{exit !(\$5==30)}' $R/03_kb_count/kb_count_summary.tsv"
check "second run has nothing to do"   "snakemake -n --configfile $OUT/config.yaml --directory $OUT --snakefile $REPO/workflow/Snakefile 2>&1 | grep -q 'Nothing to be done'"
if [[ $fail -eq 0 ]]; then
    echo "PASS: the pipeline works. Results are in $R"
else
    echo "FAIL - see $OUT/snakemake.log and $R/logs/"; exit 1
fi
