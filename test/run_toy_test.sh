#!/usr/bin/env bash
# =============================================================================
# run_toy_test.sh  --  1-minute check that your installation works
#
# Builds a tiny fake "genome" (2 chromosomes, 2 genes) and a fake RepBase
# archive (4 TEs, one of them duplicated on purpose), runs build_index.sh on
# them and checks that every gene and TE ends up in the index.
#
# Usage (from the repository folder, with the conda env active):
#   bash test/run_toy_test.sh            # writes to test/toy_out/
# =============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$REPO/test/toy_out"
rm -rf "$OUT"; mkdir -p "$OUT/RepBaseToy.fasta"
cd "$OUT"

# Reproducible random DNA of a given length, wrapped at 70 columns.
rnd() { awk -v n="$1" -v seed="$2" 'BEGIN{srand(seed); s="ACGT";
        for(i=1;i<=n;i++){printf "%s", substr(s,int(rand()*4)+1,1); if(i%70==0) printf "\n"}
        if(n%70) printf "\n"}'; }

# --- fake RepBase archive (TAB-separated headers, like the real thing) -------
{ printf '>MMERVK10C\tERVK\tMus musculus\n'; rnd 600 1
  printf '>L1Md_A\tL1\tMus musculus\n';      rnd 900 2; } > RepBaseToy.fasta/mousub.ref
{ printf '>RLTR10C\tERVK\tMus musculus\n';   rnd 400 3
  printf '>MER/ODN\tDNA transposon\tRodentia\n'; rnd 300 4
  printf '>RLTR10C\tERVK\tMus musculus\n';   rnd 400 3; } > RepBaseToy.fasta/rodrep.ref   # duplicate
tar -czf RepBaseToy.fasta.tar.gz RepBaseToy.fasta

# --- fake genome + GTF --------------------------------------------------------
{ echo '>chr1'; rnd 3000 5; echo '>chr2'; rnd 2000 6; } | gzip > genome.fa.gz
{
printf 'chr1\ttoy\tgene\t101\t1100\t.\t+\t.\tgene_id "G1"; gene_name "Gapdh";\n'
printf 'chr1\ttoy\ttranscript\t101\t1100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Gapdh";\n'
printf 'chr1\ttoy\texon\t101\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Gapdh";\n'
printf 'chr1\ttoy\texon\t701\t1100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Gapdh";\n'
printf 'chr2\ttoy\tgene\t201\t900\t.\t-\t.\tgene_id "G2"; gene_name "Actb";\n'
printf 'chr2\ttoy\ttranscript\t201\t900\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "Actb";\n'
printf 'chr2\ttoy\texon\t201\t900\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "Actb";\n'
} | gzip > genes.gtf.gz

# --- run the real pipeline ----------------------------------------------------
bash "$REPO/scripts/index/build_index.sh" \
    --species custom --name toy \
    --repbase RepBaseToy.fasta.tar.gz --te-libs "mousub.ref rodrep.ref" \
    --genome genome.fa.gz --gtf genes.gtf.gz \
    --outdir "$OUT/ref" --threads 2 > "$OUT/build.log" 2>&1 \
    || { echo "FAILED - see $OUT/build.log"; tail -20 "$OUT/build.log"; exit 1; }

# --- expectations -------------------------------------------------------------
n_tx=$(wc -l < "$OUT/ref/index/t2g.txt" | tr -d ' ')
n_te=$(grep -c $'\tTE:' "$OUT/ref/index/t2g.txt" || true)
if [[ "$n_tx" -eq 6 && "$n_te" -eq 4 && -s "$OUT/ref/index/toy_genes_TE.idx" ]]; then
    echo "PASS: index built with 2 genes + 4 TEs (duplicate RLTR10C removed)"
    echo
    cat "$OUT/ref/index/t2g.txt"
else
    echo "FAIL: expected 6 transcripts / 4 TEs, got $n_tx / $n_te - see $OUT/build.log"
    exit 1
fi
