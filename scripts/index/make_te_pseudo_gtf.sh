#!/usr/bin/env bash
# =============================================================================
# make_te_pseudo_gtf.sh
#
# Build a "pseudo-GTF" for a TE consensus library. Each TE consensus sequence
# is treated as its own small chromosome that carries exactly one gene
# spanning the whole sequence (position 1 to its length):
#
#   MMERVK10C_ERVK_Mus_musculus  RepBase  gene  1  7323  .  +  .  gene_id "MMERVK10C_ERVK_Mus_musculus"; gene_type "TE"; gene_name "TE:MMERVK10C_ERVK_Mus_musculus";
#
# When this GTF is appended to the genome annotation, `kb ref` turns each TE
# gene into a single transcript/exon, so every TE consensus ends up in the
# kallisto index and in t2g.txt next to the normal genes. The "TE:" prefix
# in gene_name makes TEs easy to pick out later (e.g. grep "TE:").
#
# Usage:
#   make_te_pseudo_gtf.sh TE_library.fa TE_library.pseudo.gtf
# =============================================================================
set -euo pipefail

if [[ $# -ne 2 ]]; then
    awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"
    exit 1
fi
TE_FA="$1"
OUT="$2"
[[ -s "$TE_FA" ]] || { echo "ERROR: TE FASTA not found: $TE_FA" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

# Sequence lengths are computed directly (same numbers as `samtools faidx`).
awk '
    function emit() {
        if (name == "") return
        printf "%s\tRepBase\tgene\t1\t%d\t.\t+\t.\tgene_id \"%s\"; gene_type \"TE\"; gene_name \"TE:%s\";\n", \
               name, len, name, name
        n++
    }
    /^>/ { emit(); name = substr($1, 2); len = 0; next }
         { gsub(/[ \t\r]/, ""); len += length($0) }
    END  { emit(); printf "  %d TE entries written to the pseudo-GTF\n", n > "/dev/stderr" }
' "$TE_FA" > "$OUT"

echo "Pseudo-GTF written to: $OUT" >&2
