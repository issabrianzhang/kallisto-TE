#!/usr/bin/env bash
# =============================================================================
# check_index.sh
#
# Sanity checks after `kb ref`:
#   - every TE in the TE library is present in t2g.txt (and in the cDNA FASTA,
#     if one is given)
#   - reports how many genes / transcripts / TEs the index contains
#
# Exits with status 1 if any TE is missing, so a broken index never goes
# unnoticed.
#
# Usage:
#   check_index.sh TE_library.fa t2g.txt [cdna.fa]
# =============================================================================
set -euo pipefail

if [[ $# -lt 2 ]]; then
    awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"
    exit 1
fi
TE_FA="$1"; T2G="$2"; CDNA="${3:-}"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

grep '^>' "$TE_FA" | sed 's/^>//; s/[[:space:]].*//' | sort -u > "$tmp/te_names"
cut -f1 "$T2G" | sort -u > "$tmp/t2g_tx"

n_te=$(wc -l < "$tmp/te_names" | tr -d ' ')
n_tx=$(wc -l < "$tmp/t2g_tx" | tr -d ' ')
n_genes=$(cut -f2 "$T2G" | sort -u | wc -l | tr -d ' ')
n_te_t2g=$(grep -c $'\tTE:' "$T2G" || true)

echo "t2g.txt: ${n_tx} transcripts, ${n_genes} genes (${n_te_t2g} TE entries)"
echo "TE library: ${n_te} TE consensus sequences"

comm -23 "$tmp/te_names" "$tmp/t2g_tx" > "$tmp/missing"
status=0
if [[ -s "$tmp/missing" ]]; then
    echo "ERROR: $(wc -l < "$tmp/missing" | tr -d ' ') TE(s) are missing from t2g.txt:" >&2
    head -20 "$tmp/missing" >&2
    status=1
else
    echo "OK: all ${n_te} TEs are present in t2g.txt"
fi

if [[ -n "$CDNA" ]]; then
    grep '^>' "$CDNA" | sed 's/^>//; s/[[:space:]].*//' | sort -u > "$tmp/cdna_tx"
    comm -23 "$tmp/te_names" "$tmp/cdna_tx" > "$tmp/missing_cdna"
    if [[ -s "$tmp/missing_cdna" ]]; then
        echo "ERROR: $(wc -l < "$tmp/missing_cdna" | tr -d ' ') TE(s) are missing from the cDNA FASTA" >&2
        status=1
    else
        echo "OK: all ${n_te} TEs are present in $(basename "$CDNA")"
    fi
fi
exit $status
