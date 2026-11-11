#!/usr/bin/env bash
# =============================================================================
# extract_repbase_species.sh
#
# Pull every RepBase consensus sequence whose species field matches a name,
# from ALL *.ref files in a RepBase folder. Useful when you work on a species
# that has no ready-made preset and you do not know which *.ref file holds it.
#
# RepBase headers look like:   >NAME<TAB>CLASS<TAB>Species name
# The match is on the 3rd (species) field, case-insensitive, as a substring,
# so "Arabidopsis" matches "Arabidopsis thaliana" and "Arabidopsis lyrata".
#
# Usage:
#   extract_repbase_species.sh -d RepBaseXX.XX.fasta -p "Species name" -o out.fa [-x]
#
#   -d  extracted RepBase folder (contains *.ref files)
#   -p  species pattern, e.g. "Arabidopsis thaliana", "Danio rerio", "Rodentia"
#   -o  output FASTA (raw RepBase headers; clean it with prepare_te_library.sh)
#   -x  exact match on the species field instead of substring
#   -l  only list how many matching sequences each *.ref file contains
#
# Examples:
#   extract_repbase_species.sh -d RepBase28.04.fasta -p "Arabidopsis thaliana" -l
#   extract_repbase_species.sh -d RepBase28.04.fasta -p "Arabidopsis thaliana" -o ath_raw.fa
# =============================================================================
set -euo pipefail
usage() { awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"; exit "${1:-0}"; }

DIR="" PATTERN="" OUT="" EXACT=0 LIST=0
while getopts ":d:p:o:xlh" opt; do
    case "$opt" in
        d) DIR="$OPTARG" ;; p) PATTERN="$OPTARG" ;; o) OUT="$OPTARG" ;;
        x) EXACT=1 ;; l) LIST=1 ;; h) usage 0 ;; *) usage 1 ;;
    esac
done
[[ -d "$DIR" && -n "$PATTERN" ]] || usage 1
[[ $LIST -eq 1 || -n "$OUT" ]] || usage 1

shopt -s nullglob
refs=("$DIR"/*.ref "$DIR"/*/*.ref)
[[ ${#refs[@]} -gt 0 ]] || { echo "ERROR: no *.ref files under $DIR" >&2; exit 1; }

select_awk='
    BEGIN { FS = "\t"; p = tolower(pat) }
    /^>/ {
        sp = tolower($3)
        keep = exact ? (sp == p) : (index(sp, p) > 0)
        if (keep) n++
    }
    keep && !list { print }
    END { if (list && n > 0) printf "%6d  %s\n", n, FILENAME }
'

if [[ $LIST -eq 1 ]]; then
    echo "Matching sequences per RepBase file (species ~ \"$PATTERN\"):"
    for f in "${refs[@]}"; do
        awk -v pat="$PATTERN" -v exact="$EXACT" -v list=1 "$select_awk" "$f"
    done
    exit 0
fi

mkdir -p "$(dirname "$OUT")"
for f in "${refs[@]}"; do
    awk -v pat="$PATTERN" -v exact="$EXACT" -v list=0 "$select_awk" "$f"
done > "$OUT"
n=$(grep -c '^>' "$OUT" || true)
echo "$n sequences matching \"$PATTERN\" written to $OUT" >&2
echo "Note: the same TE can appear in more than one *.ref file; prepare_te_library.sh removes identical duplicates." >&2
[[ $n -gt 0 ]] || { echo "ERROR: nothing matched" >&2; exit 1; }
