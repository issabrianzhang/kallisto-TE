#!/usr/bin/env bash
# =============================================================================
# prepare_te_library.sh
#
# Merge one or more TE consensus FASTA files (e.g. RepBase *.ref files or a
# Dfam export) into a single, kallisto-compatible TE library.
#
#   1. Concatenate the input FASTA files.
#   2. Clean the headers so every TE name contains only letters, numbers,
#      "_", "-" and ".":
#      - RepBase-style headers (fields separated by TABs) keep all fields:
#           >MMERVK10C<TAB>ERV2<TAB>Mus musculus  ->  >MMERVK10C_ERV2_Mus_musculus
#      - any other header (Dfam, RepeatModeler2, EDTA, ...) keeps only the
#        first word, i.e. the usual FASTA sequence name:
#           >rnd-1_family-1#LTR/Gypsy ( Recon Family Size = 25 )
#                                                 ->  >rnd-1_family-1_LTR-Gypsy
#      Then: spaces/tabs -> "_", "/" -> "-", any other symbol (#, :, |) -> "_".
#   3. Handle duplicated names:
#         - identical name + identical sequence  -> keep the first copy only
#         - identical name + different sequence  -> keep both, rename the
#           later copy to NAME_dup1, NAME_dup2, ...
#   4. Stop with an error if any record ends up with an empty sequence.
#
# Usage:
#   prepare_te_library.sh -o TE_library.fa input1.fa [input2.fa ...]
#
# Example (mouse, RepBase 25.05):
#   prepare_te_library.sh -o mmu_RepBase25.05.fa \
#       RepBase25.05.fasta/mousub.ref RepBase25.05.fasta/rodrep.ref
# =============================================================================
set -euo pipefail

usage() { awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"; exit "${1:-0}"; }

OUT=""
while getopts ":o:h" opt; do
    case "$opt" in
        o) OUT="$OPTARG" ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done
shift $((OPTIND - 1))

[[ -n "$OUT" && $# -ge 1 ]] || usage 1
for f in "$@"; do
    [[ -s "$f" ]] || { echo "ERROR: input file not found or empty: $f" >&2; exit 1; }
done
mkdir -p "$(dirname "$OUT")"

# Read every input (gzip or plain), clean headers, de-duplicate, write output.
for f in "$@"; do
    case "$f" in
        *.gz) gzip -dc "$f" ;;
        *)    cat "$f" ;;
    esac
    echo  # make sure each file ends with a newline before the next one starts
done | awk '
    # Clean header lines only; sequence lines are left untouched.
    { sub(/\r$/, "") }
    /^>/ {
        h = substr($0, 2)
        if (h !~ /\t/) sub(/[ ].*$/, "", h)   # not RepBase-style: keep first word
        gsub(/ /,  "_", h)
        gsub(/\t/, "_", h)
        gsub(/\//, "-", h)
        gsub(/[^A-Za-z0-9_.-]/, "_", h)    # any other symbol (#, :, |, ...) -> "_"
        $0 = ">" h
    }
    NF == 0 { next }                       # drop blank lines

    function flush(   key, i, n, newname) {
        if (name == "") return
        if (name in seq) {
            if (seq[name] == s) {
                printf "  dropped identical duplicate:   %s\n", substr(name, 2) > "/dev/stderr"
                ndrop++
            } else {
                n = ++dup[name]
                newname = name "_dup" n
                printf "  renamed different duplicate:   %s -> %s\n", substr(name, 2), substr(newname, 2) > "/dev/stderr"
                order[++k] = newname; seq[newname] = s; raw[newname] = r
                nren++
            }
        } else {
            order[++k] = name; seq[name] = s; raw[name] = r
        }
    }
    /^>/ { flush(); name = $0; s = ""; r = ""; next }
         { s = s toupper($0); r = r $0 "\n" }
    END {
        flush()
        nempty = 0
        for (i = 1; i <= k; i++) {
            if (length(seq[order[i]]) == 0) {
                printf "ERROR: record has an empty sequence: %s\n", substr(order[i], 2) > "/dev/stderr"
                nempty++
            }
            printf "%s\n%s", order[i], raw[order[i]]
        }
        printf "  %d TE consensus sequences written (%d identical duplicates dropped, %d renamed)\n", \
               k, ndrop + 0, nren + 0 > "/dev/stderr"
        if (nempty > 0) exit 2
    }
' > "$OUT"

echo "TE library written to: $OUT" >&2
