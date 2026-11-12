#!/usr/bin/env bash
# =============================================================================
# run_kallisto_quant.sh  --  count genes + TEs in bulk RNA-seq with kallisto quant
#
# For every sample in a folder of FASTQ files this script:
#   1. finds the read files (paired-end R1 + R2, or single-end; several lanes
#      of one sample are combined)
#   2. works out the library strandedness from the first sample (or uses --strand)
#   3. runs kallisto quant against the genes + TEs index from step 1
#   4. writes kallisto's results per sample, plus a summary (reads used and
#      % pseudoaligned per sample). Load the results into R with tximport
#      (docs/04_bulk_kallisto_quant.md, section 7).
#
# Usage:
#   run_kallisto_quant.sh -i FASTQ_DIR --index INDEX.idx -o OUT_DIR [options]
#
# Required:
#   -i, --input        folder with the reads: the trimmed/ folder from step 2,
#                      raw FASTQs, or one sub-folder per sample
#       --index        the .idx file from step 1, e.g. ref/index/mm10_genes_TE.idx
#   -o, --outdir       output folder (created if needed)
#
# Options:
#       --strand S     auto (default), unstranded, rf or fr  (see NOTE 1)
#   -l, --fragment-length N   single-end only: mean fragment length (default 200)
#   -s, --fragment-sd N       single-end only: its standard deviation (default 20)
#   -b, --bootstraps N        bootstrap rounds (default 30; see NOTE 2)
#   -t, --threads N           threads (default 8)
#       --r1-tag TAG   text that marks read 1 in file names (default: auto)
#       --single-end   treat every file as single-end, whatever its name
#       --force        re-run samples that already have results
#   -h, --help         show this help
#
# How samples are named:
#   one folder of files : Ctrl1_R1.trimmed.fq.gz + Ctrl1_R2.trimmed.fq.gz -> "Ctrl1"
#                         lanes Ctrl1_S1_L001_R1_001 + _L002_ ...          -> "Ctrl1"
#   one folder / sample : raw/Ctrl1/<any FASTQ files>                       -> "Ctrl1"
#
# Outputs:
#   OUT_DIR/<sample>/abundance.tsv   per-target counts + TPM (genes' transcripts and TEs)
#   OUT_DIR/<sample>/abundance.h5    the same + bootstraps (for sleuth)
#   OUT_DIR/<sample>/run_info.json   reads processed / pseudoaligned
#   OUT_DIR/quant_summary.tsv            one line per sample
#   OUT_DIR/samples_used.tsv             which files went into which sample
#
# NOTE 1 - strandedness. Most modern kits (Illumina TruSeq Stranded, NEB Ultra II
#   Directional, any "dUTP" protocol) are reverse-stranded = "rf". With
#   --strand auto the script quantifies 200,000 reads of the first sample both
#   ways and picks the one that fits; it stops and asks you if the answer is
#   unclear. Strand matters for TEs: many are transcribed antisense.
# NOTE 2 - bootstraps are only used by sleuth / fishpond. For DESeq2 / edgeR
#   they are not needed, and -b 0 makes the run faster and smaller.
# =============================================================================
set -euo pipefail

usage() { awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"; exit "${1:-0}"; }
die()   { echo "ERROR: $*" >&2; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IN="" INDEX="" OUT="" STRAND=auto FRAG_L=200 FRAG_S=20 BOOT=30 THREADS=8
R1TAG="" SINGLE=0 FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)            IN="$2"; shift 2 ;;
        --index)               INDEX="$2"; shift 2 ;;
        -o|--outdir)           OUT="$2"; shift 2 ;;
        --strand)              STRAND="$2"; shift 2 ;;
        -l|--fragment-length)  FRAG_L="$2"; shift 2 ;;
        -s|--fragment-sd)      FRAG_S="$2"; shift 2 ;;
        -b|--bootstraps)       BOOT="$2"; shift 2 ;;
        -t|--threads)          THREADS="$2"; shift 2 ;;
        --r1-tag)              R1TAG="$2"; shift 2 ;;
        --single-end)          SINGLE=1; shift ;;
        --force)               FORCE=1; shift ;;
        -h|--help)             usage 0 ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done
[[ -n "$IN" && -n "$INDEX" && -n "$OUT" ]] || usage 1
[[ -d "$IN" ]]    || die "input folder not found: $IN"
[[ -s "$INDEX" ]] || die "index file not found: $INDEX (it is the .idx file in the index/ folder from step 1)"
case "$STRAND" in auto|unstranded|rf|fr) ;; *) die "--strand must be auto, unstranded, rf or fr" ;; esac
command -v python3 >/dev/null 2>&1 || die "python3 not found. Activate the environment first: conda activate kallisto_te"

# kallisto: use the copy bundled with kb-python (the one that built the index).
KALLISTO=""
if command -v kb >/dev/null 2>&1; then
    KALLISTO=$(kb info 2>/dev/null | sed -n 's/^kallisto: [^ ]* (\(.*\))$/\1/p' | head -1 || true)
fi
[[ -n "$KALLISTO" && -x "$KALLISTO" ]] || KALLISTO=$(command -v kallisto || true)
[[ -n "$KALLISTO" ]] || die "kallisto not found. Activate the environment first: conda activate kallisto_te"

mkdir -p "$OUT/logs"
OUT="$(cd "$OUT" && pwd)"
IN="$(cd "$IN" && pwd)"
INDEX="$(cd "$(dirname "$INDEX")" && pwd)/$(basename "$INDEX")"

# shellcheck source=../lib/fastq_pairs.sh
source "$HERE/../lib/fastq_pairs.sh"

# ----------------------------- find the samples ------------------------------
collect_samples "$IN" "$SINGLE" "$R1TAG"
[[ ${#S_NAME[@]} -gt 0 ]] || die "no FASTQ files (.fq.gz/.fastq.gz/.fq/.fastq) in $IN or its sub-folders"

if [[ ${#NOT_PAIRS[@]} -gt 0 ]]; then
    {
        echo "NOTE: these files look like pairs by name, but their read names differ, so"
        echo "      they are NOT read 1 + read 2 of one sample. Each is used as its own"
        echo "      single-end sample (typical for replicates named Ctrl_1 / Ctrl_2):"
        printf '        %s\n' "${NOT_PAIRS[@]}"
    } >&2
fi

# Refuse droplet single-cell reads: they need kb count, not kallisto quant.
for ((i = 0; i < ${#S_NAME[@]}; i++)); do
    [[ "${S_LAYOUT[i]}" == paired ]] || continue
    IFS=$'\t' read -r r1 r2 _ <<<"${S_FILES[i]}"
    l1=$(max_len "$r1"); l2=$(max_len "$r2")
    if (( l1 > 0 && l1 <= 30 && l2 >= 2 * l1 )); then
        die "sample ${S_NAME[i]} looks like droplet single-cell data (read 1: ${l1} bp barcode/UMI, read 2: ${l2} bp).
       Use scripts/quant/run_kb_count.sh for single-cell reads (docs/05_singlecell_kb_count.md)."
    fi
done

if [[ ${#SKIPPED[@]} -gt 0 ]]; then
    echo "NOTE: ignoring ${#SKIPPED[@]} Illumina index-read file(s) (*_I1_001 / *_I2_001): they hold sample barcodes, not RNA" >&2
fi
n_pe=0; n_se=0
for l in "${S_LAYOUT[@]}"; do [[ "$l" == paired ]] && n_pe=$((n_pe+1)) || n_se=$((n_se+1)); done
echo "Found ${#S_NAME[@]} sample(s): $n_pe paired-end, $n_se single-end" >&2
for ((i = 0; i < ${#S_NAME[@]}; i++)); do
    files=""
    IFS=$'\t' read -r -a fa <<<"${S_FILES[i]}"
    for f in "${fa[@]}"; do files+="${files:+ }$(basename "$f")"; done
    printf '  %-28s %-7s %s\n' "${S_NAME[i]}" "${S_LAYOUT[i]}" "$files" >&2
done

# ----------------------------- strandedness ----------------------------------
# Quantify a subset of the first sample as "fr" and as "rf". A stranded library
# keeps nearly all reads in one orientation; an unstranded one splits ~50/50.
run_kallisto() {   # out-dir layout strand-flag files...
    local o="$1" layout="$2" flag="$3"; shift 3
    local a=(quant -i "$INDEX" -o "$o" -t "$THREADS")
    [[ -n "$flag" ]] && a+=("$flag")
    [[ "$layout" == single ]] && a+=(--single -l "$FRAG_L" -s "$FRAG_S")
    "$KALLISTO" "${a[@]}" "$@"
}
n_aligned() { python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["n_pseudoaligned"])' "$1/run_info.json"; }

if [[ "$STRAND" == auto ]]; then
    chk="$OUT/strand_check"; rm -rf "$chk"; mkdir -p "$chk"
    IFS=$'\t' read -r -a fa <<<"${S_FILES[0]}"
    sub=("$chk/sub_1.fq")
    cat_fastq "${fa[0]}" | head -n 800000 > "$chk/sub_1.fq" || true
    if [[ "${S_LAYOUT[0]}" == paired ]]; then
        cat_fastq "${fa[1]}" | head -n 800000 > "$chk/sub_2.fq" || true
        sub+=("$chk/sub_2.fq")
    fi
    echo "Checking strandedness on the first 200,000 reads of ${S_NAME[0]} ..." >&2
    for m in fr rf; do
        # kallisto exits with an error when no read fits, which is expected
        # for the wrong orientation of a stranded library: only a missing
        # run_info.json means it really failed.
        run_kallisto "$chk/$m" "${S_LAYOUT[0]}" "--$m-stranded" "${sub[@]}" > "$chk/$m.log" 2>&1 || true
        [[ -s "$chk/$m/run_info.json" ]] || die "kallisto failed during the strand check - see $chk/$m.log"
    done
    fr=$(n_aligned "$chk/fr"); rf=$(n_aligned "$chk/rf")
    read -r STRAND verdict < <(python3 - "$fr" "$rf" <<'EOF'
import sys
fr, rf = int(sys.argv[1]), int(sys.argv[2])
tot = fr + rf
if tot < 1000:
    print("unknown", f"only {tot} reads pseudoaligned - is this the right index for these reads?")
else:
    p = rf / tot
    pct = f"{100 * p:.0f}% of strand-matched reads fit 'rf', {100 - 100 * p:.0f}% fit 'fr'"
    if p >= 0.8:   print("rf", f"reverse-stranded ({pct})")
    elif p <= 0.2: print("fr", f"forward-stranded ({pct})")
    elif 0.35 <= p <= 0.65: print("unstranded", f"unstranded ({pct})")
    else:          print("unknown", f"unclear ({pct})")
EOF
)
    {
        echo "  fr-stranded: $fr reads pseudoaligned"
        echo "  rf-stranded: $rf reads pseudoaligned"
        echo "  -> $verdict"
    } | tee "$chk/result.txt" >&2
    [[ "$STRAND" != unknown ]] || die "could not tell the strandedness. Check your library kit and re-run with
       --strand rf (most stranded kits), --strand fr, or --strand unstranded."
    rm -f "$chk"/sub_*.fq
else
    echo "Strandedness: $STRAND (given with --strand)" >&2
fi
case "$STRAND" in rf) SFLAG=--rf-stranded ;; fr) SFLAG=--fr-stranded ;; *) SFLAG="" ;; esac

# ----------------------------- kallisto quant --------------------------------
printf 'sample\tlayout\tstrand\tfiles\n' > "$OUT/samples_used.tsv"
for ((i = 0; i < ${#S_NAME[@]}; i++)); do
    name="${S_NAME[i]}" layout="${S_LAYOUT[i]}"
    IFS=$'\t' read -r -a fa <<<"${S_FILES[i]}"
    printf '%s\t%s\t%s\t%s\n' "$name" "$layout" "$STRAND" "$(IFS=,; echo "${fa[*]}")" >> "$OUT/samples_used.tsv"
    if [[ $FORCE -eq 0 && -s "$OUT/$name/abundance.tsv" && -s "$OUT/$name/run_info.json" ]]; then
        echo "  [skip] $name (already quantified; use --force to redo)" >&2
        continue
    fi
    echo "  [kallisto quant] $name ($layout)" >&2
    rm -rf "${OUT:?}/$name"
    # kallisto also exits with an error when no read pseudoaligns; that sample
    # is reported at 0% below. Only a missing run_info.json is a real failure.
    run_kallisto "$OUT/$name" "$layout" "$SFLAG" -b "$BOOT" "${fa[@]}" > "$OUT/logs/$name.log" 2>&1 || true
    [[ -s "$OUT/$name/run_info.json" ]] || die "kallisto quant failed for $name - see $OUT/logs/$name.log"
done

# ----------------------------- summary ---------------------------------------
# One line per sample, from kallisto's run_info.json.
python3 - "$OUT" <<'EOF'
import csv, json, sys
out = sys.argv[1]
with open(f"{out}/samples_used.tsv") as fh:
    samples = list(csv.DictReader(fh, delimiter="\t"))
with open(f"{out}/quant_summary.tsv", "w") as fh:
    fh.write("sample\tlayout\tstrand\treads_processed\treads_pseudoaligned\tpercent_pseudoaligned\n")
    for s in samples:
        info = json.load(open(f"{out}/{s['sample']}/run_info.json"))
        fh.write(f"{s['sample']}\t{s['layout']}\t{s['strand']}\t{info['n_processed']}"
                 f"\t{info['n_pseudoaligned']}\t{info['p_pseudoaligned']:.1f}\n")
EOF

if command -v multiqc >/dev/null 2>&1; then
    multiqc -q -f -o "$OUT/multiqc" "$OUT/logs" > "$OUT/logs/multiqc.log" 2>&1 \
        || echo "WARNING: MultiQC failed - see $OUT/logs/multiqc.log" >&2
fi

echo >&2
echo "Summary (also in $OUT/quant_summary.tsv):" >&2
column -t -s $'\t' "$OUT/quant_summary.tsv" >&2 || cat "$OUT/quant_summary.tsv" >&2
echo >&2
low=$(awk -F'\t' 'NR > 1 && $6 + 0 < 30 { print "  " $1 " (" $6 "%)" }' "$OUT/quant_summary.tsv")
if [[ -n "$low" ]]; then
    {
        echo "WARNING: under 30% of reads pseudoaligned in:"
        echo "$low"
        echo "  Common causes: index from another species, heavy rRNA / genomic DNA"
        echo "  contamination, or total RNA with many intronic reads (see docs/04_bulk_kallisto_quant.md)."
        echo
    } >&2
fi
echo "Per-sample results : $OUT/<sample>/abundance.tsv" >&2
echo "Next               : load them into R with tximport + DESeq2 (docs/04_bulk_kallisto_quant.md, section 7)" >&2
