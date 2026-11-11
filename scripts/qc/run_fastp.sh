#!/usr/bin/env bash
# =============================================================================
# run_fastp.sh  --  trim adapters / low-quality bases from bulk RNA-seq reads
#
# For every sample in a folder of FASTQ files this script:
#   1. finds the read files (paired-end R1 + R2, or single-end)
#   2. runs fastp (adapter + quality trimming, with fastp's default settings)
#   3. writes one HTML + JSON QC report per sample
#   4. combines all reports into one MultiQC report
#   5. prints a summary table (reads kept, % Q30, duplication) per sample
#
# Usage:
#   run_fastp.sh -i FASTQ_DIR -o OUT_DIR [options]
#
# Required:
#   -i, --input        folder with the raw FASTQ files (.fq.gz / .fastq.gz / .fq / .fastq)
#   -o, --outdir       output folder (created if needed)
#
# Options:
#   -t, --threads      threads per fastp run (default: 8; fastp uses at most 16)
#       --r1-tag TAG   text that marks read 1 in file names, e.g. "_R1" or "_1"
#                      (default: auto-detect "_R1_001", "_R1" or "_1")
#       --fastqc       also run FastQC on the raw reads (optional; fastp's own
#                      report already covers the usual checks)
#       --fastp-args "ARGS"   extra fastp options, e.g. "--detect_adapter_for_pe"
#       --force        re-run samples whose trimmed files already exist
#   -h, --help         show this help
#
# File names it understands (sample name = everything before the tag):
#   paired-end :  Sample1_R1.fq.gz + Sample1_R2.fq.gz
#                 Sample1_R1_001.fastq.gz + Sample1_R2_001.fastq.gz  (Illumina)
#                 Sample1_1.fq.gz + Sample1_2.fq.gz                  (SRA / ENA)
#   single-end :  Sample1_R1.fq.gz with no R2 file, or Sample1.fq.gz
#
# Outputs:
#   OUT_DIR/trimmed/<sample>_R1.trimmed.fq.gz  (+ _R2 for paired-end)
#   OUT_DIR/reports/<sample>.fastp.html / .json
#   OUT_DIR/multiqc/multiqc_report.html        (all samples in one page)
#   OUT_DIR/fastp_summary.tsv                  (one line per sample)
#
# NOTE: for droplet single-cell RNA-seq (10x Genomics etc.) do NOT trim.
#   Read 1 holds the cell barcode + UMI; `kb count` needs the raw FASTQs.
# =============================================================================
set -euo pipefail

usage() { awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"; exit "${1:-0}"; }
die()   { echo "ERROR: $*" >&2; exit 1; }

IN="" OUT="" THREADS=8 R1TAG="" FASTQC=0 FASTP_ARGS="" FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)    IN="$2"; shift 2 ;;
        -o|--outdir)   OUT="$2"; shift 2 ;;
        -t|--threads)  THREADS="$2"; shift 2 ;;
        --r1-tag)      R1TAG="$2"; shift 2 ;;
        --fastqc)      FASTQC=1; shift ;;
        --fastp-args)  FASTP_ARGS="$2"; shift 2 ;;
        --force)       FORCE=1; shift ;;
        -h|--help)     usage 0 ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done
[[ -n "$IN" && -n "$OUT" ]] || usage 1
[[ -d "$IN" ]] || die "input folder not found: $IN"
command -v fastp >/dev/null 2>&1 || die "fastp not found. Activate the environment first: conda activate kallisto_te"

mkdir -p "$OUT"/{trimmed,reports}
OUT="$(cd "$OUT" && pwd)"
IN="$(cd "$IN" && pwd)"

# ----------------------------- find the files --------------------------------
EXTS=".fastq.gz .fq.gz .fastq .fq"
shopt -s nullglob
all=()
for e in $EXTS; do for f in "$IN"/*"$e"; do all+=("$f"); done; done
[[ ${#all[@]} -gt 0 ]] || die "no FASTQ files (.fq.gz/.fastq.gz/.fq/.fastq) in $IN"

ext_of() { for e in $EXTS; do [[ "$1" == *"$e" ]] && { echo "$e"; return; }; done; }
mate_tag() { case "$1" in *R1_001) echo "${1%R1_001}R2_001" ;; *R1) echo "${1%R1}R2" ;; *1) echo "${1%1}2" ;; esac; }

# Auto-detect the read-1 tag: the first of these that matches any file.
if [[ -z "$R1TAG" ]]; then
    for t in _R1_001 _R1 _1; do
        for f in "${all[@]}"; do
            e=$(ext_of "$f"); b=$(basename "$f" "$e")
            [[ "$b" == *"$t" ]] && { R1TAG="$t"; break 2; }
        done
    done
fi
R2TAG=""
[[ -n "$R1TAG" ]] && R2TAG=$(mate_tag "$R1TAG")

# Build the sample list: "name<TAB>R1<TAB>R2" (R2 empty for single-end).
samples=()
for f in "${all[@]}"; do
    e=$(ext_of "$f"); b=$(basename "$f" "$e")
    if [[ -n "$R1TAG" && "$b" == *"$R1TAG" ]]; then
        name="${b%"$R1TAG"}"
        r2="$IN/${name}${R2TAG}${e}"
        [[ -f "$r2" ]] || r2=""
        samples+=("$name"$'\t'"$f"$'\t'"$r2")
    elif [[ -n "$R2TAG" && "$b" == *"$R2TAG" && -f "$IN/${b%"$R2TAG"}${R1TAG}${e}" ]]; then
        continue                                  # handled together with its R1
    else
        samples+=("$b"$'\t'"$f"$'\t')             # single-end, no tag
    fi
done
[[ ${#samples[@]} -gt 0 ]] || die "could not find any samples in $IN"

n_pe=0; n_se=0
for s in "${samples[@]}"; do [[ -n "$(cut -f3 <<<"$s")" ]] && n_pe=$((n_pe+1)) || n_se=$((n_se+1)); done
echo "Found ${#samples[@]} sample(s) in $IN: $n_pe paired-end, $n_se single-end (read-1 tag: ${R1TAG:-none})" >&2
for s in "${samples[@]}"; do
    IFS=$'\t' read -r name r1 r2 <<<"$s"
    printf '  %-30s %s%s\n' "$name" "$(basename "$r1")" "${r2:+  +  $(basename "$r2")}" >&2
done

# ----------------------------- optional FastQC -------------------------------
if [[ $FASTQC -eq 1 ]]; then
    command -v fastqc >/dev/null 2>&1 || die "--fastqc given but fastqc is not installed (mamba install -c bioconda fastqc)"
    mkdir -p "$OUT/fastqc_raw"
    echo "Running FastQC on the raw reads" >&2
    fastqc -q -t "$THREADS" -o "$OUT/fastqc_raw" "${all[@]}"
fi

# ----------------------------- fastp per sample ------------------------------
FT=$(( THREADS > 16 ? 16 : THREADS ))
for s in "${samples[@]}"; do
    IFS=$'\t' read -r name r1 r2 <<<"$s"
    o1="$OUT/trimmed/${name}_R1.trimmed.fq.gz"
    o2="$OUT/trimmed/${name}_R2.trimmed.fq.gz"
    [[ -n "$r2" ]] || o1="$OUT/trimmed/${name}.trimmed.fq.gz"
    if [[ $FORCE -eq 0 && -s "$o1" && -s "$OUT/reports/${name}.fastp.json" ]]; then
        echo "  [skip] $name (already trimmed; use --force to redo)" >&2
        continue
    fi
    echo "  [fastp] $name" >&2
    args=(--thread "$FT" -i "$r1" -o "$o1"
          --html "$OUT/reports/${name}.fastp.html"
          --json "$OUT/reports/${name}.fastp.json"
          --report_title "$name")
    [[ -n "$r2" ]] && args+=(-I "$r2" -O "$o2")
    # shellcheck disable=SC2086
    fastp "${args[@]}" $FASTP_ARGS 2> "$OUT/reports/${name}.fastp.log" \
        || die "fastp failed for $name - see $OUT/reports/${name}.fastp.log"
done

# ----------------------------- MultiQC ---------------------------------------
if command -v multiqc >/dev/null 2>&1; then
    echo "Combining reports with MultiQC" >&2
    mq_in=("$OUT/reports")
    [[ $FASTQC -eq 1 ]] && mq_in+=("$OUT/fastqc_raw")
    multiqc -q -f -o "$OUT/multiqc" "${mq_in[@]}" > "$OUT/multiqc.log" 2>&1 \
        || echo "WARNING: MultiQC failed - see $OUT/multiqc.log" >&2
else
    echo "NOTE: multiqc not installed; skipping the combined report" >&2
fi

# ----------------------------- summary table ---------------------------------
python3 - "$OUT/reports" "$OUT/fastp_summary.tsv" <<'EOF'
import json, sys, pathlib
rep, out = pathlib.Path(sys.argv[1]), sys.argv[2]
rows = ["sample\tlayout\treads_before\treads_after\tpercent_kept\tq30_after\tduplication\tadapter_trimmed_reads"]
for j in sorted(rep.glob("*.fastp.json")):
    d = json.load(open(j))
    s = d["summary"]
    before, after = s["before_filtering"]["total_reads"], s["after_filtering"]["total_reads"]
    layout = "paired" if "read2_before_filtering" in d else "single"
    rows.append("\t".join(map(str, [
        j.name[:-len(".fastp.json")], layout, before, after,
        f"{100 * after / before:.1f}" if before else "NA",
        f"{100 * s['after_filtering']['q30_rate']:.1f}",
        f"{100 * d.get('duplication', {}).get('rate', 0):.1f}",
        d.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)])))
open(out, "w").write("\n".join(rows) + "\n")
EOF

echo >&2
echo "Summary (also in $OUT/fastp_summary.tsv):" >&2
column -t -s $'\t' "$OUT/fastp_summary.tsv" >&2 || cat "$OUT/fastp_summary.tsv" >&2
echo >&2
echo "Trimmed reads : $OUT/trimmed/" >&2
echo "QC report     : $OUT/multiqc/multiqc_report.html  (open in a web browser)" >&2
