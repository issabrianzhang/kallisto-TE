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
#       --single-end   treat every file as its own single-end sample, whatever
#                      its name (e.g. replicates named Ctrl_1 / Ctrl_2)
#       --fastqc       also run FastQC on the raw reads (optional; fastp's own
#                      report already covers the usual checks)
#       --fastp-args "ARGS"   extra fastp options, e.g. "--detect_adapter_for_pe"
#       --force        re-run samples whose trimmed files already exist
#       --force-sc     trim even if the reads look like single-cell data
#                      (see NOTE below; almost never what you want)
#   -h, --help         show this help
#
# File names it understands (sample name = everything before the tag):
#   paired-end :  Sample1_R1.fq.gz + Sample1_R2.fq.gz
#                 Sample1_R1_001.fastq.gz + Sample1_R2_001.fastq.gz  (Illumina)
#                 Sample1_1.fq.gz + Sample1_2.fq.gz                  (SRA / ENA)
#   single-end :  Sample1_R1.fq.gz with no R2 file, or Sample1.fq.gz
#
# Pairs are checked: read 1 and read 2 of a real pair carry the same read
# names. Two files that only LOOK like a pair (e.g. single-end replicates
# Ctrl_1.fq.gz + Ctrl_2.fq.gz) are treated as two single-end samples.
#
# Outputs:
#   OUT_DIR/trimmed/<sample>_R1.trimmed.fq.gz  (+ _R2 for paired-end)
#   OUT_DIR/reports/<sample>.fastp.html / .json
#   OUT_DIR/multiqc/multiqc_report.html        (all samples in one page)
#   OUT_DIR/fastp_summary.tsv                  (one line per sample)
#
# NOTE: for droplet single-cell RNA-seq (10x Genomics etc.) do NOT trim.
#   Read 1 holds the cell barcode + UMI; `kb count` needs the raw FASTQs.
#   The script checks this: if read 1 of a paired sample is short (<= 30 bp)
#   while read 2 is at least twice as long, which is the layout of 10x-style
#   barcode + cDNA reads, it stops before trimming anything.
# =============================================================================
set -euo pipefail

usage() { awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"; exit "${1:-0}"; }
die()   { echo "ERROR: $*" >&2; exit 1; }

IN="" OUT="" THREADS=8 R1TAG="" FASTQC=0 FASTP_ARGS="" FORCE=0 FORCE_SC=0 SINGLE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)    IN="$2"; shift 2 ;;
        -o|--outdir)   OUT="$2"; shift 2 ;;
        -t|--threads)  THREADS="$2"; shift 2 ;;
        --r1-tag)      R1TAG="$2"; shift 2 ;;
        --single-end)  SINGLE=1; shift ;;
        --fastqc)      FASTQC=1; shift ;;
        --fastp-args)  FASTP_ARGS="$2"; shift 2 ;;
        --force)       FORCE=1; shift ;;
        --force-sc)    FORCE_SC=1; shift ;;
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
if [[ $SINGLE -eq 1 ]]; then
    R1TAG=""                                      # --single-end: no pairing at all
elif [[ -z "$R1TAG" ]]; then
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

# Check that each pair really is read 1 + read 2 of one library: mates carry
# the same read names (first word of the header, ignoring a /1 or /2 suffix).
read_ids() {
    { case "$1" in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac; } 2>/dev/null \
        | head -n 20 | awk 'NR % 4 == 1 { id = $1; sub(/^@/, "", id); sub(/\/[12]$/, "", id); print id }' || true
}
checked=() not_pairs=()
for s in "${samples[@]}"; do
    IFS=$'\t' read -r name r1 r2 <<<"$s"
    if [[ -n "$r2" && "$(read_ids "$r1")" != "$(read_ids "$r2")" ]]; then
        for f in "$r1" "$r2"; do
            checked+=("$(basename "$f" "$(ext_of "$f")")"$'\t'"$f"$'\t')
        done
        not_pairs+=("$(basename "$r1") + $(basename "$r2")")
    else
        checked+=("$s")
    fi
done
samples=("${checked[@]}")
if [[ ${#not_pairs[@]} -gt 0 ]]; then
    {
        echo "NOTE: these files look like pairs by name, but their read names differ, so"
        echo "      they are NOT read 1 + read 2 of one sample. Each is trimmed as its own"
        echo "      single-end sample (typical for replicates named Ctrl_1 / Ctrl_2):"
        printf '        %s\n' "${not_pairs[@]}"
    } >&2
fi

n_pe=0; n_se=0
for s in "${samples[@]}"; do [[ -n "$(cut -f3 <<<"$s")" ]] && n_pe=$((n_pe+1)) || n_se=$((n_se+1)); done
tag_info="read-1 tag: ${R1TAG:-none}"
[[ $SINGLE -eq 1 ]] && tag_info="--single-end"
echo "Found ${#samples[@]} sample(s) in $IN: $n_pe paired-end, $n_se single-end ($tag_info)" >&2
for s in "${samples[@]}"; do
    IFS=$'\t' read -r name r1 r2 <<<"$s"
    printf '  %-30s %s%s\n' "$name" "$(basename "$r1")" "${r2:+  +  $(basename "$r2")}" >&2
done

# ------------------- refuse droplet single-cell reads ------------------------
# Longest read among the first 2,000 reads of a FASTQ (gzip or plain).
max_len() {
    { case "$1" in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac; } 2>/dev/null \
        | head -n 8000 | awk 'NR % 4 == 2 { if (length($0) > m) m = length($0) } END { print m + 0 }' || true
}
sc_like=()
for s in "${samples[@]}"; do
    IFS=$'\t' read -r name r1 r2 <<<"$s"
    [[ -n "$r2" ]] || continue
    l1=$(max_len "$r1"); l2=$(max_len "$r2")
    if (( l1 > 0 && l1 <= 30 && l2 >= 2 * l1 )); then
        sc_like+=("$name (read 1: ${l1} bp, read 2: ${l2} bp)")
    fi
done
if [[ ${#sc_like[@]} -gt 0 ]]; then
    {
        echo
        echo "These samples look like droplet single-cell RNA-seq (10x Genomics, Drop-seq, ...):"
        printf '  %s\n' "${sc_like[@]}"
        echo "Read 1 is short because it holds the cell barcode + UMI. Trimming would break"
        echo "the barcodes, so single-cell reads must NOT be trimmed: give the raw FASTQ"
        echo "files straight to 'kb count' (see docs/03_qc_trimming.md)."
        echo
    } >&2
    if [[ $FORCE_SC -eq 1 ]]; then
        echo "WARNING: --force-sc given, trimming them anyway" >&2
    else
        die "stopped before trimming. If these really are bulk reads, re-run with --force-sc"
    fi
fi

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
# fastp finds adapters in single-end reads by looking for them; if it found
# none, the reads were only quality-trimmed. Say so, with the fix.
no_adapter=()
for s in "${samples[@]}"; do
    IFS=$'\t' read -r name r1 r2 <<<"$s"
    [[ -z "$r2" ]] || continue
    grep -q "No adapter detected" "$OUT/reports/${name}.fastp.log" 2>/dev/null && no_adapter+=("$name")
done
if [[ ${#no_adapter[@]} -gt 0 ]]; then
    {
        echo "WARNING: fastp found no adapter in these single-end samples, so no adapter was cut:"
        printf '  %s\n' "${no_adapter[@]}"
        echo "  If the library has adapters (most Illumina RNA-seq does), re-run with the adapter given:"
        echo "    --force --fastp-args \"--adapter_sequence AGATCGGAAGAGC\"   (start of the Illumina TruSeq adapter; check your library kit)"
        echo
    } >&2
fi

echo "Trimmed reads : $OUT/trimmed/" >&2
echo "QC report     : $OUT/multiqc/multiqc_report.html  (open in a web browser)" >&2
