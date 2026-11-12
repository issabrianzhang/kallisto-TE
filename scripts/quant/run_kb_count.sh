#!/usr/bin/env bash
# =============================================================================
# run_kb_count.sh  --  count genes + TEs per cell in droplet single-cell RNA-seq
#
# For every sample in a folder of RAW (untrimmed) FASTQ files this script:
#   1. finds read 1 (barcode + UMI) and read 2 (cDNA) of every lane, skips the
#      index reads (I1 / I2), and combines the lanes of each sample
#   2. checks the reads fit the chosen technology (read 1 long enough, not
#      swapped, not trimmed)
#   3. runs kb count against the genes + TEs index from step 1
#   4. prints a summary (reads, % pseudoaligned, % valid barcodes, cells,
#      % of UMIs from TEs)
#
# Usage:
#   run_kb_count.sh -i FASTQ_DIR --index INDEX.idx -x TECH -o OUT_DIR [options]
#
# Required:
#   -i, --input        folder with the RAW FASTQ files, or one sub-folder per sample
#       --index        the .idx file from step 1, e.g. ref/index/mm10_genes_TE.idx
#   -x, --tech         single-cell technology, as in `kb --list`:
#                        10xv2  10x Chromium 3' v2      (read 1 = 26 bp)
#                        10xv3  10x Chromium 3' v3/v3.1 (read 1 = 28 bp)
#                        10xv4  10x GEM-X 3' v4         (read 1 = 28 bp)
#                        dropseq, indropsv3, bdwta, ...
#   -o, --outdir       output folder (created if needed)
#
# Options:
#       --t2g FILE     transcript-to-gene table (default: t2g.txt next to the index)
#   -t, --threads N    threads (default 8)
#   -m, --memory M     memory for sorting barcodes (default 16G)
#       --filter       also write counts_filtered/ (cells chosen by the bustools
#                      knee-plot filter; you can also filter later in Seurat)
#       --h5ad         also write an AnnData .h5ad file (for scanpy users)
#       --force        re-run samples that already have results
#   -h, --help         show this help
#
# File names it understands (sample name = everything before _S1_L001):
#   Sample1_S1_L001_R1_001.fastq.gz + Sample1_S1_L001_R2_001.fastq.gz  (+ more lanes)
#   Sample1_R1.fastq.gz + Sample1_R2.fastq.gz,  Sample1_1.fastq.gz + Sample1_2.fastq.gz
#
# Outputs:
#   OUT_DIR/<sample>/counts_unfiltered/cells_x_genes.mtx   UMI counts (cells x genes/TEs)
#   OUT_DIR/<sample>/counts_unfiltered/cells_x_genes.barcodes.txt
#   OUT_DIR/<sample>/counts_unfiltered/cells_x_genes.genes.txt        gene / TE IDs
#   OUT_DIR/<sample>/counts_unfiltered/cells_x_genes.genes.names.txt  names ("TE:..." for TEs)
#   OUT_DIR/<sample>/run_info.json, inspect.json, kb_info.json
#   OUT_DIR/kb_count_summary.tsv   one line per sample
#
# NOTE: give kb count the RAW reads. Never trim single-cell reads: read 1 holds
#   the cell barcode + UMI at fixed positions.
# =============================================================================
set -euo pipefail

usage() { awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"; exit "${1:-0}"; }
die()   { echo "ERROR: $*" >&2; exit 1; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IN="" INDEX="" TECH="" OUT="" T2G="" THREADS=8 MEM=16G FILTER=0 H5AD=0 FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)    IN="$2"; shift 2 ;;
        --index)       INDEX="$2"; shift 2 ;;
        -x|--tech)     TECH="$2"; shift 2 ;;
        -o|--outdir)   OUT="$2"; shift 2 ;;
        --t2g)         T2G="$2"; shift 2 ;;
        -t|--threads)  THREADS="$2"; shift 2 ;;
        -m|--memory)   MEM="$2"; shift 2 ;;
        --filter)      FILTER=1; shift ;;
        --h5ad)        H5AD=1; shift ;;
        --force)       FORCE=1; shift ;;
        -h|--help)     usage 0 ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done
[[ -n "$IN" && -n "$INDEX" && -n "$TECH" && -n "$OUT" ]] || usage 1
[[ -d "$IN" ]]    || die "input folder not found: $IN"
[[ -s "$INDEX" ]] || die "index file not found: $INDEX (it is the .idx file in the index/ folder from step 1)"
[[ -n "$T2G" ]] || T2G="$(dirname "$INDEX")/t2g.txt"
[[ -s "$T2G" ]]   || die "t2g.txt not found next to the index; give it with --t2g FILE"
command -v kb >/dev/null 2>&1 || die "kb not found. Activate the environment first: conda activate kallisto_te"
TECH_UC=$(tr '[:lower:]' '[:upper:]' <<<"$TECH")
TECH_LIST=$(kb --list 2>&1 || true)
awk '{print $1}' <<<"$TECH_LIST" | grep -qx "$TECH_UC" \
    || die "unknown technology '$TECH'. Run 'kb --list' to see the supported names (e.g. 10xv3)"
case "$TECH_UC" in BULK|SMARTSEQ2|SMARTSEQ3|STORMSEQ)
    die "'$TECH' is not a droplet technology; for bulk RNA-seq use scripts/quant/run_kallisto_quant.sh" ;;
esac

mkdir -p "$OUT/logs"
OUT="$(cd "$OUT" && pwd)"
IN="$(cd "$IN" && pwd)"
INDEX="$(cd "$(dirname "$INDEX")" && pwd)/$(basename "$INDEX")"
T2G="$(cd "$(dirname "$T2G")" && pwd)/$(basename "$T2G")"

# shellcheck source=../lib/fastq_pairs.sh
source "$HERE/../lib/fastq_pairs.sh"

# ----------------------------- find the samples ------------------------------
collect_samples "$IN" 0 ""
[[ ${#S_NAME[@]} -gt 0 ]] || die "no FASTQ files (.fq.gz/.fastq.gz/.fq/.fastq) in $IN or its sub-folders"
if [[ ${#SKIPPED[@]} -gt 0 ]]; then
    echo "NOTE: ignoring ${#SKIPPED[@]} Illumina index-read file(s) (*_I1_001 / *_I2_001): they hold sample barcodes, not RNA" >&2
fi
[[ ${#NOT_PAIRS[@]} -eq 0 ]] || die "these files look like read 1 + read 2 but their read names differ:
$(printf '         %s\n' "${NOT_PAIRS[@]}")
       Single-cell samples need matching read 1 + read 2 files from the same run."

echo "Found ${#S_NAME[@]} sample(s), technology $TECH_UC:" >&2
for ((i = 0; i < ${#S_NAME[@]}; i++)); do
    IFS=$'\t' read -r -a fa <<<"${S_FILES[i]}"
    files=""
    for f in "${fa[@]}"; do files+="${files:+ }$(basename "$f")"; done
    printf '  %-24s %s\n' "${S_NAME[i]}" "$files" >&2
done

# ----------------------------- check the reads -------------------------------
# Read 1 must be at least barcode + UMI long for the chosen technology.
case "$TECH_UC" in
    10XV2) need=26 ;; 10XV3|10XV4|VISIUM) need=28 ;; DROPSEQ) need=20 ;; *) need=0 ;;
esac
for ((i = 0; i < ${#S_NAME[@]}; i++)); do
    name="${S_NAME[i]}"
    [[ "${S_LAYOUT[i]}" == paired ]] || die "sample $name has no read 2 file. Single-cell data needs read 1 (barcode + UMI)
       and read 2 (cDNA) of every lane; files found: $(tr '\t' '\n' <<<"${S_FILES[i]}" | xargs -n1 basename | tr '\n' ' ')"
    [[ "${S_FILES[i]}" != *".trimmed."* ]] || die "sample $name uses trimmed files (*.trimmed.*). Single-cell reads must NOT be
       trimmed - give kb count the raw FASTQ files."
    IFS=$'\t' read -r r1 r2 _ <<<"${S_FILES[i]}"
    l1=$(max_len "$r1"); l2=$(max_len "$r2")
    if (( l2 > 0 && l2 <= 30 && l1 >= 2 * l2 )); then
        die "sample $name: read 1 is ${l1} bp but read 2 is only ${l2} bp. Read 1 and read 2 look swapped:
       read 1 should be the short barcode + UMI read."
    fi
    if (( l1 < need )); then
        die "sample $name: read 1 is only ${l1} bp, but $TECH_UC needs at least $need bp (cell barcode + UMI).
       Check the chemistry (10x v2 = 26 bp, v3/v4 = 28 bp) and that the reads were not trimmed."
    fi
done

# ----------------------------- kb count --------------------------------------
for ((i = 0; i < ${#S_NAME[@]}; i++)); do
    name="${S_NAME[i]}"
    IFS=$'\t' read -r -a fa <<<"${S_FILES[i]}"
    if [[ $FORCE -eq 0 && -s "$OUT/$name/counts_unfiltered/cells_x_genes.mtx" ]]; then
        echo "  [skip] $name (already counted; use --force to redo)" >&2
        continue
    fi
    echo "  [kb count] $name (${#fa[@]} files)" >&2
    rm -rf "${OUT:?}/$name"
    args=(count -i "$INDEX" -g "$T2G" -x "$TECH_UC" -o "$OUT/$name" -t "$THREADS" -m "$MEM")
    [[ $FILTER -eq 1 ]] && args+=(--filter bustools)
    [[ $H5AD -eq 1 ]] && args+=(--h5ad)
    kb "${args[@]}" "${fa[@]}" > "$OUT/logs/$name.log" 2>&1 \
        || die "kb count failed for $name - see $OUT/logs/$name.log"
done

# ----------------------------- summary ---------------------------------------
python3 - "$OUT" "$T2G" "${S_NAME[@]}" <<'EOF'
import json, sys
from pathlib import Path
out, t2g, names = Path(sys.argv[1]), sys.argv[2], sys.argv[3:]
te_ids = set()
for line in open(t2g):
    f = line.rstrip("\n").split("\t")
    if len(f) > 2 and f[2].startswith("TE:"):
        te_ids.add(f[1])
rows = ["sample\treads_processed\tpercent_pseudoaligned\tpercent_valid_barcode\tbarcodes"
        "\tbarcodes_500plus_UMIs\tmedian_UMIs_500plus\tUMIs_total\tpercent_UMIs_TE"]
for n in names:
    d = out / n
    info = json.load(open(d / "run_info.json"))
    try:
        onlist = f'{json.load(open(d / "inspect.json"))["percentageReadsOnOnlist"]:.1f}'
    except (OSError, KeyError, ValueError):
        onlist = "NA"
    cu = d / "counts_unfiltered"
    genes = [g.strip() for g in open(cu / "cells_x_genes.genes.txt")]
    is_te = [g in te_ids for g in genes]
    per_cell, total, te = {}, 0.0, 0.0
    with open(cu / "cells_x_genes.mtx") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            break                                   # the size line
        for line in fh:                             # cell gene count (1-based)
            c, g, v = line.split()
            v = float(v)
            per_cell[c] = per_cell.get(c, 0.0) + v
            total += v
            if is_te[int(g) - 1]:
                te += v
    big = sorted(v for v in per_cell.values() if v >= 500)
    med = big[len(big) // 2] if big else 0
    rows.append("\t".join(map(str, [
        n, info["n_processed"], f'{info["p_pseudoaligned"]:.1f}', onlist, len(per_cell), len(big),
        f"{med:.0f}", f"{total:.0f}", f"{100 * te / total:.2f}" if total else "NA"])))
open(out / "kb_count_summary.tsv", "w").write("\n".join(rows) + "\n")
EOF

echo >&2
echo "Summary (also in $OUT/kb_count_summary.tsv):" >&2
column -t -s $'\t' "$OUT/kb_count_summary.tsv" >&2 || cat "$OUT/kb_count_summary.tsv" >&2
echo >&2
low=$(awk -F'\t' 'NR > 1 && $4 != "NA" && $4 + 0 < 70 { print "  " $1 " (" $4 "%)" }' "$OUT/kb_count_summary.tsv")
if [[ -n "$low" ]]; then
    {
        echo "WARNING: under 70% of reads carry a valid cell barcode in:"
        echo "$low"
        echo "  Usually the wrong --tech (e.g. 10xv3 vs 10xv4 use different barcode lists)."
        echo "  See docs/05_singlecell_kb_count.md, section 2."
        echo
    } >&2
fi
echo "Count matrices: $OUT/<sample>/counts_unfiltered/  (load into Seurat - see docs/05_singlecell_kb_count.md)" >&2
