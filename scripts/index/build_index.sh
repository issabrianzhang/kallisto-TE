#!/usr/bin/env bash
# =============================================================================
# build_index.sh  --  build a kallisto | bustools index of genes + TEs
#
# Combines a reference genome + gene annotation (GTF) with a library of
# transposable element (TE) consensus sequences, so that one pseudoalignment
# quantifies canonical genes AND TE subfamilies together.
#
#   genome.fa  +  TE_library.fa       ->  genome_plus_TE.fa
#   genes.gtf  +  TE_library.pseudo.gtf ->  genes_plus_TE.gtf
#   kb ref (genome_plus_TE.fa, genes_plus_TE.gtf) -> index.idx, t2g.txt, cdna.fa
#
# -----------------------------------------------------------------------------
# PRESETS (genome + GTF are downloaded from GENCODE automatically)
#   mouse : GRCm38/mm10, GENCODE M25 primary assembly annotation
#           RepBase libraries: mousub.ref + rodrep.ref
#   human : GRCh38/hg38, GENCODE v40 basic annotation
#           RepBase libraries: humrep.ref + humsub.ref
#
# CUSTOM (any species): supply your own genome, GTF and TE library
#   --species custom --genome X.fa[.gz] --gtf X.gtf[.gz] \
#       ( --te-fasta TEs.fa [--te-fasta more.fa]
#         | --repbase RepBaseXX.fasta.tar.gz --te-libs "abc.ref xyz.ref" )
#
# -----------------------------------------------------------------------------
# Usage:
#   build_index.sh --species mouse|human|custom --outdir DIR [options]
#
# Required:
#   -s, --species      mouse | human | custom
#   -o, --outdir       output directory (created if needed)
#
# TE library (one of):
#   -r, --repbase      RepBase archive (RepBaseXX.XX.fasta.tar.gz) or the
#                      folder it extracts to
#   -l, --te-libs      RepBase *.ref files to use, quoted & space separated
#                      (defaults: mouse "mousub.ref rodrep.ref",
#                                 human "humrep.ref humsub.ref")
#   -e, --te-fasta     any TE consensus FASTA (repeatable); use instead of,
#                      or in addition to, --repbase
#
# Genome / annotation (required for custom; optional overrides for presets):
#   -g, --genome       genome FASTA (.fa or .fa.gz)
#   -a, --gtf          gene annotation GTF (.gtf or .gtf.gz)
#
# Options:
#   -n, --name         prefix for output files (default: species name)
#   -t, --threads      threads for kb ref (default: 8)
#   -d, --d-list       none | genome  (default: none; see docs/01_build_index.md)
#       --tmp          temporary directory for kb ref (must NOT exist yet;
#                      default: OUTDIR/tmp_kb_<pid>)
#   -h, --help         show this help
#
# Examples:
#   build_index.sh -s mouse -r RepBase25.05.fasta.tar.gz -o ref/mm10_TE -t 16
#   build_index.sh -s human -r RepBase28.04.fasta.tar.gz -o ref/hg38_TE -t 16
#   build_index.sh -s custom -n danRer11 -g danRer11.fa.gz -a danRer11.gtf.gz \
#       -e zebrafish_TEs.fa -o ref/danRer11_TE
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
usage() { awk 'NR>2 && /^# ====/{exit} NR>2{sub(/^# ?/,""); print}' "$0"; exit "${1:-0}"; }
log()   { printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()   { echo "ERROR: $*" >&2; exit 1; }

# ------------------------------- arguments ----------------------------------
SPECIES="" OUTDIR="" REPBASE="" TE_LIBS="" GENOME="" GTF="" NAME=""
THREADS=8 DLIST="none" TMPDIR_KB=""
TE_FASTAS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--species)  SPECIES="$2"; shift 2 ;;
        -o|--outdir)   OUTDIR="$2"; shift 2 ;;
        -r|--repbase)  REPBASE="$2"; shift 2 ;;
        -l|--te-libs)  TE_LIBS="$2"; shift 2 ;;
        -e|--te-fasta) TE_FASTAS+=("$2"); shift 2 ;;
        -g|--genome)   GENOME="$2"; shift 2 ;;
        -a|--gtf)      GTF="$2"; shift 2 ;;
        -n|--name)     NAME="$2"; shift 2 ;;
        -t|--threads)  THREADS="$2"; shift 2 ;;
        -d|--d-list)   DLIST="$2"; shift 2 ;;
        --tmp)         TMPDIR_KB="$2"; shift 2 ;;
        -h|--help)     usage 0 ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

[[ -n "$SPECIES" && -n "$OUTDIR" ]] || usage 1
NAME="${NAME:-$SPECIES}"

# ------------------------------- presets ------------------------------------
GENCODE_GENOME_URL="" GENCODE_GTF_URL=""
case "$SPECIES" in
    mouse)
        TE_LIBS="${TE_LIBS:-mousub.ref rodrep.ref}"
        GENCODE_GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz"
        GENCODE_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz"
        ;;
    human)
        TE_LIBS="${TE_LIBS:-humrep.ref humsub.ref}"
        GENCODE_GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz"
        GENCODE_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
        ;;
    custom)
        [[ -n "$GENOME" && -n "$GTF" ]] || die "--species custom needs --genome and --gtf"
        if [[ -n "$REPBASE" && -z "$TE_LIBS" ]]; then
            die "--species custom with --repbase needs --te-libs (e.g. \"zebrep.ref\")"
        fi
        ;;
    *) die "--species must be mouse, human or custom" ;;
esac
[[ -n "$REPBASE" || ${#TE_FASTAS[@]} -gt 0 ]] || die "give a TE library with --repbase and/or --te-fasta"
[[ "$DLIST" == "none" || "$DLIST" == "genome" ]] || die "--d-list must be none or genome"
command -v kb >/dev/null 2>&1 || die "kb not found. Activate the environment first: conda activate kallisto_te"

mkdir -p "$OUTDIR"/{te,reference,combined,index}
OUTDIR="$(cd "$OUTDIR" && pwd)"
TE_FA="$OUTDIR/te/${NAME}_TE_library.fa"
TE_GTF="$OUTDIR/te/${NAME}_TE_library.pseudo.gtf"
COMB_FA="$OUTDIR/combined/${NAME}_genome_plus_TE.fa"
COMB_GTF="$OUTDIR/combined/${NAME}_genes_plus_TE.gtf"
INDEX="$OUTDIR/index/${NAME}_genes_TE.idx"
T2G="$OUTDIR/index/t2g.txt"
CDNA="$OUTDIR/index/cdna.fa"
TMPDIR_KB="${TMPDIR_KB:-$OUTDIR/tmp_kb_$$}"

{
    echo "kallisto-TE build_index.sh"
    echo "date:     $(date)"
    echo "species:  $SPECIES  (name: $NAME)"
    echo "repbase:  ${REPBASE:-<none>}  libs: ${TE_LIBS:-<none>}"
    echo "te-fasta: ${TE_FASTAS[*]:-<none>}"
    echo "genome:   ${GENOME:-$GENCODE_GENOME_URL}"
    echo "gtf:      ${GTF:-$GENCODE_GTF_URL}"
    echo "d-list:   $DLIST"
    echo "kb:       $(kb --version 2>&1 | grep -i kb_python || true)"
} | tee "$OUTDIR/build_info.txt" >&2

# ---------------------- step 1: assemble the TE library ----------------------
log "Step 1/5  Preparing the TE consensus library"
TE_INPUTS=()
if [[ -n "$REPBASE" ]]; then
    if [[ -f "$REPBASE" ]]; then
        log "  extracting $REPBASE"
        mkdir -p "$OUTDIR/te/repbase_raw"
        tar -xzf "$REPBASE" -C "$OUTDIR/te/repbase_raw"
        RB_DIR="$OUTDIR/te/repbase_raw"
    elif [[ -d "$REPBASE" ]]; then
        RB_DIR="$REPBASE"
    else
        die "--repbase not found: $REPBASE"
    fi
    for lib in $TE_LIBS; do
        hit=$(find "$RB_DIR" -type f -name "$lib" | head -1)
        [[ -n "$hit" ]] || die "RepBase library '$lib' not found under $RB_DIR (ls it to see the available *.ref files)"
        echo "  using $hit" >&2
        TE_INPUTS+=("$hit")
    done
fi
if [[ ${#TE_FASTAS[@]} -gt 0 ]]; then
    TE_INPUTS+=("${TE_FASTAS[@]}")
fi
bash "$SCRIPT_DIR/prepare_te_library.sh" -o "$TE_FA" "${TE_INPUTS[@]}"

# ------------------------ step 2: TE pseudo-GTF ------------------------------
log "Step 2/5  Writing the TE pseudo-GTF"
bash "$SCRIPT_DIR/make_te_pseudo_gtf.sh" "$TE_FA" "$TE_GTF"

# ------------------- step 3: genome + annotation -----------------------------
log "Step 3/5  Getting the genome FASTA and gene GTF"
fetch() {  # fetch URL DEST
    if [[ -s "$2" ]]; then echo "  already downloaded: $(basename "$2")" >&2; return; fi
    echo "  downloading $1" >&2
    if command -v wget >/dev/null 2>&1; then wget -q -O "$2.part" "$1"
    else curl -fsSL -o "$2.part" "$1"; fi
    mv "$2.part" "$2"
}
if [[ -z "$GENOME" ]]; then
    GENOME="$OUTDIR/reference/$(basename "$GENCODE_GENOME_URL")"
    fetch "$GENCODE_GENOME_URL" "$GENOME"
fi
if [[ -z "$GTF" ]]; then
    GTF="$OUTDIR/reference/$(basename "$GENCODE_GTF_URL")"
    fetch "$GENCODE_GTF_URL" "$GTF"
fi
[[ -s "$GENOME" ]] || die "genome not found: $GENOME"
[[ -s "$GTF" ]]    || die "GTF not found: $GTF"

# ------------------- step 4: combine genome/GTF with TEs ---------------------
log "Step 4/5  Appending TE consensus sequences to the genome and GTF"
catx() { case "$1" in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac; }

# TE names must not collide with chromosome/contig names in the genome.
catx "$GENOME" | grep '^>' | sed 's/^>//; s/[[:space:]].*//' | sort -u > "$OUTDIR/combined/.genome_contigs"
grep '^>' "$TE_FA" | sed 's/^>//' | sort -u > "$OUTDIR/combined/.te_names"
clash=$(comm -12 "$OUTDIR/combined/.genome_contigs" "$OUTDIR/combined/.te_names" | head -5)
[[ -z "$clash" ]] || die "TE names also used as genome contig names: $clash"
rm -f "$OUTDIR/combined/.genome_contigs" "$OUTDIR/combined/.te_names"

catx "$GENOME" > "$COMB_FA"
[[ "$(tail -c1 "$COMB_FA" | od -An -c | tr -d ' ')" == '\n' ]] || echo >> "$COMB_FA"
cat "$TE_FA" >> "$COMB_FA"

catx "$GTF" > "$COMB_GTF"
[[ "$(tail -c1 "$COMB_GTF" | od -An -c | tr -d ' ')" == '\n' ]] || echo >> "$COMB_GTF"
cat "$TE_GTF" >> "$COMB_GTF"
echo "  $COMB_FA" >&2
echo "  $COMB_GTF" >&2

# ------------------------- step 5: kb ref ------------------------------------
log "Step 5/5  Building the kallisto index with kb ref (this can take a few hours)"
[[ ! -e "$TMPDIR_KB" ]] || die "temporary directory already exists: $TMPDIR_KB (kb ref needs a fresh one)"
DLIST_ARGS=()
[[ "$DLIST" == "none" ]] && DLIST_ARGS=(--d-list None)

kb ref --verbose --overwrite \
    -t "$THREADS" \
    --tmp "$TMPDIR_KB" \
    "${DLIST_ARGS[@]+"${DLIST_ARGS[@]}"}" \
    -i "$INDEX" \
    -g "$T2G" \
    -f1 "$CDNA" \
    "$COMB_FA" \
    "$COMB_GTF" 2>&1 | tee "$OUTDIR/index/kb_ref.log" >&2

log "Checking that every TE made it into the index"
bash "$SCRIPT_DIR/check_index.sh" "$TE_FA" "$T2G" "$CDNA" | tee -a "$OUTDIR/build_info.txt" >&2

log "Done. Use these files for quantification:"
echo "  index : $INDEX"   >&2
echo "  t2g   : $T2G"     >&2
echo "  (the combined genome FASTA can be deleted to save space: $COMB_FA)" >&2
