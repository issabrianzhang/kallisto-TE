# shellcheck shell=bash
# =============================================================================
# fastq_pairs.sh  --  shared helpers: find FASTQ files and pair read 1 with read 2
#
# Sourced by scripts/qc/run_fastp.sh, scripts/quant/run_kallisto_quant.sh and
# scripts/quant/run_kb_count.sh. Not meant to be run on its own.
#
# find_fastq_units DIR SINGLE R1TAG
#   DIR     folder with FASTQ files (not searched recursively)
#   SINGLE  1 = every file is its own single-end unit, 0 = pair by name
#   R1TAG   read-1 tag ("_R1", "_1", ...); empty = try _R1_001, _R1, _1 per file
# sets these global arrays / variables:
#   FASTQ_FILES      every FASTQ file in DIR
#   FASTQ_UNITS      "name<TAB>R1<TAB>R2" per unit (R2 empty for single-end)
#   FASTQ_NOT_PAIRS  "a.fq.gz + b.fq.gz" for files that looked like a pair by
#                    name but whose read names differ (so are NOT mates)
#   FASTQ_R1TAG      the read-1 tag(s) used for pairs, comma-separated ("" if none)
#   FASTQ_SKIPPED    Illumina index-read files (*_I1_001.*, *_I2_001.*), ignored
#
# collect_samples DIR SINGLE R1TAG
#   groups the units of find_fastq_units into samples (needs a die() function):
#     - FASTQs directly in DIR: lanes of one sample are combined
#       (Ctrl_S1_L001_R1_001 + Ctrl_S1_L002_R1_001 -> sample "Ctrl")
#     - no FASTQs directly in DIR: each sub-folder is one sample
#   sets parallel arrays (bash 3.2 has no associative arrays):
#     S_NAME[i]   sample name      S_LAYOUT[i]  paired | single
#     S_FILES[i]  its files in kallisto order, TAB-separated (R1 R2 R1 R2 ...)
#     NOT_PAIRS   / SKIPPED: FASTQ_NOT_PAIRS / FASTQ_SKIPPED of every folder
#
# File names understood (unit name = everything before the tag):
#   Sample1_R1.fq.gz + Sample1_R2.fq.gz
#   Sample1_R1_001.fastq.gz + Sample1_R2_001.fastq.gz      (Illumina)
#   Sample1_1.fq.gz + Sample1_2.fq.gz                      (SRA / ENA)
#   Sample1_R1.trimmed.fq.gz + Sample1_R2.trimmed.fq.gz    (output of step 2)
#   Sample1.fq.gz, or an R1 file without its R2            (single-end)
# =============================================================================

# Longest suffix first, so "x_R1.trimmed.fq.gz" gives the name "x_R1".
FASTQ_EXTS=".trimmed.fastq.gz .trimmed.fq.gz .trimmed.fastq .trimmed.fq .fastq.gz .fq.gz .fastq .fq"

ext_of()   { local e; for e in $FASTQ_EXTS; do [[ "$1" == *"$e" ]] && { echo "$e"; return; }; done; }
mate_tag() { case "$1" in *R1_001) echo "${1%R1_001}R2_001" ;; *R1) echo "${1%R1}R2" ;; *1) echo "${1%1}2" ;; esac; }
cat_fastq() { case "$1" in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac; }

# Read names of the first 5 reads (first word of the header, without @ and
# without a /1 or /2 suffix). Mates of a real pair give identical output.
read_ids() {
    cat_fastq "$1" 2>/dev/null | head -n 20 \
        | awk 'NR % 4 == 1 { id = $1; sub(/^@/, "", id); sub(/\/[12]$/, "", id); print id }' || true
}

# Longest read among the first 2,000 reads.
max_len() {
    cat_fastq "$1" 2>/dev/null | head -n 8000 \
        | awk 'NR % 4 == 2 { if (length($0) > m) m = length($0) } END { print m + 0 }' || true
}

find_fastq_units() {
    local dir="$1" single="$2" r1tag="$3" r2tag="" e f b t name r2 s r1
    FASTQ_FILES=() FASTQ_UNITS=() FASTQ_NOT_PAIRS=() FASTQ_R1TAG="" FASTQ_SKIPPED=()

    local had_nullglob=0
    shopt -q nullglob && had_nullglob=1
    shopt -s nullglob
    local seen=" "
    for e in $FASTQ_EXTS; do
        for f in "$dir"/*"$e"; do
            # ".fq.gz" also matches "x.trimmed.fq.gz": keep each file once
            [[ "$seen" == *" $f "* ]] && continue
            seen+="$f "
            # Illumina index reads (sample barcodes), not RNA: never use them
            if [[ "$(basename "$f" "$e")" =~ _I[12]_001$ ]]; then
                FASTQ_SKIPPED+=("$f"); continue
            fi
            FASTQ_FILES+=("$f")
        done
    done
    [[ $had_nullglob -eq 1 ]] || shopt -u nullglob
    [[ ${#FASTQ_FILES[@]} -gt 0 ]] || return 0

    # Read-1 tags to try, file by file (so one folder may mix naming styles).
    local tags=() used=" " t2 paired
    if [[ "$single" -eq 1 ]]; then
        tags=()
    elif [[ -n "$r1tag" ]]; then
        tags=("$r1tag")
    else
        tags=(_R1_001 _R1 _1)
    fi

    local units=()
    for f in "${FASTQ_FILES[@]}"; do
        e=$(ext_of "$f"); b=$(basename "$f" "$e")
        paired=0
        for t in ${tags[@]+"${tags[@]}"}; do          # read 1 with its read 2 next to it
            [[ "$b" == *"$t" ]] || continue
            name="${b%"$t"}"; r2="$dir/${name}$(mate_tag "$t")${e}"
            if [[ -f "$r2" ]]; then
                units+=("$name"$'\t'"$f"$'\t'"$r2"); used+="$t "; paired=1; break
            fi
        done
        [[ $paired -eq 1 ]] && continue
        for t in ${tags[@]+"${tags[@]}"}; do          # read 2: handled with its read 1
            t2=$(mate_tag "$t")
            [[ "$b" == *"$t2" && -f "$dir/${b%"$t2"}${t}${e}" ]] && continue 2
        done
        # single-end. A lone read-1 file loses an unambiguous tag (Sample_R1 ->
        # Sample); "_1" is kept, as it often numbers replicates.
        name="$b"
        for t in ${tags[@]+"${tags[@]}"}; do
            [[ "$t" == _1 && -z "$r1tag" ]] && continue
            [[ "$b" == *"$t" ]] && { name="${b%"$t"}"; break; }
        done
        units+=("$name"$'\t'"$f"$'\t')
    done
    used="${used# }"; used="${used% }"
    FASTQ_R1TAG=$(tr ' ' '\n' <<<"$used" | awk 'NF && !seen[$0]++' | paste -sd, -)

    # Mates of a real pair carry the same read names; look-alikes do not.
    for s in "${units[@]}"; do
        IFS=$'\t' read -r name r1 r2 <<<"$s"
        if [[ -n "$r2" && "$(read_ids "$r1")" != "$(read_ids "$r2")" ]]; then
            for f in "$r1" "$r2"; do
                FASTQ_UNITS+=("$(basename "$f" "$(ext_of "$f")")"$'\t'"$f"$'\t')
            done
            FASTQ_NOT_PAIRS+=("$(basename "$r1") + $(basename "$r2")")
        else
            FASTQ_UNITS+=("$s")
        fi
    done
}

S_NAME=() S_LAYOUT=() S_FILES=() NOT_PAIRS=() SKIPPED=()

_add_unit() {   # sample-name R1 [R2]
    local name="$1" r1="$2" r2="${3:-}" layout=single i
    [[ -n "$r2" ]] && layout=paired
    for ((i = 0; i < ${#S_NAME[@]}; i++)); do
        if [[ "${S_NAME[i]}" == "$name" ]]; then
            [[ "${S_LAYOUT[i]}" == "$layout" ]] \
                || die "sample $name has both paired-end and single-end files; put them in separate folders"
            S_FILES[i]+=$'\t'"$r1${r2:+$'\t'$r2}"
            return 0
        fi
    done
    S_NAME+=("$name"); S_LAYOUT+=("$layout"); S_FILES+=("$r1${r2:+$'\t'$r2}")
}

_collect_dir() {   # DIR SINGLE R1TAG [sample-name for every unit in DIR]
    local dir="$1" fixed="${4:-}" s name r1 r2
    find_fastq_units "$dir" "$2" "$3"
    [[ ${#FASTQ_SKIPPED[@]} -eq 0 ]] || SKIPPED+=("${FASTQ_SKIPPED[@]}")
    [[ ${#FASTQ_UNITS[@]} -gt 0 ]] || return 0
    [[ ${#FASTQ_NOT_PAIRS[@]} -eq 0 ]] || NOT_PAIRS+=("${FASTQ_NOT_PAIRS[@]}")
    for s in "${FASTQ_UNITS[@]}"; do
        IFS=$'\t' read -r name r1 r2 <<<"$s"
        if [[ -n "$fixed" ]]; then
            name="$fixed"
        else    # lanes of one sample: Ctrl_S1_L001 / Ctrl_S1_L002 -> Ctrl
            name=$(sed -E 's/(_S[0-9]+)?_L[0-9]{3}$//' <<<"$name")
        fi
        _add_unit "$name" "$r1" "$r2"
    done
}

collect_samples() {   # DIR SINGLE R1TAG
    local d
    S_NAME=() S_LAYOUT=() S_FILES=() NOT_PAIRS=() SKIPPED=()
    _collect_dir "$1" "$2" "$3"
    if [[ ${#S_NAME[@]} -eq 0 ]]; then        # no FASTQs at the top: one sub-folder per sample
        for d in "$1"/*/; do
            [[ -d "$d" ]] && _collect_dir "${d%/}" "$2" "$3" "$(basename "$d")"
        done
    fi
    return 0
}
