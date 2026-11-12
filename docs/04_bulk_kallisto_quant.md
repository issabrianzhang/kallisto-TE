# Step 3 — Count genes and TEs in bulk RNA-seq (kallisto quant)

This step counts, for every sample, how many reads come from each gene and
each TE subfamily. It uses the index from step 1 and the trimmed reads from
step 2. No prior experience is assumed. One command does everything:

```bash
bash scripts/quant/run_kallisto_quant.sh \
    --input  /path/to/project/01_trimmed/trimmed \
    --index  /path/to/references/mm10_TE/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx \
    --outdir /path/to/project/02_quant \
    --threads 8
```

> **Single-cell RNA-seq (10x Genomics, Drop-seq, ...)?** Use
> [step 4 (`kb count`)](05_singlecell_kb_count.md) instead. This script
> recognises droplet reads (a short barcode read 1 with a long read 2) and
> stops with a pointer to step 4.

---

## Contents

1. [What this step does](#1-what-this-step-does)
2. [What you need](#2-what-you-need)
3. [Three settings to know](#3-three-settings-to-know)
4. [Run it](#4-run-it)
5. [What you get](#5-what-you-get)
6. [Read the results](#6-read-the-results)
7. [Next: differential expression with DESeq2](#7-next-differential-expression-with-deseq2)
8. [If something looks wrong](#8-if-something-looks-wrong)
9. [For your methods section](#9-for-your-methods-section)

---

## 1. What this step does

[kallisto](https://pachterlab.github.io/kallisto/manual) does not align reads
base by base. It looks up the short words (k-mers) of each read in the index
and finds the set of transcripts the read is compatible with
("pseudoalignment"). This is fast: a sample takes minutes, not hours.

Your index contains every annotated transcript **plus one consensus sequence
per TE subfamily**, so a read from an `MMERVK10C` copy anywhere in the genome
is counted towards `MMERVK10C`. A read that fits several transcripts or
several related TE subfamilies is shared between them by kallisto's
expectation-maximisation (EM) step, in proportion to how much evidence each
one has.

For each sample kallisto reports, per transcript and per TE:

| Column | Meaning |
|--------|---------|
| `est_counts` | estimated number of reads (can be fractional, because shared reads are split) |
| `tpm` | transcripts per million: expression corrected for length and sequencing depth; use it to compare genes *within* a sample |
| `length` / `eff_length` | sequence length, and the length usable by fragments |

The script then adds up the transcripts of each gene, so you get one row per
gene and one row per TE subfamily.

## 2. What you need

| What | Where it comes from |
|------|---------------------|
| Reads | `01_trimmed/trimmed/` from [step 2](03_qc_trimming.md). Raw FASTQ files also work |
| Index | `<outdir>/index/<name>_genes_TE.idx` from [step 1](01_build_index.md) |
| `t2g.txt` | the transcript-to-gene table, found automatically next to the index |
| Software | the `kallisto_te` conda environment (same kb-python version that built the index) |

Use **one index for every sample of a project**. Counts from different
indices cannot be compared.

## 3. Three settings to know

### Strandedness (`--strand`)

Most RNA-seq kits keep track of which DNA strand the RNA came from. This
matters for TEs, because many are transcribed in the antisense direction
inside genes. Telling kallisto the right strand keeps those reads apart.

| Your library kit | `--strand` |
|------------------|------------|
| Illumina TruSeq Stranded mRNA / Total RNA, NEBNext Ultra II **Directional**, KAPA Stranded, any "dUTP" method | `rf` (reverse-stranded; the most common) |
| Lexogen QuantSeq FWD, Illumina ScriptSeq | `fr` |
| NEBNext Ultra II (non-directional), Smart-seq2, Illumina TruSeq RNA v2 (non-stranded) | `unstranded` |

Not sure? Leave the default, **`--strand auto`**. The script takes the
first 200,000 reads of the first sample, quantifies them once as `fr` and once
as `rf`, and compares how many reads fit each way (example output):

```
Checking strandedness on the first 200,000 reads of Ctrl_Rep1 ...
  fr-stranded: 2112 reads pseudoaligned
  rf-stranded: 151630 reads pseudoaligned
  -> reverse-stranded (99% of strand-matched reads fit 'rf', 1% fit 'fr')
```

* ≥ 80% fit one way: a stranded library, and that orientation is used.
* About 50/50: unstranded.
* Anything in between: the script stops and asks you to choose with
  `--strand`. Check the library kit in your sequencing facility's report or
  the GEO / SRA record.

All samples of one run are assumed to share the same library kit. If you
combine projects made with different kits, run them separately.

### Fragment length, single-end reads only (`-l`, `-s`)

With paired-end reads kallisto measures the fragment length itself. With
single-end reads it cannot, so you give an estimate: `-l` is the average
length of the RNA fragment (without adapters) and `-s` its standard
deviation. The defaults, **`-l 200 -s 20`**, suit most standard libraries.
If your facility reports a library size from a Bioanalyzer / TapeStation
(e.g. 320 bp), subtract the adapters (about 120 bp for Illumina TruSeq-style
adapters) to get `-l`.

### Bootstraps (`-b`)

Bootstraps re-run the estimate many times on resampled reads to measure how
certain each count is. They are only used by
[sleuth](https://pachterlab.github.io/sleuth/) and fishpond. The default
`-b 30` keeps that option open. If you will only use DESeq2 or edgeR, `-b 0`
is faster and gives smaller output files.

> Coming from older kallisto scripts? Flags such as `--pseudobam` and
> `--fusion` are not needed for counting, and kallisto 0.50+ no longer has
> them.

## 4. Run it

**Interactively** (a laptop, or an interactive cluster session):

```bash
conda activate kallisto_te
cd /path/to/kallisto-TE

bash scripts/quant/run_kallisto_quant.sh \
    --input   /path/to/project/01_trimmed/trimmed \
    --index   /path/to/references/mm10_TE/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx \
    --outdir  /path/to/project/02_quant \
    --threads 8
```

Single-end reads with a known fragment size, no bootstraps, strand given:

```bash
bash scripts/quant/run_kallisto_quant.sh \
    --input   /path/to/project/01_trimmed/trimmed \
    --index   /path/to/references/mm10_TE/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx \
    --outdir  /path/to/project/02_quant \
    --strand rf -l 200 -s 20 -b 0
```

**On a SLURM cluster**: edit the paths in
[`slurm/run_kallisto_quant.sbatch`](../slurm/run_kallisto_quant.sbatch), then

```bash
mkdir -p logs
sbatch slurm/run_kallisto_quant.sbatch
```

If a job stops halfway, submit it again: samples that already have results
are skipped (use `--force` to redo them).

### How the script finds your samples

The script lists what it found before it starts. Check that list.

| Your files | Sample(s) |
|------------|-----------|
| `Ctrl1_R1.trimmed.fq.gz` + `Ctrl1_R2.trimmed.fq.gz` (step 2 output) | `Ctrl1`, paired-end |
| `Ctrl1.trimmed.fq.gz` | `Ctrl1`, single-end |
| `Ctrl1_S1_L001_R1_001.fastq.gz`, `Ctrl1_S1_L002_R1_001.fastq.gz` + their R2 | `Ctrl1`: **lanes are combined** |
| a folder per sample: `raw/Ctrl1/SRR123_1.fastq` + `SRR123_2.fastq` | `Ctrl1` (the folder name): all files inside are one sample |
| `Ctrl_1.fq.gz` + `Ctrl_2.fq.gz` whose read names differ | two single-end samples (replicates named like a pair; the same check as in step 2) |

Index-read files (`*_I1_001.fastq.gz`, `*_I2_001.fastq.gz`) are ignored. Use
`--single-end` if every file should be its own single-end sample.

## 5. What you get

```
02_quant/
├── Ctrl_Rep1/
│   ├── abundance.tsv            # per transcript and per TE: est_counts, tpm
│   ├── abundance.h5             # the same + bootstraps (for sleuth)
│   └── run_info.json            # reads processed / pseudoaligned
├── Ctrl_Rep2/ ...
├── matrices/
│   ├── gene_TE_counts.tsv       # ← one row per gene or TE, one column per sample
│   ├── gene_TE_tpm.tsv          #   the same, in TPM
│   └── TE_counts.tsv            #   only the TE rows
├── quant_summary.tsv            # ← one line per sample (also printed at the end)
├── samples_used.tsv             # which files went into which sample, layout, strand
├── strand_check/result.txt      # the strandedness test (with --strand auto)
├── logs/<sample>.log            # kallisto's messages per sample
└── multiqc/multiqc_report.html  # pseudoalignment rates of all samples in one page
```

The tables open in Excel, R or Python. Each row of `gene_TE_counts.tsv` is
labelled `gene` or `TE` in the `type` column:

```
gene_id                      gene_name                       type  Ctrl_Rep1  Ctrl_Rep2  KO_Rep1
ENSMUSG00000057666.18        Gapdh                           gene  51234.000  49876.000  50311.000
MMERVK10C_ERVK_Mus_musculus  TE:MMERVK10C_ERVK_Mus_musculus  TE    812.000    790.512    4410.000
```

(numbers made up for illustration). A TE's ID is its name in the TE library:
`NAME_CLASS_Species` for RepBase, `NAME_Type-Subtype` for the bundled Dfam
libraries (e.g. `MMERVK10C_LTR-ERVK`).

## 6. Read the results

`quant_summary.tsv` has one line per sample:

| Column | Meaning | What to look for |
|--------|---------|------------------|
| `reads_processed` | reads (or read pairs) kallisto saw | should match the "reads after" column of step 2 (halved for pairs) |
| `percent_pseudoaligned` | % that matched a gene or TE | often 60–90% for poly(A) mRNA-seq, lower for total RNA (intronic reads). Similar across samples is what matters |
| `TE_reads` / `percent_TE` | estimated reads assigned to TE subfamilies, and their % of all assigned reads | varies with species, tissue and library type; compare samples with each other, not with a fixed number |

The script warns about samples under 30% pseudoaligned. The most common
reasons are listed in [section 8](#8-if-something-looks-wrong).

Quick look at the most highly expressed TEs of the first sample:

```bash
sort -t $'\t' -k4,4gr 02_quant/matrices/TE_counts.tsv | head -20 | cut -f2,4-
```

## 7. Next: differential expression with DESeq2

This tutorial stops at the count tables. For differential expression we
recommend [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html),
reading kallisto's per-sample files with
[tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html).
tximport keeps kallisto's length information, which is more accurate than
rounding `gene_TE_counts.tsv`. Genes and TEs go into DESeq2 together, so they
are normalised the same way.

Make a small sample sheet, `samples.tsv` (TAB-separated; the names must
match the folders in `02_quant/`):

```
sample	condition
Ctrl_Rep1	Ctrl
Ctrl_Rep2	Ctrl
KO_Rep1	KO
KO_Rep2	KO
```

Then, in R (install once with
`BiocManager::install(c("tximport", "DESeq2"))`):

```r
library(tximport)
library(DESeq2)

# 1. which sample is which
samples <- read.delim("samples.tsv")
samples$condition <- factor(samples$condition, levels = c("Ctrl", "KO"))   # first = reference

# 2. kallisto results, one folder per sample
files <- file.path("02_quant", samples$sample, "abundance.tsv")
names(files) <- samples$sample

# 3. transcript -> gene table from step 1 (each TE maps to itself)
t2g <- read.delim("/path/to/references/mm10_TE/index/t2g.txt", header = FALSE, quote = "")
txi <- tximport(files, type = "kallisto", tx2gene = t2g[, 1:2],
                dropInfReps = TRUE)   # bootstraps are not needed for DESeq2

# 4. DESeq2: genes and TEs are normalised and tested together
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "KO", "Ctrl"))

# 5. add names and a gene / TE label, then save
id2name <- unique(t2g[, 2:3])
res$gene_name <- id2name$V3[match(rownames(res), id2name$V2)]
res$type <- ifelse(startsWith(res$gene_name, "TE:"), "TE", "gene")
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "DESeq2_KO_vs_Ctrl_genes_and_TEs.csv")
write.csv(as.data.frame(res[res$type == "TE", ]), "DESeq2_KO_vs_Ctrl_TEs_only.csv")
```

From here, follow the
[DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
for quality plots (PCA, sample distances), shrunken fold changes and
more complex designs. For a transcript-level analysis with the bootstraps,
see [sleuth](https://pachterlab.github.io/sleuth/).

A TE result is about a **whole subfamily**: a change in `MMERVK10C` means
reads from all its copies together went up or down, not one particular
locus.

## 8. If something looks wrong

| Message or symptom | Likely cause | What to do |
|--------------------|--------------|------------|
| `could not tell the strandedness` | the strand test gave neither a clear stranded nor a 50/50 answer | find the library kit and re-run with `--strand rf`, `fr` or `unstranded` (table in section 3) |
| `only N reads pseudoaligned - is this the right index` | reads and index don't match | check the species of the index and the reads |
| under 30% pseudoaligned | index from another species; heavy rRNA or genomic DNA contamination; total RNA with many intronic reads | check the species, look at the step 2 MultiQC report (GC content, duplication); for total RNA a lower rate is expected |
| very different % pseudoaligned between samples | a sample failed in the lab, or samples come from different kits | compare with the step 2 report; quantify different kits separately |
| `looks like droplet single-cell data` | 10x-style reads | use [step 4](05_singlecell_kb_count.md) |
| `index file not found` | wrong path | the file ends in `.idx` and sits in the `index/` folder of step 1 |
| kallisto log says the index is incompatible / cannot be read | index built with a different kallisto / kb-python version | use the `kallisto_te` environment for both steps (kb-python 0.30.2) |
| `has both paired-end and single-end files` | a folder mixes layouts for one sample | put the files of each sample in their own folder |
| two samples were combined that should not be | names differ only by the lane part (`_S1_L001`) or sit in the same per-sample folder | rename the files, or give each sample its own folder |
| `targets ... are not in t2g.txt` | results from another index are in the output folder | use a fresh `--outdir` for each index |

## 9. For your methods section

Adjust the versions (`kb --version`) and the settings you used:

> Transcript and TE subfamily abundances were estimated with kallisto
> v0.51.1 (via kb-python v0.30.2) against an index combining the GENCODE
> M25 transcriptome with RepBase 25.05 TE consensus sequences, with the genome
> as D-list (de Villarreal *et al.*, 2023). Reads were quantified as
> reverse-stranded (`--rf-stranded`) with 30 bootstraps; single-end samples
> used an estimated fragment length of 200 ± 20 bp. Transcript-level estimates
> were summarised to genes and TE subfamilies with tximport, and differential
> expression was tested with DESeq2.

Next: single-cell data in [step 4](05_singlecell_kb_count.md), or go
straight to DESeq2 (section 7).
