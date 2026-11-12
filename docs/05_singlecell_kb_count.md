# Step 4 — Count genes and TEs per cell (single-cell RNA-seq, kb count)

This step turns droplet single-cell RNA-seq reads (10x Genomics, Drop-seq,
...) into a **cell × (gene + TE) count matrix** that you can load into
Seurat or scanpy. It uses the same index as bulk RNA-seq (step 1). No prior
experience is assumed. One command does everything:

```bash
bash scripts/quant/run_kb_count.sh \
    --input  /path/to/project/00_fastq_10x \
    --index  /path/to/references/mm10_TE/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx \
    --tech   10xv3 \
    --outdir /path/to/project/03_kb_count \
    --threads 8
```

> **Use the raw FASTQ files. Do not trim single-cell reads** (skip step 2).
> Read 1 holds the cell barcode and UMI at fixed positions; trimming shifts
> or shortens them. The script refuses files from step 2 (`*.trimmed.*`) and
> read 1 files too short for the chosen technology.

This tutorial stops at the count matrix. For clustering, cell types and
differential expression between cells, [section 6](#6-next-seurat) shows how
to load the matrix into Seurat and links to the standard tutorials.

---

## Contents

1. [What this step does](#1-what-this-step-does)
2. [Before you start: files and technology](#2-before-you-start-files-and-technology)
3. [Run it](#3-run-it)
4. [What you get](#4-what-you-get)
5. [Read the summary](#5-read-the-summary)
6. [Next: Seurat](#6-next-seurat)
7. [If something looks wrong](#7-if-something-looks-wrong)
8. [For your methods section](#8-for-your-methods-section)

---

## 1. What this step does

In droplet single-cell RNA-seq every read pair carries:

| Read | Contains | Example (10x v3) |
|------|----------|------------------|
| Read 1 (`R1`) | **cell barcode** (which droplet / cell) + **UMI** (which RNA molecule) | 16 bp barcode + 12 bp UMI = 28 bp |
| Read 2 (`R2`) | the **cDNA**: the part that tells you the gene (or TE) | usually 90 bp |
| Index 1 / 2 (`I1`, `I2`) | sample index used for demultiplexing | not needed here; ignored |

[`kb count`](https://www.kallistobus.tools/) runs kallisto to pseudoalign
read 2 against the genes + TE index, and bustools to correct the cell
barcodes against the technology's barcode list and count each UMI once. Each
TE subfamily is a "gene" in the result, so TEs appear as extra rows (named
`TE:...`) next to the genes.

One difference from bulk: by default `kb count` **drops** a read whose UMI
fits several genes at once, e.g. a read shared by two related TE
subfamilies, or by a TE and a gene that contains a copy of it. Bulk kallisto
splits such reads instead. The script keeps the kb default, so single-cell TE
counts are conservative.

## 2. Before you start: files and technology

### Your files

Put the raw FASTQ files of one experiment in one folder, as they come from
the sequencing facility or `cellranger mkfastq` / `bcl2fastq`:

```
00_fastq_10x/
├── Mouse1_S1_L001_I1_001.fastq.gz    # index read: ignored
├── Mouse1_S1_L001_R1_001.fastq.gz    # barcode + UMI
├── Mouse1_S1_L001_R2_001.fastq.gz    # cDNA
├── Mouse1_S1_L002_R1_001.fastq.gz    # a second lane of the same sample
├── Mouse1_S1_L002_R2_001.fastq.gz
├── Mouse2_S2_L001_R1_001.fastq.gz
└── Mouse2_S2_L001_R2_001.fastq.gz
```

Everything before `_S1_L001` is the sample name, and **all lanes of a sample
are combined**. A folder per sample (`00_fastq_10x/Mouse1/*.fastq.gz`) also
works: then the folder name is the sample name.

### Which technology (`--tech`)?

| Kit | `--tech` | Read 1 length |
|-----|----------|---------------|
| 10x Chromium Single Cell 3′ **v2** | `10xv2` | 26 bp |
| 10x Chromium Single Cell 3′ **v3 / v3.1** | `10xv3` | 28 bp |
| 10x Chromium **GEM-X** Single Cell 3′ **v4** | `10xv4` | 28 bp |
| Drop-seq | `dropseq` | 20 bp |
| others (inDrops, BD Rhapsody, SPLiT-seq, ...) | see `kb --list` | |

The kit name is on the sequencing facility's report, the GEO record, or the
Cell Ranger web summary ("Chemistry"). To see the read 1 length yourself:

```bash
gzip -dc Mouse1_S1_L001_R1_001.fastq.gz | head -2 | tail -1 | awk '{ print length }'
```

A read 1 longer than listed (e.g. 150 bp from a paired 150 bp run) is fine:
kb only uses the first 26 / 28 bases. **10x v3 and v4 have the same read
layout but different barcode lists**, so read length cannot tell them apart.
A wrong choice shows up as few reads with a valid barcode (section 5).

## 3. Run it

**Interactively** (needs about as much memory as the index plus a few GB):

```bash
conda activate kallisto_te
cd /path/to/kallisto-TE

bash scripts/quant/run_kb_count.sh \
    --input   /path/to/project/00_fastq_10x \
    --index   /path/to/references/mm10_TE/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx \
    --tech    10xv3 \
    --outdir  /path/to/project/03_kb_count \
    --threads 8 \
    --memory  16G
```

**On a SLURM cluster**: edit the paths and `TECH` in
[`slurm/run_kb_count.sbatch`](../slurm/run_kb_count.sbatch), then

```bash
mkdir -p logs
sbatch slurm/run_kb_count.sbatch
```

If a job stops halfway, submit it again: samples that already have a matrix
are skipped (use `--force` to redo them).

Options:

| Option | What it does |
|--------|--------------|
| `--memory 16G` | memory bustools may use for sorting barcodes (keep it below the job's memory) |
| `--filter` | also writes `counts_filtered/`: only barcodes that look like cells (bustools knee-plot filter). You can instead filter in Seurat, as in section 6 |
| `--h5ad` | also writes an AnnData `.h5ad` file, for [scanpy](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html) users |

## 4. What you get

```
03_kb_count/
├── Mouse1/
│   ├── counts_unfiltered/
│   │   ├── cells_x_genes.mtx               # ← UMI counts: rows = barcodes, columns = genes/TEs
│   │   ├── cells_x_genes.barcodes.txt      #   one barcode per row of the matrix
│   │   ├── cells_x_genes.genes.txt         #   gene / TE IDs, one per column
│   │   └── cells_x_genes.genes.names.txt   #   gene names; TEs start with "TE:"
│   ├── counts_filtered/                    #   (only with --filter) the same, cells only
│   ├── run_info.json                       #   reads processed / pseudoaligned
│   ├── inspect.json                        #   barcode statistics
│   └── output.bus, matrix.ec, ...          #   intermediate files kb needs
├── Mouse2/ ...
├── kb_count_summary.tsv                    # ← one line per sample (also printed at the end)
└── logs/<sample>.log
```

`counts_unfiltered` holds **every barcode seen**, including empty droplets
with a handful of UMIs. That is normal: cells are picked in the next step.

## 5. Read the summary

`kb_count_summary.tsv` has one line per sample:

| Column | Meaning | What to look for |
|--------|---------|------------------|
| `reads_processed` | read pairs kb saw | close to the read count from the facility |
| `percent_pseudoaligned` | % of reads whose cDNA matched a gene or TE | similar across samples; low values point to the wrong species or poor libraries |
| `percent_valid_barcode` | % of reads with a barcode on the technology's list | usually well above 90% with the right `--tech`; the script warns below 70% |
| `barcodes` | all barcodes with at least one UMI | many more than cells: mostly empty droplets |
| `barcodes_500plus_UMIs` | barcodes with ≥ 500 UMIs: a rough cell count | in the range of the number of cells you expected to capture |
| `median_UMIs_500plus` | median UMIs of those barcodes | depends on sequencing depth and cell type |
| `percent_UMIs_TE` | % of all UMIs counted on TE subfamilies | compare samples with each other |

## 6. Next: Seurat

[Seurat](https://satijalab.org/seurat/articles/install) is the most widely
used R package for single-cell analysis. Load one sample like this (it
works with `counts_filtered` too):

```r
library(Seurat)
library(Matrix)

d <- "03_kb_count/Mouse1/counts_unfiltered"
counts <- readMM(file.path(d, "cells_x_genes.mtx"))         # kb writes cells x genes ...
counts <- t(as(counts, "CsparseMatrix"))                   # ... Seurat wants genes x cells
colnames(counts) <- readLines(file.path(d, "cells_x_genes.barcodes.txt"))
rownames(counts) <- make.unique(readLines(file.path(d, "cells_x_genes.genes.names.txt")))

# keep barcodes with at least 200 detected genes / TEs: drops empty droplets
seu <- CreateSeuratObject(counts, project = "Mouse1", min.cells = 3, min.features = 200)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")   # "^MT-" for human
seu[["percent.TE"]] <- PercentageFeatureSet(seu, pattern = "^TE:")   # share of UMIs from TEs
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.TE"), ncol = 4)
```

Things to know:

* **TE names change slightly.** Seurat does not allow `_` in feature names
  and replaces it with `-` (it prints a warning), so
  `TE:MMERVK10C_ERV2_Mus_musculus` becomes `TE:MMERVK10C-ERV2-Mus-musculus`.
  The `TE:` prefix stays, so `grep("^TE:", rownames(seu), value = TRUE)`
  lists all TEs.
* **TEs are ordinary features.** They are normalised, scaled and used for
  clustering like genes. Plot one with
  `FeaturePlot(seu, "TE:MMERVK10C-ERV2-Mus-musculus")`, or compare clusters with
  `FindMarkers()`.
* **Several samples**: load each one as above, then combine them with
  `merge()` or integrate them.

Continue with these tutorials:

* [Seurat: guided clustering tutorial (PBMC 3k)](https://satijalab.org/seurat/articles/pbmc3k_tutorial):
  quality control, normalisation, clustering, marker genes. Start here.
* [Seurat: integrating several samples](https://satijalab.org/seurat/articles/integration_introduction)
* [Single-cell best practices](https://www.sc-best-practices.org/): a free
  online book covering the whole analysis (Python / scanpy based).
* [Orchestrating Single-Cell Analysis with Bioconductor (OSCA)](https://bioconductor.org/books/release/OSCA/),
  including empty-droplet detection with
  [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html).
* Python users: run with `--h5ad` and follow the
  [scanpy clustering tutorial](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html).

## 7. If something looks wrong

| Message or symptom | Likely cause | What to do |
|--------------------|--------------|------------|
| `uses trimmed files (*.trimmed.*)` | reads from step 2 | use the raw FASTQ files |
| `read 1 is only N bp, but ... needs at least M bp` | wrong `--tech`, or trimmed reads | check the kit (section 2) and use the raw files |
| `Read 1 and read 2 look swapped` | file names swapped, e.g. by a download tool | read 1 must be the short barcode + UMI read; rename the files |
| `has no read 2 file` | a lane is missing its R2, or the files are named `_1` / `_2` / `_3` (some SRA downloads) | make sure each lane has `_R1` and `_R2` files; for SRA, check which file is the barcode read and rename them |
| `under 70% of reads carry a valid cell barcode` | wrong chemistry, most often 10x v3 vs v4 | re-run with the other `--tech` and `--force` |
| `unknown technology` | name not in `kb --list` | use one of the names in section 2 |
| very few barcodes with ≥ 500 UMIs | low-quality run, or wrong `--tech` | check `percent_valid_barcode` first |
| job killed / out of memory | job memory too low | raise the job's memory (`--mem` in the sbatch file); keep `--memory` well below it |

## 8. For your methods section

Adjust the versions (`kb --version`) and the technology:

> Single-cell reads were pseudoaligned with kb-python v0.30.2 (kallisto
> v0.51.1, bustools v0.45.1) using the 10x Genomics v3 technology setting
> against an index combining the GENCODE M25 transcriptome with RepBase 25.05
> TE consensus sequences, with the genome as D-list (de Villarreal *et al.*,
> 2023). Gene and TE subfamily UMI count matrices were analysed with Seurat.
