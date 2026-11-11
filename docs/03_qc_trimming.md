# Step 2 — Trim and check your reads (bulk RNA-seq)

Before quantifying, remove adapter sequence and low-quality bases from the
reads and check that each sample looks healthy. This page assumes no prior
experience. One command does everything:

```bash
bash scripts/qc/run_fastp.sh --input /path/to/fastq_folder --outdir /path/to/01_trimmed --threads 8
```

> **Single-cell RNA-seq (10x Genomics, Drop-seq, ...): skip this step.**
> In droplet data, read 1 holds the cell barcode and UMI. Trimming it
> breaks the barcodes. `kb count` takes the raw FASTQ files directly.
> As a safety net, the script recognises this layout (a short read 1 of
> ≤ 30 bp with a read 2 at least twice as long) and stops before trimming
> anything.

---

## Contents

1. [What this step does](#1-what-this-step-does)
2. [Look at your FASTQ files first](#2-look-at-your-fastq-files-first)
3. [Run fastp](#3-run-fastp)
4. [What you get](#4-what-you-get)
5. [Read the results](#5-read-the-results)
6. [If something looks wrong](#6-if-something-looks-wrong)
7. [Optional extras](#7-optional-extras)
8. [For your methods section](#8-for-your-methods-section)

---

## 1. What this step does

[fastp](https://github.com/OpenGene/fastp) reads each FASTQ file once and,
with its default settings:

| What | Why |
|------|-----|
| **Cuts adapter sequence** from the 3' end of reads | When the RNA fragment is shorter than the read, the sequencer reads into the adapter. Adapter bases match nothing in the index |
| **Removes low-quality reads** (more than 40% of bases below quality 15, or more than 5 `N`s) | Sequencing errors make reads harder to assign |
| **Removes very short reads** (under 15 bp after trimming) | Too short to place reliably |
| **Trims poly-G tails** on NovaSeq / NextSeq data (switched on automatically) | Two-colour sequencers report "G" when there is no signal |
| **Writes a QC report** per sample | Read quality, GC content, duplication, adapter content, insert size |

[MultiQC](https://multiqc.info/) then puts every sample's report into one
web page, so you can compare samples side by side.

Both tools are in the `kallisto_te` environment
(`mamba env create -f envs/kallisto_te.yaml`).

---

## 2. Look at your FASTQ files first

Put all raw FASTQ files of one experiment in one folder, then look at the
names:

```bash
ls /path/to/fastq_folder
```

The script recognises these naming styles automatically:

| Layout | Example file names | Sample name it uses |
|--------|--------------------|---------------------|
| Paired-end | `Ctrl_Rep1_R1.fq.gz` + `Ctrl_Rep1_R2.fq.gz` | `Ctrl_Rep1` |
| Paired-end (Illumina) | `Ctrl_S1_L001_R1_001.fastq.gz` + `Ctrl_S1_L001_R2_001.fastq.gz` | `Ctrl_S1_L001` |
| Paired-end (SRA / ENA downloads) | `SRR1234567_1.fastq.gz` + `SRR1234567_2.fastq.gz` | `SRR1234567` |
| Single-end | `Ctrl_Rep1.fq.gz` (or `Ctrl_Rep1_R1.fq.gz` with no R2) | `Ctrl_Rep1` |

> **Paired-end** means two files per sample, read 1 and read 2, from the two
> ends of the same RNA fragment. **Single-end** means one file per sample.
> Your sequencing facility's report tells you which you have (e.g. "PE150"
> = paired-end, 150 bp).

Peek inside one file to see what a read looks like. Each read takes 4
lines: name, sequence, `+`, quality.

```bash
gzip -dc /path/to/fastq_folder/Ctrl_Rep1_R1.fq.gz | head -8
```

```
@A00123:45:HXXXXDSXY:1:1101:1234:1000 1:N:0:ACGTACGT
NTGCAGGTCAGCTTGACCATGGAGCTGCAGGCCTTCAAGGAGCTGGAGGCGGAGCTG...
+
#FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF...
```

If your samples were sequenced on several lanes (`_L001_`, `_L002_`, ...),
merge the lanes of each sample first. Plain `cat` works on gzipped files:

```bash
cat Ctrl_S1_L001_R1_001.fastq.gz Ctrl_S1_L002_R1_001.fastq.gz > Ctrl_R1.fastq.gz
cat Ctrl_S1_L001_R2_001.fastq.gz Ctrl_S1_L002_R2_001.fastq.gz > Ctrl_R2.fastq.gz
```

---

## 3. Run fastp

### On your computer or an interactive cluster session

```bash
conda activate kallisto_te
cd /path/to/kallisto-TE

bash scripts/qc/run_fastp.sh \
    --input   /path/to/project/00_fastq \
    --outdir  /path/to/project/01_trimmed \
    --threads 8
```

The script first lists the samples it found. Check that every sample is
there and that pairs are matched correctly:

```
Found 4 sample(s) in /path/to/project/00_fastq: 4 paired-end, 0 single-end (read-1 tag: _R1)
  Ctrl_Rep1                      Ctrl_Rep1_R1.fq.gz  +  Ctrl_Rep1_R2.fq.gz
  Ctrl_Rep2                      Ctrl_Rep2_R1.fq.gz  +  Ctrl_Rep2_R2.fq.gz
  KD_Rep1                        KD_Rep1_R1.fq.gz  +  KD_Rep1_R2.fq.gz
  KD_Rep2                        KD_Rep2_R1.fq.gz  +  KD_Rep2_R2.fq.gz
  [fastp] Ctrl_Rep1
  ...
```

### As a SLURM job (recommended for many samples)

Edit the `EDIT` lines in [`slurm/run_fastp.sbatch`](../slurm/run_fastp.sbatch)
(the two folders and your partition), then:

```bash
mkdir -p logs
sbatch slurm/run_fastp.sbatch
squeue -u $USER                    # is it running?
tail -f logs/fastp_trim_*.out      # follow its progress (Ctrl+C to stop watching)
```

If a job stops halfway, submit it again. Samples that are already trimmed
are skipped (use `--force` to redo them).

---

## 4. What you get

```
01_trimmed/
├── trimmed/                         # ← use these files for quantification (step 3)
│   ├── Ctrl_Rep1_R1.trimmed.fq.gz
│   ├── Ctrl_Rep1_R2.trimmed.fq.gz
│   └── ...
├── reports/
│   ├── Ctrl_Rep1.fastp.html         # one QC report per sample
│   ├── Ctrl_Rep1.fastp.json         # same, machine-readable (used by MultiQC)
│   └── Ctrl_Rep1.fastp.log
├── multiqc/
│   └── multiqc_report.html          # ← all samples in one page: open this
└── fastp_summary.tsv                # one line per sample (also printed at the end)
```

---

## 5. Read the results

### The summary table

Printed at the end of the run and saved as `fastp_summary.tsv`:

```
sample     layout  reads_before  reads_after  percent_kept  q30_after  duplication  adapter_trimmed_reads
Ctrl_Rep1  paired  61234880      60512344     98.8          93.1       35.2         4123456
KD_Rep1    paired  58877102      58190016     98.8          92.8       33.9         3998721
```

| Column | Meaning | Usually fine (rough guide) |
|--------|---------|----------------------------|
| `reads_before` / `reads_after` | reads in / out (both mates counted for paired-end) | similar across samples of one experiment |
| `percent_kept` | share of reads that passed the filters | > 90% |
| `q30_after` | % of bases with quality ≥ 30 (1 error in 1,000) after trimming | > 85% on current Illumina machines |
| `duplication` | % of duplicated reads | 20–60% is **normal for RNA-seq**: highly expressed genes give many identical reads. Do not remove duplicates |
| `adapter_trimmed_reads` | reads that had adapter cut off | anything; depends on fragment size. Higher in degraded or short-insert libraries |

The most useful check is whether **one sample looks very different from
the others**. That points to a problem with that library, not with the
analysis.

### The MultiQC report

`multiqc/multiqc_report.html` opens in any web browser. On a cluster, first
copy it to your own computer (FileZilla, Cyberduck, or `scp`):

```bash
# run this on YOUR computer, not on the cluster
scp username@cluster.address:/path/to/project/01_trimmed/multiqc/multiqc_report.html .
```

Worth looking at:

* **General Statistics** (top table): one row per sample, the same numbers
  as above. Sort by a column to spot outliers.
* **Filtered Reads**: how many reads were removed and why (low quality, too
  many `N`, too short). Most reads should be "Passed Filter".
* **Sequence Quality**: quality along the read should stay high (green).
  A drop at the very end is normal.
* **Insert Sizes** (paired-end only): the RNA fragment lengths. A peak far
  below the read length means many short fragments, which is common for
  degraded RNA.
* **GC Content**: samples should overlap. A second peak can mean rRNA or
  contamination.

---

## 6. If something looks wrong

| You see | Likely cause | What to do |
|---------|--------------|------------|
| `no FASTQ files ... in` | wrong folder, or files are not `.fq.gz` / `.fastq.gz` / `.fq` / `.fastq` | check the path with `ls` |
| `These samples look like droplet single-cell RNA-seq` ... `stopped before trimming` | read 1 is ≤ 30 bp and read 2 much longer: the barcode + cDNA layout of 10x / Drop-seq data | single-cell data: don't trim, use the raw files with `kb count`. If they really are bulk reads with a short read 1, re-run with `--force-sc` |
| Pairs listed as two single-end samples | read 1 / read 2 are named differently than expected | tell the script the read-1 tag, e.g. `--r1-tag _read1` (read 2 must then be `_read2`) |
| `percent_kept` far below 90% in one sample | low-quality or failed library | look at that sample's "Filtered Reads" and quality plots; ask the sequencing core |
| `adapter_trimmed_reads` is 0 for **single-end** data | fastp could not detect the adapter automatically | give it explicitly: `--fastp-args "--adapter_sequence AGATCGGAAGAGC"` (standard Illumina TruSeq adapter start) |
| `fastp not found` | environment not active | `conda activate kallisto_te` |

---

## 7. Optional extras

You do **not** need these for kallisto-TE. fastp's report already covers
the usual quality checks.

* **FastQC on the raw reads**, if you want the classic FastQC plots as well:

  ```bash
  mamba install -n kallisto_te -c bioconda fastqc
  bash scripts/qc/run_fastp.sh --input ... --outdir ... --fastqc
  ```

  The FastQC results go to `fastqc_raw/` and are included in the MultiQC
  report.
* **Extra fastp settings** are passed straight through with
  `--fastp-args`, e.g. `--fastp-args "--detect_adapter_for_pe"` (also
  search for adapter sequences in paired-end data) or
  `--fastp-args "--length_required 25"`. Use the same settings for every
  sample in a project.
* **Contamination screening** (FastQ Screen, Kraken2) is only worth doing
  if MultiQC shows something odd, e.g. an unexpected GC peak.

---

## 8. For your methods section

`run_fastp.sh` uses fastp's default settings. A sentence like this covers
it (use the versions from `fastp --version` and `multiqc --version`):

> Adapters and low-quality bases were removed with fastp v1.4.0 using
> default parameters, and quality reports were summarised with MultiQC
> v1.35.

Next: quantification of genes and TEs with the trimmed reads (step 3,
coming soon), using the files in `01_trimmed/trimmed/`.
