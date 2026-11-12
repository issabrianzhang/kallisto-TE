# Step 5 — Run a project with one command (Snakemake)

Once you have a genes + TE index (step 1), steps 2–4 can run one after
another from one short settings file, with
[Snakemake](https://snakemake.readthedocs.io/). Each step runs the same
script as in the step-by-step guides, so the results are identical. Snakemake
adds three things:

* **One command** per project: trimming → bulk counts and/or single-cell
  counts. Bulk, single-cell or both: you choose by which reads you give it.
* **Picks up where it stopped.** If a run is interrupted (power cut, cluster
  time limit), run the same command again: finished samples are kept.
* **Keeps results consistent.** Change a setting and the affected steps are
  redone for every sample, so old and new settings never mix.

New to this? Do the step-by-step guides once ([steps 1–4](../README.md#tutorial));
then this page is all you need for the next project.

---

## Contents

1. [Build the index first (once per species)](#1-build-the-index-first-once-per-species)
2. [Install and test](#2-install-and-test)
3. [Write your settings file](#3-write-your-settings-file)
4. [Run it on your computer](#4-run-it-on-your-computer)
5. [Run it on a SLURM cluster](#5-run-it-on-a-slurm-cluster)
6. [What you get](#6-what-you-get)
7. [Re-running, adding samples, changing settings](#7-re-running-adding-samples-changing-settings)
8. [If something goes wrong](#8-if-something-goes-wrong)

---

## 1. Build the index first (once per species)

The pipeline does **not** build an index; it needs one that already exists.
Build it once per species and genome / annotation / TE library release, with
[step 1](01_build_index.md) (or [step 1b](02_custom_species.md) for other
species), and keep it in a shared place, for example:

```
references/
├── mm10_GencodeM25_RepBase25.05/     # build_index.sh --outdir
│   ├── build_info.txt                # how it was built: inputs, versions, checks
│   └── index/
│       ├── mm10_GencodeM25_RepBase25.05_genes_TE.idx
│       └── t2g.txt
└── hg38_GencodeV40_RepBase28.04/
```

Then every project of that species uses the same index:

* **Comparable**: all your analyses count against exactly the same
  reference.
* **Reproducible**: `build_info.txt` records the genome, annotation, TE
  library and software versions for your methods section.
* **Fast**: building takes hours and about 32 GB of memory for a mammalian
  genome; each project then starts straight away.

Build and quantify with the same kb-python version (the `kallisto_te`
environment): index formats change between kallisto releases.

## 2. Install and test

Snakemake is part of the `kallisto_te` environment (`envs/kallisto_te.yaml`).
If you created the environment before step 5 existed, add it:

```bash
conda activate kallisto_te
mamba install -c conda-forge -c bioconda snakemake-minimal=9.27.1 snakemake-executor-plugin-slurm=2.8.0
```

Then check that everything works (about 1 minute; it builds a toy index from
a fake genome and runs all steps on simulated reads):

```bash
cd /path/to/kallisto-TE
bash test/run_pipeline_test.sh
```

It should end with `PASS: the pipeline works`.

## 3. Write your settings file

All settings live in one YAML file. [`config/config.yaml`](../config/config.yaml)
lists every option with a comment. **Don't edit that file**: make a short
file for your project with only what differs. Everything you leave out keeps
its default.

**Bulk RNA-seq only:**

```yaml
outdir: /path/to/project/results
index: /path/to/references/mm10_GencodeM25_RepBase25.05/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx

bulk:
  fastq_dir: /path/to/project/00_fastq
```

**Single-cell RNA-seq only:**

```yaml
outdir: /path/to/project/results
index: /path/to/references/mm10_GencodeM25_RepBase25.05/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx

singlecell:
  fastq_dir: /path/to/project/00_fastq_10x
  tech: 10xv3
```

**Both, with a few settings changed:**

```yaml
outdir: /path/to/project/results
index: /path/to/references/mm10_GencodeM25_RepBase25.05/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx

bulk:
  fastq_dir: /path/to/project/00_fastq_bulk
  fastp_args: "--adapter_sequence AGATCGGAAGAGC"   # if step 2 warned "no adapter"
  bootstraps: 0                                    # only using DESeq2

singlecell:
  fastq_dir: /path/to/project/00_fastq_10x
  tech: 10xv4
```

The main options:

| Section | Option | Meaning |
|---------|--------|---------|
| (top) | `outdir` | where all results go |
| | `index` | the `.idx` file from step 1 (required) |
| | `t2g` | only if `t2g.txt` is not next to the index |
| `bulk` | `fastq_dir` | folder of raw bulk FASTQ files; empty = no bulk |
| | `trim` | `true` = run fastp first (step 2) |
| | `fastp_args`, `single_end` | as in [step 2](03_qc_trimming.md) |
| | `strand`, `fragment_length`, `fragment_sd`, `bootstraps` | as in [step 3](04_bulk_kallisto_quant.md#3-three-settings-to-know) |
| `singlecell` | `fastq_dir`, `tech` | raw single-cell FASTQ folder and technology; empty = no single-cell |
| | `memory`, `filter`, `h5ad` | as in [step 4](05_singlecell_kb_count.md#3-run-it) |
| `resources` | `threads`, `mem_mb`, `runtime` | per step; used on a cluster |

Paths can be absolute (safest) or relative to the folder you run
`snakemake` from.

## 4. Run it on your computer

From the repository folder, first a **dry run**: it lists what would run
and checks the settings file, without running anything:

```bash
conda activate kallisto_te
cd /path/to/kallisto-TE
snakemake -n --configfile /path/to/project/my_project.yaml
```

Then the real run (`--cores` = how many CPU threads to use):

```bash
snakemake --cores 8 --configfile /path/to/project/my_project.yaml
```

A computer with 8 cores and 16 GB of memory is enough for steps 2–4 with a
mouse or human index (counting used under 2 GB in our tests); for
single-cell data set `singlecell: memory: 4G`. It just takes longer than a
cluster.

Snakemake keeps its own bookkeeping in a `.snakemake/` folder where you run
it. To keep that next to your project instead of in the repository, run from
the project folder and point at the pipeline file:

```bash
cd /path/to/project
snakemake --cores 8 --configfile my_project.yaml --snakefile /path/to/kallisto-TE/workflow/Snakefile
```

## 5. Run it on a SLURM cluster

[`profiles/slurm/config.yaml`](../profiles/slurm/config.yaml) makes every
step a SLURM job. Edit two things once:

* `slurm_partition`: your cluster's partition (`sinfo` lists them)
* `slurm_account`: only if your cluster requires one

Each step's threads, memory and time limit come from the `resources:`
section of the settings file (defaults in `config/config.yaml`). Each step is
**one job for all samples**, so for a large project raise its `runtime`
(minutes): for example bulk quantification takes roughly 20 minutes per
sample of 30 million reads. Then:

```bash
conda activate kallisto_te
cd /path/to/kallisto-TE
snakemake -n --profile profiles/slurm --configfile /path/to/project/my_project.yaml   # dry run
snakemake --profile profiles/slurm --configfile /path/to/project/my_project.yaml
```

Snakemake itself keeps running while the jobs are queued and running, so
start it where it survives you logging out: inside `tmux` or `screen` on the
login node, or as a small SLURM job of its own (1 CPU, 2 GB, the full
runtime), if your cluster prefers that.

## 6. What you get

```
results/
├── 01_trimmed/           # step 2: trimmed reads, fastp reports, MultiQC
├── 02_quant/             # step 3: <sample>/abundance.tsv, quant_summary.tsv
├── 03_kb_count/          # step 4: <sample>/counts_unfiltered/, kb_count_summary.tsv
└── logs/
    ├── 2_trim_bulk.log
    ├── 3_quant_bulk.log
    └── 4_kb_count.log
```

Each folder is exactly what the matching step-by-step guide describes, so
the downstream code there works unchanged:
[DESeq2 for bulk](04_bulk_kallisto_quant.md#7-next-differential-expression-with-deseq2),
[Seurat for single-cell](05_singlecell_kb_count.md#6-next-seurat).

Each step folder also has a `pipeline_settings.json`: the settings (and
index) that step last ran with. Snakemake uses it to notice changes
(section 7), and it documents how the results were made.

## 7. Re-running, adding samples, changing settings

Run the same command again in every case; Snakemake works out what to do.

| What happened | What the next run does |
|---------------|------------------------|
| nothing changed | nothing ("Nothing to be done") |
| the run was interrupted | finishes the samples that were not done yet |
| you added FASTQ files for new samples | trims and quantifies only the new samples |
| you changed a setting (e.g. `strand`, `fastp_args`, `tech`) | redoes that step, and the steps after it, for **every** sample |
| you pointed `index` at another index, or rebuilt it | redoes all quantification and counting (trimming is kept) |

To redo a step from scratch anyway, delete its folder (e.g.
`results/02_quant/`) and run again.

## 8. If something goes wrong

| Message or symptom | What to do |
|--------------------|------------|
| `kallisto-TE config: ...` | the settings file has a problem; the message says which option. Nothing has run yet |
| `set 'index' to the .idx file ...` / `index file not found` | build the index first (section 1), then give the path of its `.idx` file |
| `'index' must be the path of an .idx file` | `index:` is a single path, not a section with `species:` etc.; the pipeline does not build indexes |
| `nothing to do: set bulk.fastq_dir ... and/or singlecell.fastq_dir` | give at least one folder of reads |
| `Error in rule quant_bulk` (or another step) | open that step's log in `results/logs/`. It is the same output as running the script by hand, so the troubleshooting table of that step's guide applies |
| `bulk reads are in one sub-folder per sample, but trimming needs ...` | fastp needs all files in one folder: move them, or set `bulk.trim: false` to quantify the raw reads |
| `snakemake: command not found` | `conda activate kallisto_te`, or install Snakemake (section 2) |
| `LockException` / "Directory cannot be locked" | an earlier run was killed. Make sure no other run is using the folder, then `snakemake --unlock --configfile ...` |
| `IncompleteFilesException` | an earlier run stopped mid-step. Add `--rerun-incomplete`; finished samples are still kept |
| SLURM jobs fail at once | check `slurm_partition` / `slurm_account` in `profiles/slurm/config.yaml`; SLURM's own messages are under `.snakemake/slurm_logs/` |
| a job is killed for time or memory | raise `runtime` (minutes) or `mem_mb` for that step under `resources:` and run again; finished samples are kept |

For your methods section, cite Snakemake (Mölder *et al.*, 2021) along with
the tools of each step, and the index details from `build_info.txt`.
