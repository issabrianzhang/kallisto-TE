# Step 5 — Run everything with one command (Snakemake)

Steps 1–4 can run one after another from one settings file, with
[Snakemake](https://snakemake.readthedocs.io/). Each step runs the same
script as in the step-by-step guides, so the results are identical. Snakemake
adds three things:

* **One command** for the whole project: index → trimming → bulk counts →
  single-cell counts. Steps you don't need are skipped.
* **Picks up where it stopped.** If a run is interrupted (power cut, cluster
  time limit), run the same command again: finished samples are kept.
* **Keeps results consistent.** Change a setting and the affected steps are
  redone for every sample, so old and new settings never mix.

New to this? Do the step-by-step guides once ([steps 1–4](../README.md#tutorial));
then this page is all you need for the next project.

---

## Contents

1. [Install and test](#1-install-and-test)
2. [Write your settings file](#2-write-your-settings-file)
3. [Run it on your computer](#3-run-it-on-your-computer)
4. [Run it on a SLURM cluster](#4-run-it-on-a-slurm-cluster)
5. [What you get](#5-what-you-get)
6. [Re-running, adding samples, changing settings](#6-re-running-adding-samples-changing-settings)
7. [If something goes wrong](#7-if-something-goes-wrong)

---

## 1. Install and test

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

## 2. Write your settings file

All settings live in one YAML file. [`config/config.yaml`](../config/config.yaml)
lists every option with a comment. **Don't edit that file**: make a short
file for your project with only what differs. Everything you leave out keeps
its default.

**Mouse bulk RNA-seq, building a new index from the bundled Dfam library:**

```yaml
outdir: /path/to/project/results

index:
  species: mouse
  name: mm10_GencodeM25_Dfam4.0
  dfam: true                       # genome + GTF are downloaded (GENCODE M25)

bulk:
  fastq_dir: /path/to/project/00_fastq
```

**Re-using an index you already built, bulk + single-cell:**

```yaml
outdir: /path/to/project/results

index:
  existing_index: /path/to/references/mm10_TE/index/mm10_GencodeM25_RepBase25.05_genes_TE.idx

bulk:
  fastq_dir: /path/to/project/00_fastq_bulk
  fastp_args: "--adapter_sequence AGATCGGAAGAGC"   # if step 2 warned "no adapter"

singlecell:
  fastq_dir: /path/to/project/00_fastq_10x
  tech: 10xv3
```

**Another species** (`species: custom`): give the genome, the GTF and a TE
library, as in [docs/02_custom_species.md](02_custom_species.md):

```yaml
index:
  species: custom
  name: GRCz11_Dfam4.0
  dfam: false
  te_fasta: [te_libraries/dfam/zebrafish_Dfam4.0_curated.fa.gz]
  genome: /path/to/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
  gtf: /path/to/Danio_rerio.GRCz11.113.gtf.gz
```

The main options:

| Section | Option | Meaning |
|---------|--------|---------|
| (top) | `outdir` | where all results go |
| `index` | `existing_index` | an index from an earlier build; the rest of `index` is then ignored |
| | `species`, `name` | `mouse`, `human` or `custom`; a name for the index file |
| | `dfam` / `repbase` / `te_fasta` | the TE library: bundled Dfam (mouse / human), a RepBase archive, your own FASTA files |
| | `genome`, `gtf` | required for `custom`; for mouse / human empty = download the preset |
| `bulk` | `fastq_dir` | folder of raw bulk FASTQ files; empty = no bulk |
| | `trim` | `true` = run fastp first (step 2) |
| | `strand`, `fragment_length`, `fragment_sd`, `bootstraps` | as in [step 3](04_bulk_kallisto_quant.md#3-three-settings-to-know) |
| `singlecell` | `fastq_dir`, `tech` | raw single-cell FASTQ folder and technology; empty = no single-cell |
| `resources` | `threads`, `mem_mb`, `runtime` | per step; used on a cluster |

Paths can be absolute (safest) or relative to the folder you run
`snakemake` from.

## 3. Run it on your computer

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

Snakemake keeps its own bookkeeping in a `.snakemake/` folder where you run
it. To keep that next to your project instead of in the repository, run from
the project folder and point at the pipeline file:

```bash
cd /path/to/project
snakemake --cores 8 --configfile my_project.yaml --snakefile /path/to/kallisto-TE/workflow/Snakefile
```

## 4. Run it on a SLURM cluster

[`profiles/slurm/config.yaml`](../profiles/slurm/config.yaml) makes every
step a SLURM job. Edit two things once:

* `slurm_partition`: your cluster's partition (`sinfo` lists them)
* `slurm_account`: only if your cluster requires one

Each step's threads, memory and time limit come from the `resources:`
section of the settings file (defaults in `config/config.yaml`). Then:

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

## 5. What you get

```
results/
├── reference/            # step 1 (only when the index is built here)
│   └── index/<name>_genes_TE.idx, t2g.txt, ...
├── 01_trimmed/           # step 2: trimmed reads, fastp reports, MultiQC
├── 02_quant/             # step 3: <sample>/abundance.tsv, quant_summary.tsv
├── 03_kb_count/          # step 4: <sample>/counts_unfiltered/, kb_count_summary.tsv
└── logs/
    ├── 1_build_index.log
    ├── 2_trim_bulk.log
    ├── 3_quant_bulk.log
    └── 4_kb_count.log
```

Each folder is exactly what the matching step-by-step guide describes, so
the downstream code there works unchanged:
[DESeq2 for bulk](04_bulk_kallisto_quant.md#7-next-differential-expression-with-deseq2),
[Seurat for single-cell](05_singlecell_kb_count.md#6-next-seurat).

Each step folder also has a `pipeline_settings.json`: the settings that step
last ran with (Snakemake uses it to notice changes, section 6).

## 6. Re-running, adding samples, changing settings

Run the same command again in every case; Snakemake works out what to do.

| What happened | What the next run does |
|---------------|------------------------|
| nothing changed | nothing ("Nothing to be done") |
| the run was interrupted | finishes the samples that were not done yet |
| you added FASTQ files for new samples | trims and quantifies only the new samples |
| you changed a setting (e.g. `strand`, `fastp_args`, `tech`) | redoes that step, and the steps after it, for **every** sample |
| you switched to another index | redoes all quantification |

To redo a step from scratch anyway, delete its folder (e.g.
`results/02_quant/`) and run again.

## 7. If something goes wrong

| Message or symptom | What to do |
|--------------------|------------|
| `kallisto-TE config: ...` | the settings file has a problem; the message says which option. Nothing has run yet |
| `Error in rule quant_bulk` (or another step) | open that step's log in `results/logs/`. It is the same output as running the script by hand, so the troubleshooting table of that step's guide applies |
| `bulk reads are in one sub-folder per sample, but trimming needs ...` | fastp needs all files in one folder: move them, or set `bulk.trim: false` to quantify the raw reads |
| `snakemake: command not found` | `conda activate kallisto_te`, or install Snakemake (section 1) |
| `LockException` / "Directory cannot be locked" | an earlier run was killed. Make sure no other run is using the folder, then `snakemake --unlock --configfile ...` |
| `IncompleteFilesException` | an earlier run stopped mid-step. Add `--rerun-incomplete`; finished samples are still kept |
| SLURM jobs fail at once | check `slurm_partition` / `slurm_account` in `profiles/slurm/config.yaml`; SLURM's own messages are under `.snakemake/slurm_logs/` |
| a job is killed for time or memory | raise `runtime` (minutes) or `mem_mb` for that step under `resources:` and run again; finished samples are kept |

For your methods section, cite Snakemake (Mölder *et al.*, 2021) along with
the tools of each step.
