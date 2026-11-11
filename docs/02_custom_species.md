# Step 1b — Build an index for your own species

The mouse and human presets are only shortcuts. The method works for **any
organism with a genome assembly, a gene annotation and a set of TE consensus
sequences**: plants, fish, insects, worms, fungi, ...

This page uses ***Arabidopsis thaliana* (TAIR10)** as the worked example.
Every step says what to change for other species. Read
[01_build_index.md](01_build_index.md) first for the concepts.

---

## Contents

1. [The three ingredients](#1-the-three-ingredients)
2. [Get the genome and gene annotation](#2-get-the-genome-and-gene-annotation)
3. [Get TE consensus sequences](#3-get-te-consensus-sequences)
4. [Build the index](#4-build-the-index)
5. [Check the result](#5-check-the-result)
6. [Tips for non-model and plant genomes](#6-tips-for-non-model-and-plant-genomes)
7. [Add your species as a preset (optional)](#7-add-your-species-as-a-preset-optional)

---

## 1. The three ingredients

| Ingredient | Format | Arabidopsis example | Requirement |
|------------|--------|---------------------|-------------|
| Genome | FASTA (`.fa` / `.fa.gz`) | Ensembl Plants TAIR10 `dna.toplevel` | Chromosome names must match the GTF |
| Gene annotation | **GTF** (`.gtf` / `.gtf.gz`) | Ensembl Plants TAIR10 GTF | Needs `exon` lines with `gene_id` and `transcript_id` |
| TE consensus library | FASTA | RepBase, Dfam or EDTA/RepeatModeler2 | One sequence per TE family/subfamily |

Make a working folder:

```bash
mkdir -p ~/references/TAIR10 && cd ~/references/TAIR10
```

---

## 2. Get the genome and gene annotation

Good sources, by group of organisms:

| Organisms | Where |
|-----------|-------|
| Vertebrates | [Ensembl](https://www.ensembl.org/info/data/ftp/index.html), [GENCODE](https://www.gencodegenes.org/) (human/mouse) |
| Plants | [Ensembl Plants](https://plants.ensembl.org/info/data/ftp/index.html), [Phytozome](https://phytozome-next.jgi.doe.gov/) |
| Insects, worms, other animals | [Ensembl Metazoa](https://metazoa.ensembl.org/info/data/ftp/index.html), [FlyBase](https://flybase.org/), [WormBase](https://wormbase.org/) |
| Fungi, protists | [Ensembl Fungi](https://fungi.ensembl.org/) / [Ensembl Protists](https://protists.ensembl.org/) |
| Anything else | [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/) |

**Arabidopsis (Ensembl Plants release 63; check for a newer release at
<https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/>):**

```bash
REL=63
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${REL}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${REL}/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${REL}.gtf.gz
```

Use the plain `dna` genome, **not** the repeat-masked `dna_rm` or `dna_sm`
versions. Masking would remove the very TE sequences you want to measure.

### Only have a GFF3?

Many sources (NCBI, Phytozome, TAIR) ship GFF3. Convert it to GTF with
[gffread](https://github.com/gpertea/gffread):

```bash
mamba install -n kallisto_te -c bioconda gffread
gffread -T annotation.gff3 -o annotation.gtf
```

### Check that the names match

The first word of each FASTA header must match column 1 of the GTF:

```bash
gzip -dc Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz | grep '>' | cut -d' ' -f1 | head
#   >1  >2  >3  >4  >5  >Mt  >Pt
gzip -dc Arabidopsis_thaliana.TAIR10.${REL}.gtf.gz | grep -v '^#' | cut -f1 | sort -u | head
#   1  2  3  4  5  Mt  Pt
```

If one says `chr1` and the other `1`, fix one of them before building.
Genes on contigs missing from the FASTA are skipped by `kb ref`.

---

## 3. Get TE consensus sequences

Pick **one** of these routes. You can also combine several (§4).

### Route A: RepBase (curated; needs a free academic registration)

RepBase covers many plants and animals. You don't need to know which
`*.ref` file holds your species: `extract_repbase_species.sh` searches the
species field of every RepBase header.

```bash
tar -xzf RepBase28.04.fasta.tar.gz

# 1) How many sequences exist for my species, and in which files?
bash scripts/index/extract_repbase_species.sh \
    -d RepBase28.04.fasta -p "Arabidopsis thaliana" -l

# 2) Extract them
bash scripts/index/extract_repbase_species.sh \
    -d RepBase28.04.fasta -p "Arabidopsis thaliana" -o ath_RepBase28.04_raw.fa
```

Notes:

* The match is case-insensitive and a substring, so `-p "Arabidopsis"` also
  picks up *A. lyrata* entries. Use `-x` for an exact species match.
* Families shared across a lineage are often filed under a higher group,
  just as the mouse preset adds the general `rodrep.ref`. Check with `-l`
  whether entries exist for your genus, family or a broader group (e.g.
  `-p "Brassicaceae"`), and add them if you want them counted.
* The same entry can appear in several `*.ref` files. Exact duplicates are
  removed automatically in the next step.

### Route B: Dfam (open access)

[Dfam](https://www.dfam.org/) is an open database of TE families (HMMs and
consensus sequences). Export your species' families as a consensus FASTA
with Dfam's `famdb.py` tool (see the Dfam documentation for the current
command). Dfam is especially complete for vertebrates and model animals.

### Route C: Build a de novo library (any genome, no database needed)

For non-model species, or when databases are incomplete, build a consensus
library directly from the genome:

| Tool | Output consensus FASTA | Good for |
|------|------------------------|----------|
| [EDTA](https://github.com/oushujun/EDTA) | `<genome>.mod.EDTA.TElib.fa` | **plants** (widely used) |
| [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler) | `<db>-families.fa` | any eukaryote |

These headers look like `>rnd-1_family-1#LTR/Gypsy ( Recon Family Size = 25 )`.
`prepare_te_library.sh` keeps only the first word and cleans symbols, so this
becomes `rnd-1_family-1_LTR-Gypsy`. The TE class stays in the name.

> **Not suitable:** TE *locus* annotations, such as TAIR10's
> `TAIR10_Transposable_Elements.txt` or a RepeatMasker `.out` table, list
> genomic coordinates, not consensus sequences. kallisto-TE needs one
> representative sequence per family.

---

## 4. Build the index

```bash
conda activate kallisto_te

bash scripts/index/build_index.sh \
    --species custom \
    --name    TAIR10_Ensembl63_RepBase28.04 \
    --genome  Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
    --gtf     Arabidopsis_thaliana.TAIR10.63.gtf.gz \
    --te-fasta ath_RepBase28.04_raw.fa \
    --outdir  ~/references/TAIR10/TAIR10_TE_index \
    --threads 8
```

Other ways to give the TE library:

```bash
# Directly from the RepBase archive, naming the *.ref file(s) to use
# ("athrep.ref" is an example: find the real file names with extract_repbase_species.sh -l)
    --repbase RepBase28.04.fasta.tar.gz --te-libs "athrep.ref"

# Several sources at once (e.g. curated + de novo); all are merged and de-duplicated
    --te-fasta ath_RepBase28.04_raw.fa --te-fasta TAIR10.mod.EDTA.TElib.fa
```

> When merging a curated and a de novo library, the same TE family can be
> present under two different names. kallisto then splits reads between
> the two (multi-mapping). Prefer one well-curated source, or remove
> redundant families first (e.g. with `cd-hit-est`).

On a cluster, edit and submit
[`slurm/build_index_custom_species.sbatch`](../slurm/build_index_custom_species.sbatch).

**Resources.** Memory scales with the transcriptome + TE library size,
not the genome size alone: Arabidopsis needs ~8–16 GB and well under an
hour. Mammals need 32–64 GB. Very large plant genomes (maize, wheat) with
large de novo TE libraries may need more.

---

## 5. Check the result

The build ends with an automatic check:

```
OK: all <N> TEs are present in t2g.txt
OK: all <N> TEs are present in cdna.fa
```

Then look at a few entries:

```bash
cd ~/references/TAIR10/TAIR10_TE_index
grep -c 'TE:' index/t2g.txt                 # number of TE families in the index
grep 'TE:' index/t2g.txt | head -3
grep -v 'TE:' index/t2g.txt | head -3       # normal genes look as usual
cat build_info.txt                          # inputs and versions, for your methods section
```

---

## 6. Tips for non-model and plant genomes

* **Name your index after all its versions**, e.g.
  `TAIR10_Ensembl63_RepBase28.04`. You will need them for the methods
  section.
* **Organelle genomes** (`Mt`, `Pt`) can stay in. They carry real genes and
  no TE consensus will clash with them.
* **Many "Unknown" families** in de novo libraries are normal. They are
  kept, and you can filter them during analysis.
* **Polyploids / very repetitive genomes:** consensus-level counting
  already collapses all copies of a family into one feature, which is the
  point of the method. Expect more multi-mapping between closely related
  subfamilies.
* **Single-end bulk data and short TEs:** many TE consensus sequences
  (e.g. solo LTRs) are only a few hundred bp. With `kallisto quant --single`,
  reads near the end of short targets are dropped unless you add
  `--single-overhang`. This is covered in the quantification step.

---

## 7. Add your species as a preset (optional)

If you build indices for a species repeatedly, add it to the `case` block
near the top of [`scripts/index/build_index.sh`](../scripts/index/build_index.sh):

```bash
    arabidopsis)
        TE_LIBS="${TE_LIBS:-athrep.ref}"   # RepBase file(s) for this species (example name)
        GENCODE_GENOME_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
        GENCODE_GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.63.gtf.gz"
        ;;
```

(The variables are called `GENCODE_*` but take any URL.) Check the exact
RepBase file name with `extract_repbase_species.sh -l` first. Pull requests
adding tested presets are welcome!
