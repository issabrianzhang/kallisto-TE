# Ready-to-use Dfam TE libraries

TE consensus libraries for kallisto-TE, downloaded from
[Dfam](https://www.dfam.org) **release 4.0 (May 2026)**. Dfam is open access
and released under **CC0 1.0** ("no copyright"), so these files can be
shared freely and used without any registration. If you have no access to
RepBase, start here.

Each library contains every Dfam family that applies to the species: families
defined for the species itself, its ancestor groups (e.g. *Rodentia*,
*Mammalia*, *Eutheria*) and its descendants.

| Species | File | Families | Curated? | Largest classes |
|---------|------|---------:|----------|-----------------|
| Mouse (*Mus musculus*) | `mouse_Dfam4.0_curated.fa.gz` | 1,362 | curated | LTR 638, DNA 292, LINE 251 |
| Human (*Homo sapiens*) | `human_Dfam4.0_curated.fa.gz` | 1,403 | curated | LTR 617, DNA 329, LINE 194 |
| Zebrafish (*Danio rerio*) | `zebrafish_Dfam4.0_curated.fa.gz` | 1,965 | curated | LTR 770, DNA 757, LINE 282 |
| Fruit fly (*Drosophila melanogaster*) | `fly_Dfam4.0_curated.fa.gz` | 390 | curated | LTR 282, LINE 43, DNA 38 |
| Worm (*Caenorhabditis elegans*) | `worm_Dfam4.0_curated.fa.gz` | 187 | curated | DNA 93, LTR 42, LINE 13 |
| Arabidopsis (*A. thaliana*) | `arabidopsis_Dfam4.0_uncurated.fa.gz` | 503 | **uncurated** (de novo) | DNA 154, LTR 140, Unknown 103 |

## Use

```bash
# Mouse or human: --dfam picks the bundled file automatically
bash scripts/index/build_index.sh --species mouse --dfam \
    --name mm10_GencodeM25_Dfam4.0 --outdir /path/to/references/mm10_Dfam_TE --threads 16

# Any other species: pass the file with --te-fasta
bash scripts/index/build_index.sh --species custom \
    --genome GRCz11.fa.gz --gtf GRCz11.gtf.gz \
    --te-fasta te_libraries/dfam/zebrafish_Dfam4.0_curated.fa.gz \
    --name GRCz11_Dfam4.0 --outdir /path/to/references/GRCz11_TE
```

## Files

* **`<species>_Dfam4.0_<curated|uncurated>.fa.gz`**: consensus sequences.
  Names are `NAME_Type-Subtype`, e.g. `MMERVK10C_LTR-ERVK`,
  `IAPEz-int_LTR-ERVK`, `L1MdA_I_5end_LINE-L1`, `B1_Mm_SINE-Alu`. These
  appear as TE "genes" (`TE:MMERVK10C_LTR-ERVK`) in `t2g.txt` and in your
  count tables.
* **`<species>_Dfam4.0_*.annotation.tsv`**: one row per TE, for grouping and
  filtering in R/Python:

  | Column | Example |
  |--------|---------|
  | `te_id` | `MMERVK10C_LTR-ERVK` (the name used in the index) |
  | `dfam_accession`, `dfam_version` | `DF000001942`, `1` |
  | `name` | `MMERVK10C` |
  | `repeat_type`, `repeat_subtype` | `LTR`, `ERVK` |
  | `classification` | full Dfam classification path |
  | `clades` | taxon the family is defined for, e.g. `Mus`, `Rodentia`, `Mammalia` |
  | `length` | consensus length (bp), e.g. `7581` |
  | `curated` | `yes` (DF accession) / `no` (DR accession) |

  ```r
  te_annot <- read.delim("te_libraries/dfam/mouse_Dfam4.0_curated.annotation.tsv")
  # e.g. join to a count table whose rownames are "TE:<te_id>", then sum by repeat_type
  ```

## Things to know

* **Dfam and RepBase names differ** (e.g. Dfam `MMERVK10C_LTR-ERVK` vs RepBase
  `MMERVK10C_ERVK_Mus_musculus`), and the family sets are not identical
  (mouse: 1,362 Dfam families vs 1,740 in our RepBase 25.05 library).
  Results from a Dfam index and a RepBase index are therefore not directly
  interchangeable. Pick one library per project.
* **Small RNA and satellite families are included**, because Dfam (like
  RepBase) lists them as repeats: tRNA, 5S rRNA, 7SK, Y RNA, satellites.
  Some are also GENCODE genes, so reads can be split between the gene and
  the Dfam entry. To drop them, filter with the annotation table:

  ```bash
  awk -F'\t' 'NR>1 && $5 !~ /^(tRNA|rRNA|snRNA|scRNA|srpRNA|Satellite)$/ {print $1}' \
      te_libraries/dfam/mouse_Dfam4.0_curated.annotation.tsv > keep.txt
  gzip -dc te_libraries/dfam/mouse_Dfam4.0_curated.fa.gz \
      | awk 'NR==FNR{k[$1]; next} /^>/{p = (substr($1,2) in k)} p' keep.txt - > mouse_Dfam4.0_TE_only.fa
  ```

  | Library | small-RNA families | satellites |
  |---------|------:|------:|
  | mouse | 72 | 8 |
  | human | 71 | 36 |
  | zebrafish | 70 | 39 |
  | fly | 7 | 11 |
  | worm | 6 | 17 |
  | arabidopsis | 3 | 6 |

* **Arabidopsis has no curated Dfam families.** The Arabidopsis file holds
  Dfam's 503 *uncurated*, de novo families. They have no family name (the
  Dfam accession is used, e.g. `DR003462015_LTR-Gypsy`), and many are
  classified only as `Unknown`. For plants, a curated RepBase library or an
  EDTA library ([docs/02_custom_species.md](../../docs/02_custom_species.md))
  is usually better.
* Consensus sequences contain some IUPAC ambiguity codes (e.g. `R`, `Y`,
  `N`). kallisto replaces these when it builds the index. This is normal.

## Other species / regenerating

Any species in Dfam can be fetched with the same script (Python standard
library only):

```bash
python3 scripts/te_libraries/fetch_dfam.py --species "Rattus norvegicus" \
    --out te_libraries/dfam/rat_Dfam4.0_curated

# species without curated families (most plants, many non-model animals)
python3 scripts/te_libraries/fetch_dfam.py --species "Oryza sativa" \
    --include-uncurated --out rice_Dfam4.0_uncurated
```

The bundled files were made with:

```bash
python3 scripts/te_libraries/fetch_dfam.py --species "Mus musculus"            --out te_libraries/dfam/mouse_Dfam4.0_curated
python3 scripts/te_libraries/fetch_dfam.py --species "Homo sapiens"            --out te_libraries/dfam/human_Dfam4.0_curated
python3 scripts/te_libraries/fetch_dfam.py --species "Danio rerio"             --out te_libraries/dfam/zebrafish_Dfam4.0_curated
python3 scripts/te_libraries/fetch_dfam.py --species "Drosophila melanogaster" --out te_libraries/dfam/fly_Dfam4.0_curated
python3 scripts/te_libraries/fetch_dfam.py --species "Caenorhabditis elegans"  --out te_libraries/dfam/worm_Dfam4.0_curated
python3 scripts/te_libraries/fetch_dfam.py --species "Arabidopsis thaliana" --include-uncurated --out te_libraries/dfam/arabidopsis_Dfam4.0_uncurated
```

The same Dfam release always produces byte-identical files. MD5 checksums
are in [`MD5SUMS`](MD5SUMS).

## Citation

Storer J, Hubley R, Rosen J, Wheeler TJ, Smit AF. The Dfam community resource
of transposable element families, sequence models, and genome annotations.
*Mob DNA* 12, 2 (2021).
[doi:10.1186/s13100-020-00230-y](https://doi.org/10.1186/s13100-020-00230-y)
