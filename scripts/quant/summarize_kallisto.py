#!/usr/bin/env python3
"""Combine per-sample kallisto quant results into gene + TE tables.

Usage: summarize_kallisto.py QUANT_DIR T2G

Reads QUANT_DIR/samples_used.tsv (written by run_kallisto_quant.sh) and each
QUANT_DIR/<sample>/abundance.tsv + run_info.json. Writes:

  QUANT_DIR/matrices/gene_TE_counts.tsv  estimated counts summed per gene / TE
  QUANT_DIR/matrices/gene_TE_tpm.tsv     TPM summed per gene / TE
  QUANT_DIR/matrices/TE_counts.tsv       only the TE rows of gene_TE_counts.tsv
  QUANT_DIR/quant_summary.tsv            one line per sample

A TE is a target whose gene name in t2g.txt starts with "TE:" (that is how
build_index.sh labels them). These tables are for looking at the data and for
quick plots; for differential expression, import the per-sample abundance.tsv
files with tximport (see docs/04_bulk_kallisto_quant.md).
Only the Python standard library is used.
"""
import csv
import json
import sys
from pathlib import Path


def main(quant_dir: str, t2g_path: str) -> None:
    qd = Path(quant_dir)

    # transcript -> (gene_id, gene_name); keep genes in t2g order
    tx2gene, gene_name, gene_order = {}, {}, []
    with open(t2g_path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 2:
                continue
            tx, gid = row[0], row[1]
            name = row[2] if len(row) > 2 and row[2] else gid
            tx2gene[tx] = gid
            if gid not in gene_name:
                gene_name[gid] = name
                gene_order.append(gid)
    is_te = {g: gene_name[g].startswith("TE:") for g in gene_order}

    with open(qd / "samples_used.tsv") as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))

    counts, tpm, summary = {}, {}, []
    for s in samples:
        name = s["sample"]
        ab = qd / name / "abundance.tsv"
        info = json.load(open(qd / name / "run_info.json"))
        c, t = {}, {}
        missing = 0
        with open(ab) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gid = tx2gene.get(row["target_id"])
                if gid is None:
                    missing += 1
                    continue
                c[gid] = c.get(gid, 0.0) + float(row["est_counts"])
                t[gid] = t.get(gid, 0.0) + float(row["tpm"])
        if missing:
            sys.exit(f"ERROR: {missing} targets in {ab} are not in {t2g_path}. "
                     "Was this sample quantified with a different index?")
        counts[name], tpm[name] = c, t
        assigned = sum(c.values())
        te = sum(v for g, v in c.items() if is_te[g])
        summary.append([
            name, s["layout"], s["strand"],
            info["n_processed"], info["n_pseudoaligned"], f'{info["p_pseudoaligned"]:.1f}',
            f"{te:.0f}", f"{100 * te / assigned:.2f}" if assigned else "NA",
        ])

    out = qd / "matrices"
    out.mkdir(exist_ok=True)
    names = [s["sample"] for s in samples]

    def write(path, table, keep):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["gene_id", "gene_name", "type"] + names)
            for g in gene_order:
                if keep(g):
                    w.writerow([g, gene_name[g], "TE" if is_te[g] else "gene"]
                               + [f"{table[n].get(g, 0.0):.3f}" for n in names])

    write(out / "gene_TE_counts.tsv", counts, lambda g: True)
    write(out / "gene_TE_tpm.tsv", tpm, lambda g: True)
    write(out / "TE_counts.tsv", counts, lambda g: is_te[g])

    with open(qd / "quant_summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "layout", "strand", "reads_processed", "reads_pseudoaligned",
                    "percent_pseudoaligned", "TE_reads", "percent_TE"])
        w.writerows(summary)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
