#!/usr/bin/env python3
"""Simulate small FASTQ files from the toy index's cDNA, for test/run_pipeline_test.sh.

Usage: simulate_toy_reads.py CDNA_FASTA OUTDIR

Writes:
  OUTDIR/bulk/{Ctrl1,Ctrl2,KO1,KO2}_R1.fq.gz + _R2.fq.gz
      paired-end, reverse-stranded ("rf"), 5,000 pairs each; in KO the TE
      L1Md_A is 4x higher
  OUTDIR/singlecell/Mouse1_S1_L00{1,2}_R{1,2}_001.fastq.gz (+ I1)
      10x Genomics v3 layout: 30 cells x 200 UMIs over two lanes
Only the Python standard library is needed, plus kb-python's ngs_tools for
the 10x v3 barcode list (installed with kb-python).
"""
import gzip
import os
import random
import sys

cdna, out = sys.argv[1], sys.argv[2]
seqs, name = {}, None
for line in open(cdna):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs[name] = []
    elif name:
        seqs[name].append(line.upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
names = sorted(seqs)
comp = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(comp)[::-1]


def weights(cond):
    # genes (T1, T2) 3x more reads than each TE
    return [(4 if cond == "KO" and n.startswith("L1Md_A") else 1) * (3 if n in ("T1", "T2") else 1)
            for n in names]


def bulk(sample, cond, n, rng, d):
    os.makedirs(d, exist_ok=True)
    w = weights(cond)
    with gzip.open(f"{d}/{sample}_R1.fq.gz", "wt") as f1, gzip.open(f"{d}/{sample}_R2.fq.gz", "wt") as f2:
        for i in range(n):
            s = seqs[rng.choices(names, w)[0]]
            L = min(len(s), max(100, int(rng.gauss(200, 20))))
            st = rng.randint(0, len(s) - L)
            frag = s[st:st + L]
            # reverse-stranded: read 1 is antisense, read 2 is sense
            f1.write(f"@{sample}:{i} 1:N:0\n{rc(frag)[:50]}\n+\n{'F' * 50}\n")
            f2.write(f"@{sample}:{i} 2:N:0\n{frag[:50]}\n+\n{'F' * 50}\n")


def singlecell(d, rng):
    import ngs_tools
    os.makedirs(d, exist_ok=True)
    wl = ngs_tools.chemistry.get_chemistry("10xv3").whitelist_path
    with gzip.open(wl, "rt") as fh:
        barcodes = [next(fh).strip() for _ in range(30)]
    fh = {L: [gzip.open(f"{d}/Mouse1_S1_L00{L}_{r}_001.fastq.gz", "wt") for r in ("R1", "R2", "I1")]
          for L in (1, 2)}
    i = 0
    for bc in barcodes:
        for _ in range(200):
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            s = seqs[rng.choice(names)]
            st = rng.randint(0, len(s) - 90)
            r1, r2, i1 = fh[1 + i % 2]
            r1.write(f"@sc:{i} 1:N:0\n{bc}{umi}\n+\n{'F' * 28}\n")
            r2.write(f"@sc:{i} 2:N:0\n{s[st:st + 90]}\n+\n{'F' * 90}\n")
            i1.write(f"@sc:{i} 3:N:0\nACGTACGT\n+\nFFFFFFFF\n")
            i += 1
    for files in fh.values():
        for f in files:
            f.close()


rng = random.Random(1)
for sample, cond in [("Ctrl1", "Ctrl"), ("Ctrl2", "Ctrl"), ("KO1", "KO"), ("KO2", "KO")]:
    bulk(sample, cond, 5000, rng, f"{out}/bulk")
singlecell(f"{out}/singlecell", rng)
