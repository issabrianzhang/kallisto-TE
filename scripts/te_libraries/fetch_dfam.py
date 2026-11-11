#!/usr/bin/env python3
"""
fetch_dfam.py -- download a species' TE consensus library from Dfam

Dfam (https://www.dfam.org) is an open-access database of transposable
element families, released under CC0 (no copyright). This script asks the
Dfam API for every family that applies to one species: families defined
for the species itself, for its ancestor groups (e.g. Rodentia, Mammalia)
and for its descendants (e.g. strains). It writes:

  <out>.fa.gz           consensus sequences, kallisto-TE-ready names
                        (>NAME_Type-Subtype, e.g. >MMERVK10C_LTR-ERVK)
  <out>.annotation.tsv  one row per TE: Dfam accession, type, subtype,
                        full classification, clade(s), length

Use the FASTA directly with build_index.sh:
  bash scripts/index/build_index.sh --species mouse \
      --te-fasta te_libraries/dfam/mouse_Dfam4.0_curated.fa.gz ...

Usage:
  python3 scripts/te_libraries/fetch_dfam.py --species "Mus musculus" \
      --out te_libraries/dfam/mouse_Dfam4.0_curated
  python3 scripts/te_libraries/fetch_dfam.py --species "Arabidopsis thaliana" \
      --include-uncurated --out arabidopsis_Dfam4.0_uncurated

Options:
  --species NAME         scientific name (or use --taxid)
  --taxid ID             NCBI taxonomy id, e.g. 10090 for mouse
  --relatives MODE       ancestors | descendants | both (default: both)
  --include-uncurated    also include uncurated, de novo families (DR
                         accessions). Needed for species without curated
                         families, such as most plants. These have no
                         family name, so the Dfam accession is used.
  --out PREFIX           output prefix (".fa.gz" / ".annotation.tsv" added)

Only the Python standard library is needed.
"""
import argparse
import gzip
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request

API = "https://dfam.org/api"
PAGE = 1000
HEADERS = {"User-Agent": "kallisto-TE fetch_dfam.py", "Accept": "application/json"}


def get_json(path, params, retries=4):
    # Dfam rejects "+" for spaces, so encode them as %20
    url = f"{API}/{path}?{urllib.parse.urlencode(params, quote_via=urllib.parse.quote)}"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.load(r)
        except Exception as e:  # network hiccup or rate limit: wait and retry
            if attempt == retries - 1:
                sys.exit(f"ERROR: request failed: {url}\n  {e}")
            time.sleep(5 * (attempt + 1))


def resolve_taxid(name):
    d = get_json("taxa", {"name": name, "limit": 20})
    for t in d.get("taxa", []):
        if t["name"].lower() == name.lower():
            return t["id"], t["name"]
    hits = ", ".join(f'{t["name"]} ({t["id"]})' for t in d.get("taxa", [])[:10])
    sys.exit(f'ERROR: no exact Dfam taxon named "{name}". Close matches: {hits or "none"}')


def paged(params, total=None):
    """Yield every page of /families for these params.

    JSON pages report total_count; FASTA pages do not, so pass `total`."""
    start = 0
    while True:
        d = get_json("families", dict(params, limit=PAGE, start=start))
        yield d
        start += PAGE
        n = total if total is not None else d.get("total_count", 0)
        if start >= n:
            break


def clean(s):
    s = s.replace("/", "-").replace(" ", "_").replace("\t", "_")
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--species")
    g.add_argument("--taxid", type=int)
    ap.add_argument("--relatives", default="both", choices=["ancestors", "descendants", "both"])
    ap.add_argument("--include-uncurated", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    version = get_json("version", {})
    if a.species:
        taxid, taxname = resolve_taxid(a.species)
    else:
        taxid, taxname = a.taxid, str(a.taxid)

    params = {"clade": taxid, "clade_relatives": a.relatives}
    if a.include_uncurated:
        params["include_raw"] = "true"

    # 1) metadata (name, class, clades) for every family
    meta = {}
    for d in paged(params):
        for r in d["results"]:
            meta[r["accession"]] = r
    if not meta:
        sys.exit(f"ERROR: Dfam has no families for {taxname} (taxid {taxid}) with these settings."
                 + ("" if a.include_uncurated else " Try --include-uncurated."))

    # 2) consensus sequences (FASTA), matched to the metadata by accession
    seqs = {}
    for d in paged(dict(params, format="fasta"), total=len(meta)):
        body = d["body"] if isinstance(d, dict) and "body" in d else d
        acc, buf = None, []
        for line in io.StringIO(body):
            line = line.strip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc, buf = line[1:].split()[0].split(".")[0], []
            elif line:
                buf.append(line)
        if acc:
            seqs[acc] = "".join(buf)

    missing = sorted(set(meta) - set(seqs))
    if missing:
        sys.exit(f"ERROR: {len(missing)} families came back without a sequence, e.g. {missing[:5]}")

    # 3) write FASTA + annotation table (sorted by accession: stable output)
    used = {}
    fa_lines, tsv_lines = [], ["te_id\tdfam_accession\tdfam_version\tname\trepeat_type\trepeat_subtype"
                               "\tclassification\tclades\tlength\tcurated"]
    for acc in sorted(meta):
        r = meta[acc]
        name = r.get("name") or acc
        rtype = r.get("repeat_type_name") or "Unknown"
        sub = r.get("repeat_subtype_name")
        te_id = clean(f"{name}_{rtype}" + (f"-{sub}" if sub else ""))
        if te_id in used:                    # keep names unique
            te_id = clean(f"{te_id}_{acc}")
        used[te_id] = acc
        seq = seqs[acc].upper()
        fa_lines.append(f">{te_id}")
        fa_lines.extend(seq[i:i + 60] for i in range(0, len(seq), 60))
        clades = r.get("clades") or []
        clades = ";".join(c.split(";")[-1] for c in clades) if isinstance(clades, list) else str(clades)
        tsv_lines.append("\t".join(map(str, [
            te_id, acc, r.get("version", ""), r.get("name") or "", rtype, sub or "",
            r.get("classification", ""), clades, len(seq), "yes" if acc.startswith("DF") else "no"])))

    # Empty file name and mtime=0 in the gzip header: the same Dfam release
    # always gives byte-identical files, and no local path is stored.
    with open(a.out + ".fa.gz", "wb") as raw, \
            gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as f:
        f.write(("\n".join(fa_lines) + "\n").encode())
    with open(a.out + ".annotation.tsv", "w") as f:
        f.write("\n".join(tsv_lines) + "\n")

    n_cur = sum(acc.startswith("DF") for acc in meta)
    non_acgt = sum(len(re.sub(r"[ACGT]", "", s.upper())) for s in seqs.values())
    print(f"Dfam {version.get('dfam_version')} ({version.get('release_date')}): {taxname} (taxid {taxid}), "
          f"relatives={a.relatives}")
    print(f"  {len(meta)} families ({n_cur} curated, {len(meta) - n_cur} uncurated)")
    print(f"  {non_acgt} non-ACGT bases (IUPAC codes; kallisto replaces them when indexing)")
    print(f"  wrote {a.out}.fa.gz and {a.out}.annotation.tsv")


if __name__ == "__main__":
    main()
