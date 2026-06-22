"""Build the hg19 and hg38 fragile-sites references from HumCFS.

Source: HumCFS database (https://webs.iiitd.edu.in/raghava/humcfs/), per-chromosome
BED files in GRCh38/hg38 coordinates. This script:

  * cleans the native hg38 coordinates (chrX casing, clamps terminal-band ends that
    overshoot the chromosome, drops corrupt records) -> hg38 reference;
  * lifts the same sites over to hg19 with the UCSC hg38->hg19 chain (pyliftover).

Both are emitted to cns/utils/fragile_sites.py (the single source of truth) so they
can be attached to an Assembly the same way cytobands are. Pass --bed to also dump
4-column BEDs to ./data (e.g. for use as a pipeline/CLI mask).

Run from the repo root:  python scripts/build_fragile_sites.py [--bed]
"""

import io
import os
import sys
import urllib.request
import zipfile
from glob import glob

from pyliftover import LiftOver

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from cns.utils.genomes import hg19_chr_lengths, hg38_chr_lengths

HUMCFS_URL = "https://webs.iiitd.edu.in/raghava/humcfs/fragile_site_bed.zip"
HUMCFS_DIR = os.path.join(REPO, "scripts", "_humcfs_hg38")  # extracted HumCFS beds
DATA_DIR = os.path.join(REPO, "data")
MODULE = os.path.join(REPO, "cns", "utils", "fragile_sites.py")

# HumCFS chrX file uses lowercase "chrx"; normalize to UCSC casing.
CHROM_FIX = {"chrx": "chrX", "chry": "chrY", "chrm": "chrM"}

# Endpoints at telomeres (pos 0) or inside assembly gaps do not lift directly.
# Nudge inward in small steps to find the nearest mappable base. These are
# megabase-scale cytoband regions, so a sub-Mb boundary shift is negligible.
NUDGE_STEP = 1000
NUDGE_MAX = 600000


def fetch_humcfs():
    """Download and extract the HumCFS hg38 per-chromosome BEDs if not present."""
    if glob(os.path.join(HUMCFS_DIR, "chr*_fragile_site.bed")):
        return
    os.makedirs(HUMCFS_DIR, exist_ok=True)
    print(f"Downloading HumCFS beds from {HUMCFS_URL} ...")
    with urllib.request.urlopen(HUMCFS_URL) as resp:
        data = resp.read()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        for member in zf.namelist():
            base = os.path.basename(member)
            if base.endswith("_fragile_site.bed") and "MACOSX" not in member:
                with zf.open(member) as src, open(os.path.join(HUMCFS_DIR, base), "wb") as dst:
                    dst.write(src.read())


def chrom_sort_key(chrom):
    """Sort chromosomes naturally: chr1..chr22, chrX, chrY."""
    name = chrom[3:]
    if name.isdigit():
        return (0, int(name))
    return (1, name)


def read_humcfs():
    """Read all HumCFS hg38 beds -> list of (chrom, start, end, name)."""
    sites = []
    for path in glob(os.path.join(HUMCFS_DIR, "chr*_fragile_site.bed")):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line.lower().startswith("chr"):
                    continue
                fields = line.split("\t")
                chrom = CHROM_FIX.get(fields[0].lower(), fields[0])
                start, end, name = int(fields[1]), int(fields[2]), fields[3]
                sites.append((chrom, start, end, name))
    return sites


def clean_hg38(sites):
    """Validate native hg38 coordinates: clamp terminal overshoots, drop corrupt records."""
    cleaned, dropped = [], []
    for chrom, start, end, name in sites:
        if chrom not in hg38_chr_lengths:
            dropped.append((chrom, start, end, name, "unknown chrom"))
            continue
        hg38_len = hg38_chr_lengths[chrom]
        # Some terminal-band ends carry the hg19 length (a known HumCFS quirk); if the
        # end overshoots hg38 only up to the hg19 length, the band reaches the telomere
        # -> clamp to the hg38 chromosome length. Larger overshoots are corrupt.
        if end > hg38_len:
            if end <= hg19_chr_lengths[chrom]:
                end = hg38_len
            else:
                dropped.append((chrom, start, end, name, "end beyond chromosome"))
                continue
        if start < 0 or end <= start:
            dropped.append((chrom, start, end, name, "zero/neg length"))
            continue
        cleaned.append((chrom, start, end, name))
    return cleaned, dropped


def lift_pos(lo, chrom, pos, direction):
    """Lift a single coordinate, nudging inward (direction +1/-1) past unmappable bases."""
    for offset in range(0, NUDGE_MAX + 1, NUDGE_STEP):
        res = lo.convert_coordinate(chrom, pos + direction * offset)
        if res and res[0][0] == chrom:
            return res[0][1]
    return None


def lift(sites, lo):
    """Lift over (chrom, start, end) hg38 -> hg19, nudging unmappable endpoints inward."""
    lifted, dropped = [], []
    for chrom, start, end, name in sites:
        s_pos = lift_pos(lo, chrom, start, +1)  # start nudges forward (into the region)
        e_pos = lift_pos(lo, chrom, end, -1)    # end nudges backward (into the region)
        # Terminal bands: some HumCFS ends overshoot the hg38 chromosome (their value
        # is the hg19 length). If the end sits between the two assembly lengths, the
        # band reaches the telomere -> clamp to the hg19 chromosome length.
        if e_pos is None and chrom in hg38_chr_lengths:
            if hg38_chr_lengths[chrom] < end <= hg19_chr_lengths[chrom]:
                e_pos = hg19_chr_lengths[chrom]
        if s_pos is None or e_pos is None:
            dropped.append((chrom, start, end, name, "unmapped"))
            continue
        lo_pos, hi_pos = sorted((s_pos, e_pos))
        if hi_pos <= lo_pos:
            dropped.append((chrom, start, end, name, "zero/neg length"))
            continue
        lifted.append((chrom, lo_pos, hi_pos, name))
    return lifted, dropped


def report(label, sites, dropped):
    print(f"{label}: kept {len(sites)}, dropped {len(dropped)}.")
    for d in dropped:
        print("  dropped:", d)


def write_bed(path, sites):
    with open(path, "w", newline="\n") as f:
        for chrom, start, end, name in sites:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    print(f"Wrote {path}")


MODULE_DOCSTRING = '''"""HumCFS common fragile sites for hg19 and hg38.

This file is generated by ``scripts/build_fragile_sites.py`` -- do not edit by hand.
Regenerate with: ``python scripts/build_fragile_sites.py``

Source
------
HumCFS database (https://webs.iiitd.edu.in/raghava/humcfs/), the per-chromosome
``fragile_site_bed.zip`` release, in GRCh38/hg38 coordinates ({n_raw} sites spanning
chr1-22 and chrX; HumCFS provides no chr21 or chrY entries).

Processing
----------
hg38 ({n_hg38} sites) -- the cleaned native HumCFS coordinates:
  * the chrX file labels its records ``chrx``; normalized to UCSC ``chrX``;
  * a few terminal-band ends carry the hg19 chromosome length (a HumCFS quirk) and so
    overshoot the hg38 chromosome -- these are clamped to the hg38 chromosome length;
  * one record with a corrupt end far beyond the chromosome is dropped (FRA10D).

hg19 ({n_hg19} sites) -- the same sites lifted over with the UCSC hg38->hg19 chain
(via pyliftover):
  * endpoints landing in telomeres (position 0) or assembly gaps are nudged inward in
    1 kb steps (up to 600 kb) to the nearest mappable base;
  * terminal-band ends that overshoot hg38 are clamped to the hg19 telomere;
  * 4 sites are dropped -- FRA10D (corrupt), and FRA1F / FRA1J / FRA9F, which lie in the
    1q12 / 9q12 heterochromatin gaps and have no reliable hg38<->hg19 alignment.

Tuples are ``(chrom, start, end, name)``, 0-based half-open, matching the cytobands
module. They are attached to the corresponding ``Assembly`` (see cns/utils/assemblies.py)
and exposed per chromosome via ``cns.data_utils.load_fragile_sites(assembly)``.
"""
'''


def write_module(by_assembly, n_raw):
    with open(MODULE, "w", newline="\n") as f:
        f.write(MODULE_DOCSTRING.format(
            n_raw=n_raw,
            n_hg38=len(by_assembly["hg38"]),
            n_hg19=len(by_assembly["hg19"]),
        ))
        for assembly, sites in by_assembly.items():
            f.write(f"\n{assembly}_fragile_sites = (\n")
            for chrom, start, end, name in sites:
                f.write(f'    ("{chrom}", {start}, {end}, "{name}"),\n')
            f.write(")\n")
    print(f"Wrote {MODULE}")


def main():
    fetch_humcfs()
    sites = read_humcfs()
    print(f"Read {len(sites)} HumCFS hg38 fragile sites.\n")

    # hg38: native coordinates, cleaned.
    hg38_sites, hg38_dropped = clean_hg38(sites)
    hg38_sites.sort(key=lambda r: (chrom_sort_key(r[0]), r[1], r[2]))
    report("hg38", hg38_sites, hg38_dropped)

    # hg19: lifted over.
    print("\nLoading hg38->hg19 chain (downloads on first run)...")
    lo = LiftOver("hg38", "hg19")
    hg19_sites, hg19_dropped = lift(sites, lo)
    hg19_sites.sort(key=lambda r: (chrom_sort_key(r[0]), r[1], r[2]))
    report("hg19", hg19_sites, hg19_dropped)

    write_module({"hg19": hg19_sites, "hg38": hg38_sites}, n_raw=len(sites))

    if "--bed" in sys.argv:
        os.makedirs(DATA_DIR, exist_ok=True)
        write_bed(os.path.join(DATA_DIR, "HumCFS_fragile_sites_hg38.bed"), hg38_sites)
        write_bed(os.path.join(DATA_DIR, "HumCFS_fragile_sites_hg19.bed"), hg19_sites)


if __name__ == "__main__":
    main()
