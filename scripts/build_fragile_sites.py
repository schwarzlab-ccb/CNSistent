"""Build the hg19 fragile-sites reference from HumCFS.

Source: HumCFS database (https://webs.iiitd.edu.in/raghava/humcfs/), per-chromosome
BED files in GRCh38/hg38 coordinates. These are lifted over to hg19 with the UCSC
hg38->hg19 chain (via pyliftover) and written as a single 4-column BED to ./data,
mirroring how the gene reference BEDs (COSMIC/ENSEMBL) are stored.

It also emits cns/utils/fragile_sites.py, so the sites can be attached to an
Assembly the same way cytobands are.

Run from the repo root:  python scripts/build_fragile_sites.py
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
DATA_BED = os.path.join(REPO, "data", "HumCFS_fragile_sites.bed")
MODULE = os.path.join(REPO, "cns", "utils", "fragile_sites.py")


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


# HumCFS chrX file uses lowercase "chrx"; normalize to UCSC casing.
CHROM_FIX = {"chrx": "chrX", "chry": "chrY", "chrm": "chrM"}

# Endpoints at telomeres (pos 0) or inside assembly gaps do not lift directly.
# Nudge inward in small steps to find the nearest mappable base. These are
# megabase-scale cytoband regions, so a sub-100kb boundary shift is negligible.
NUDGE_STEP = 1000
NUDGE_MAX = 600000


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


def main():
    fetch_humcfs()

    print("Loading hg38->hg19 chain (downloads on first run)...")
    lo = LiftOver("hg38", "hg19")

    sites = read_humcfs()
    print(f"Read {len(sites)} HumCFS hg38 fragile sites.")

    lifted, dropped = lift(sites, lo)
    lifted.sort(key=lambda r: (chrom_sort_key(r[0]), r[1], r[2]))
    print(f"Lifted {len(lifted)} sites to hg19; dropped {len(dropped)}.")
    for d in dropped:
        print("  dropped:", d)

    os.makedirs(os.path.dirname(DATA_BED), exist_ok=True)
    with open(DATA_BED, "w", newline="\n") as f:
        for chrom, start, end, name in lifted:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    print(f"Wrote {DATA_BED}")

    with open(MODULE, "w", newline="\n") as f:
        f.write('"""HumCFS common fragile sites, lifted to hg19.\n\n')
        f.write("Source: HumCFS (https://webs.iiitd.edu.in/raghava/humcfs/), GRCh38 per-chromosome\n")
        f.write("BED files lifted to hg19 via the UCSC hg38->hg19 chain. Generated by\n")
        f.write("scripts/build_fragile_sites.py. Tuples are (chrom, start, end, name), 0-based\n")
        f.write('half-open, matching the cytobands module.\n"""\n\n')
        f.write("hg19_fragile_sites = (\n")
        for chrom, start, end, name in lifted:
            f.write(f'    ("{chrom}", {start}, {end}, "{name}"),\n')
        f.write(")\n")
    print(f"Wrote {MODULE}")


if __name__ == "__main__":
    main()
