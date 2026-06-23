import numpy as np
import pandas as pd

from cns.utils.assemblies import hg19
from cns.utils.canonization import get_cn_cols
from cns.utils.selection import only_aut


def calc_lengths(cns_df):
    return cns_df["end"] - cns_df["start"]


def calc_mid(cns_df):
    return cns_df["start"] + calc_lengths(cns_df) // 2


def calc_cum_mid(cns_df, assembly=hg19):
    mid = calc_mid(cns_df)
    offset = cns_df.apply(lambda x: assembly.chr_starts[x["chrom"]], axis=1)
    return mid + offset


def calc_nan_cols(cns_df, cn_columns=None):
    cn_columns = get_cn_cols(cns_df, cn_columns)
    return cns_df[cn_columns].isna().any(axis=1)


def cytobands_to_df(cytobands):
    return pd.DataFrame(cytobands, columns=["chrom", "start", "end", "name", "stain"])


def gaps_to_df(gaps):
    return pd.DataFrame(gaps, columns=["chrom", "start", "end", "type", "bridge"])


def fragile_sites_to_df(fragile_sites):
    return pd.DataFrame(fragile_sites, columns=["chrom", "start", "end", "name"])


def segments_to_cns_df(segments, sample_id="segment"):
    seg_list = []
    for chrom in sorted(segments.keys()):  # Sort the keys lexicographically
        for seg in segments[chrom]:
            seg_list.append((sample_id, chrom, seg[0], seg[1], seg[2], 2))
    res_df = pd.DataFrame(seg_list, columns=["sample_id", "chrom", "start", "end", "name", "cn"])
    return res_df


def chrom_to_sortable(chrom, aut_count = 22):
    if chrom == "chrX":
        return aut_count + 1  # Make 'chrX' sort last
    if chrom == "chrY":
        return aut_count + 2
    if chrom == "chrM":
        return aut_count + 3
    else:
        return int(chrom[3:])  # Remove 'chr' and convert to int
    

def sortable_to_chrom(sortable, aut_count = 22):
    if sortable <= aut_count:
        return "chr" + str(sortable)
    if sortable == aut_count + 1:
        return "chrX"
    if sortable == aut_count + 2:
        return "chrY"
    if sortable == aut_count + 3:
        return "chrM"


def tuples_to_segments(tuples):
    segs = {}
    if len(tuples) > 0 and len(tuples[0]) >= 3:
        for i, tuple in enumerate(tuples):
            if tuple[0] not in segs:
                segs[tuple[0]] = []
            seg_name = tuple[3] if len(tuple) >= 4 else f"seg_{i}"
            segs[tuple[0]].append((tuple[1], tuple[2], seg_name))
    return segs


def breaks_to_segments(breakpoints):
    segs = {}
    for chrom, breaks in breakpoints.items():
        segs[chrom] = []
        for i in range(len(breaks) - 1):
            segs[chrom].append((breaks[i], breaks[i + 1], f"{chrom}_{i}"))
    return segs


def segments_to_breaks(segments):
    breaks = { chrom: [] for chrom in segments }
    for chrom, segs in segments.items():
        for start, end, _ in segs:
            breaks[chrom].append(start)
            breaks[chrom].append(end)
    for chrom in breaks:
        breaks[chrom] = sorted(set(breaks[chrom]))
    return breaks


def genome_to_segments(assembly=hg19):
    segs = {}
    for chrom, len in assembly.chr_lens.items():
        segs[chrom] = [(0, len, chrom)]
    return segs


def bins_to_features(cns_df, cn_columns=None, drop_sex=True):
    cn_columns = get_cn_cols(cns_df, cn_columns)
    sel_df = only_aut(cns_df, inplace=False) if drop_sex else cns_df
    groups = sel_df.groupby("sample_id")
    columns_df = next(iter(groups))[1][['chrom', 'start', 'end', 'name']].set_index('name')
    rows_list = list(groups.groups.keys())

    # Calculate the cumulative count for each group without assigning it to the DataFrame
    cumcount = groups.cumcount()

    if len(cumcount) != len(columns_df) * len(rows_list):
        raise ValueError("The number of cumulative counts does not match the number of rows and columns. Make sure that each sample has the same number of bins.")

    arrays = []
    for cn_col in cn_columns:
        # Use the cumulative count directly in the pivot operation
        array = sel_df.pivot_table(index="sample_id", columns=cumcount, values=cn_col)
        arrays.append(array)

    stacked = np.stack(arrays, axis=0)
    # Transpose from (channels, samples, features) to (samples, channels, features)
    stacked = np.transpose(stacked, (1, 0, 2))
    # if there is only one channel, remove the channel dimension
    if stacked.shape[1] == 1:
        stacked = stacked.squeeze(1)

    return stacked, rows_list, columns_df


def values_count(values_dict):
    return sum(len(values) for values in values_dict.values())