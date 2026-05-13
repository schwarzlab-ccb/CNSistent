import numpy as np
import pandas as pd
from numba import njit
from cns.process.breakpoints import make_breaks
from cns.utils.conversions import breaks_to_segments, calc_mid, calc_cum_mid
from cns.utils.canonization import get_cn_cols
from cns.utils.assemblies import hg19
from cns.utils.logging import log_info, log_warn


@njit
def _mean_func(cns_array, weights):
    return [np.average(cns_array[:, i], weights=weights) for i in range(cns_array.shape[1])]


@njit
def _max_func(cns_array, weights):
    return [np.max(cns_array[:, i]) for i in range(cns_array.shape[1])]


@njit
def _min_func(cns_array, weights):
    return [np.min(cns_array[:, i]) for i in range(cns_array.shape[1])]


def _aggregate_regs(sample_id, chrom, values, seg_start, seg_end, seg_name, agg_func):
    row_id = 0
    seg_cns = []
    weights = []
    cns_cols_count = values.shape[1] - 2
    while row_id < len(values) and values[row_id, 0] < seg_end:
        if values[row_id, 1] > seg_start:
            row = values[row_id]
            start = max(row[0], seg_start)
            end = min(row[1], seg_end)
            weight = end - start
            if weight > 0:
                seg_cns.append(row[2:])
                weights.append(weight)
        row_id += 1
    
    if seg_cns == []:
        # insert NaN when no data found
        return [sample_id, chrom, seg_start, seg_end] + [np.nan] * cns_cols_count + [seg_name]
    sel_array = np.array(seg_cns)
    weight_array = np.array(weights, dtype=np.uint32)
    cns = agg_func(sel_array, weight_array)
    return [sample_id, chrom, seg_start, seg_end] + cns + [seg_name]


def _mask_by_regs(sample_id, chrom, values, seg_start, seg_end, seg_name):
    row_id = 0
    seg_cns = []
    while row_id < len(values) and values[row_id, 0] < seg_end:
        if values[row_id, 1] > seg_start:
            row = values[row_id]
            start = max(row[0], seg_start)
            end = min(row[1], seg_end)
            seg_cns.append([sample_id, chrom, start, end] + list(row[2:]) + [seg_name])
        row_id += 1

    return seg_cns


def _get_agg_func(how):
    if how == "" or how is None or how == "none":
        return None
    if how == "mean":
        return _mean_func
    if how == "max":
        return _max_func
    if how == "min":
        return _min_func
    raise ValueError("how must be one of ['mean', 'max', 'min', 'none', '']  got " + how)


# Add column names
def aggregate_by_segments(cns_df, segs, how="mean", cn_columns=None, print_info=True):
    """
    Aggregates CNS data by specified segments.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    segments : pandas.DataFrame
        DataFrame containing segments to use for aggregation.
    how : str, optional
        Aggregation method to use. Default is "mean".
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is True.

    Returns
    -------
    pandas.DataFrame
        DataFrame with aggregated CNS data.
    """
    agg_func = _get_agg_func(how)
    cn_columns = get_cn_cols(cns_df, cn_columns)
    sel_cols = ["sample_id", "chrom", "start", "end"] + cn_columns
    cns_df_view = cns_df[sel_cols].set_index(["sample_id", "chrom"])
    new_rows = []
    indices = cns_df_view.index.unique()
    i = 0
    for i, ((sample, chrom), group) in enumerate(cns_df_view.groupby(level=[0, 1])):
        if print_info:
            print(f"Aggregating chr ({i+1}/{len(indices)})", end="\r")
        if chrom not in segs:
            continue
        for seg_start, seg_end, seg_name in segs[chrom]:
            if agg_func != None:
                cn_segs = _aggregate_regs(sample, chrom, group.values, seg_start, seg_end, seg_name, agg_func)
                new_rows.append(cn_segs)
            else:
                cn_segs = _mask_by_regs(sample, chrom, group.values, seg_start, seg_end, seg_name)
                new_rows.extend(cn_segs)
    if print_info:
        print(f"Aggregation finished. Converting {len(new_rows)} rows...", end="\r")
    res_df = pd.DataFrame(new_rows, columns=sel_cols + ["name"]) 
    res_df["start"] = res_df["start"].astype(np.uint32)
    res_df["end"] = res_df["end"].astype(np.uint32)
    log_info(f"Aggregated into {len(new_rows)} CNS." + " " * 40, suppress=not print_info)
    return res_df


def aggregate_by_breaks(cns_df, breaks, how="mean", cn_columns=None, print_info=True):
    """
    Aggregates CNS data by specified breaks.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    breaks : list of tuples
        List of breakpoints to use for aggregation.
    how : str, optional
        Aggregation method to use. Default is "mean".
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is True.

    Returns
    -------
    pandas.DataFrame
        DataFrame with aggregated CNS data.
    """
    segments = breaks_to_segments(breaks)
    return aggregate_by_segments(cns_df, segments, how, cn_columns, print_info)


def aggregate_by_break_type(cns_df, break_type, assembly=hg19, how="mean", cn_columns=None, print_info=True):
    """
    Aggregates CNS data by break type.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    break_type : str
        Type of break to use for aggregation.
    assembly : object, optional
        Genome assembly to use. Default is hg19.
    how : str, optional
        Aggregation method to use. Default is "mean".
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is True.

    Returns
    -------
    pandas.DataFrame
        DataFrame with aggregated CNS data.
    """
    breaks = make_breaks(break_type, assembly=assembly)
    return aggregate_by_breaks(cns_df, breaks, how, cn_columns, print_info)


def add_total_cn(cns_df, cn_columns=None, remove_cn_columns=False, inplace=True):
    """
    Adds a total copy number (CN) column to the CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    remove_cn_columns : bool, optional
        If True, removes the individual CN columns after adding the total CN column. Default is False.
    inplace : bool, optional
        If True, modifies the original DataFrame. Default is True.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the added total CN column.
    """
    cn_columns = get_cn_cols(cns_df, cn_columns)
    # remove total_cn from cn_columns if it is there
    if "total_cn" in cn_columns:
        cn_columns.remove("total_cn")
    if not inplace:
        cns_df = cns_df.copy()
    cns_df["total_cn"] = cns_df[cn_columns].sum(axis=1)
    if remove_cn_columns:
        cns_df.drop(columns=cn_columns, inplace=True)
    return cns_df


def add_mid(cns_df, inplace=True):
    """
    Adds a 'mid' column to the CNS data representing the middle point of each segment.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data with 'start' and 'end' columns.
    inplace : bool, optional
        If True, modifies the original DataFrame. Default is True.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the added 'mid' column.
    """
    if not inplace:
        cns_df = cns_df.copy()
    cns_df["mid"] = calc_mid(cns_df)
    return cns_df


def add_cum_mid(cns_df, assembly=hg19, inplace=True):
    """
    Adds a 'cum_mid' column to the CNS data representing the cumulative middle position 
    across chromosomes based on the specified genome assembly.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data with 'chrom', 'start', and 'end' columns.
    assembly : object, optional
        Genome assembly to use for chromosome lengths. Default is hg19.
    inplace : bool, optional
        If True, modifies the original DataFrame. Default is True.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the added 'cum_mid' column.
    """
    if not inplace:
        cns_df = cns_df.copy()
    cns_df["cum_mid"] = calc_cum_mid(cns_df, assembly)
    return cns_df


def group_samples(cns_df, cn_columns=None, how="mean", group_name="grouped"): 
    """
    Groups CNS data by samples and aggregates the results.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    how : str, optional
        Aggregation method to use. Options are "mean", "max", "min". Default is "mean".
    group_name : str, optional
        Name for the grouped samples. Default is "grouped".

    Returns
    -------
    pandas.DataFrame
        DataFrame with grouped and aggregated CNS data.
    """
    if len(cns_df) == 0:
        log_warn("No data to group.")
        return cns_df   
    if how not in ["mean", "max", "min"]:
        raise ValueError("to group samples, how must be one of ['mean', 'max', 'min']")
    cn_columns = get_cn_cols(cns_df, cn_columns)
    grouped = cns_df.drop("sample_id", axis=1).groupby(["chrom", "start", "end"])

    # calculate mean on grouped except for chrom, where take the first value
    agg_scheme = {}
    if "name" in cns_df.columns:
        agg_scheme["name"] = "first"
    for column in cn_columns:
        agg_scheme[column] = how
    grouped = grouped.agg(agg_scheme).reset_index()
    grouped["sample_id"] = group_name
    return grouped


def stack_groups(cns_dfs, labels=None):
    """
    Stacks multiple CNS DataFrames into a single DataFrame.

    Parameters
    ----------
    cns_dfs : list of pandas.DataFrame
        List of CNS DataFrames to stack.
    labels : list of str, optional
        List of labels for the stacked DataFrames. If specified, the length of labels must be equal to the number of DataFrames.

    Returns
    -------
    pandas.DataFrame
        Stacked DataFrame.
    """
    if not isinstance(cns_dfs, list):
        cns_dfs = [cns_dfs]
    if labels is not None:
        if len(cns_dfs) != len(labels):
            raise ValueError("If specified, the length of labels must be equal to the number of dataframes.")
        for i, df in enumerate(cns_dfs):
            df["sample_id"] = labels[i]
    return pd.concat(cns_dfs)


def mean_value_per_seg(cns_df, segs, score_col):
    """
    Calculate weighted scores for segments based on overlap with copy number segments.
    
    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing copy number segments
    segs : dict
        Dictionary mapping chromosome names to lists of segments (start, end, name)
    score_col : str
        Column name containing the scores to aggregate
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with scores for each segment
    """
    scores = []
    
    # Group the data by chromosome once
    cns_by_chrom = {chrom: df.reset_index(drop=True) for chrom, df in cns_df.groupby("chrom")}
    
    for chrom, segs_list in segs.items():
        if chrom not in cns_by_chrom:
            continue
            
        chr_df = cns_by_chrom[chrom]
        
        for s_start, s_end, s_name in segs_list:
            seg_len = s_end - s_start
            
            # Find overlapping segments more efficiently
            overlaps = chr_df[(chr_df["end"] > s_start) & (chr_df["start"] < s_end)]
            
            if len(overlaps) == 0:
                scores.append([chrom, s_start, s_end, s_name, 0])
                continue
                
            # Calculate weighted score in one step
            overlap_sizes = (
                np.minimum(overlaps["end"], s_end) - 
                np.maximum(overlaps["start"], s_start)
            )
            weighted_score = np.sum(overlap_sizes * overlaps[score_col]) / seg_len
            
            scores.append([chrom, s_start, s_end, s_name, weighted_score])
            
    return pd.DataFrame(scores, columns=["chrom", "start", "end", "name", score_col])