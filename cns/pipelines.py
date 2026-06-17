"""
This module provides functions for processing Copy Number Segment (CNS) data. It includes functions to fill gaps, infer missing values, aggregate data, calculate coverage, and more.

Functions
---------
- main_align: Fills gaps in the CNS data and adds missing chromosomes.
- main_infer: Infers missing values in the CNS data.
- main_impute: Fills gaps, adds missing chromosomes, and infers missing values in the CNS data.
- main_breakage: Identifies breakpoints in CNS data.
- main_coverage: Calculates coverage statistics for CNS data.
- main_ploidy: Calculates ploidy statistics for CNS data.
- main_segment: Segments CNS data based on specified parameters.
- main_aggregate: Aggregates CNS data over specified genomic segments.
- main_seg_agg: Segments CNS data and aggregates the results.
"""

from .analyze import *
from .process import *
from .utils import *


def main_align(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    """
    Aligns all samples with the reference assembly. Adds missing segments, chromosomes, and cuts off the ends if needed.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS (Copy Number Segment) data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. If None, samples are created from `cns_df`.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    assembly : object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with filled gaps and added missing chromosomes.

    Notes
    -----
    This function performs the following steps:
    1. Adds tails to the CNS data to cover chromosome ends.
    2. Fills gaps between segments.
    3. Optionally adds missing chromosomes.
    4. Removes outlier segments.
    5. Merges neighboring segments with the same copy number.

    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    if samples_df is None:
        log_info("No samples provided, creating samples from CNS data.", suppress=not print_info)
        samples_df = samples_df_from_cns_df(cns_df)
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    cn_columns = get_cn_cols(cns_df, cn_columns)
    cns_tailed_df = add_tails(cns_df, assembly=assembly, print_info=print_info)
    cns_aligned_df = fill_gaps(cns_tailed_df, print_info=print_info)
    cns_aligned_df = add_missing(cns_aligned_df, samples_df, assembly=assembly, print_info=print_info)
    cns_cleared_df = remove_outliers(cns_aligned_df, assembly=assembly, print_info=print_info)
    res_df = merge_cns_df(cns_cleared_df, cn_columns, print_info=print_info)
    if segs is not None:
        res_df = aggregate_by_segments(res_df, segs, "none", cn_columns, print_info=print_info)
    return res_df


def main_infer(cns_df, samples_df=None, cn_columns=None, segs=None, method="extend", print_info=False):
    """
    Infers values to replace the NaNs in the CNS data. 

    NOTE: Only replaces NaNs! Usually requires main_align to be ran first.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. Required if `method` is "diploid".
    method : str, optional
        Inference method to use. Options are "extend", "diploid", or "zero". Default is "extend".
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with imputed copy number values.

    Notes
    -----
    This function performs the following steps:
    1. Infers missing values based on the specified method.
    2. Fills any remaining NaNs with zeros.
    3. Merges neighboring segments with the same copy number.

    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    cn_columns = get_cn_cols(cns_df, cn_columns)
    if samples_df is None:
        if method == "diploid":
            log_info("Diploid inference method requires samples_df, but none provided, creating samples from CNS data.", suppress=not print_info)
            samples_df = samples_df_from_cns_df(cns_df)    
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    if segs is not None:
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info=print_info)
    inferred_df = cns_infer(cns_df, samples_df, method, cn_columns=cn_columns, print_info=print_info)
    complete_df = fill_nans_with_zeros(inferred_df, cn_columns=cn_columns, print_info=print_info)
    res_df = merge_cns_df(complete_df, cn_columns=cn_columns, print_info=print_info)
    return res_df


def main_impute(
    cns_df,
    samples_df=None,
    cn_columns=None,
    segs=None,
    method="extend",
    assembly=hg19,
    print_info=False,
):
    """
    Fills gaps in the CNS data, adds missing chromosomes, and infers missing values.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS (Copy Number Segment) data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. If None, samples are created from `cns_df`.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    assembly : object, optional
        Genome assembly to use. Default is `hg19`.
    method : str, optional
        Inference method to use. Options are "extend", "diploid", or "zero". Default is "extend".
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with filled gaps, added missing chromosomes, and inferd values.
    """
    res_df = main_align(cns_df, samples_df, cn_columns, segs, assembly, print_info)
    res_df = main_infer(res_df, samples_df, cn_columns, segs, method, print_info)
    return res_df


# any: if True, based is considered as covered if any CN column has values assigned
def main_coverage(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    """
    Calculates coverage statistics for CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. If None, samples are created from `cns_df`.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    segs : segments dictionary, optional
        Dictionary of segments used for selective masking. Default is None.
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing coverage statistics for each sample.

    Notes
    -----
    The function calculates coverage metrics such as the fraction of the genome covered by CNS data.

    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    if samples_df is None:
        log_info("No samples provided, creating samples from CNS data.", suppress=not print_info)
        samples_df = samples_df_from_cns_df(cns_df)    
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    
    res_df = samples_df.copy()
    cn_columns = get_cn_cols(cns_df, cn_columns)

    if segs is not None:
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)
    norm_sizes = get_norm_sizes(segs, assembly)

    # Select the rows where copy-numbers are not Not a Number (NaN == NaN) is false
    any_nan_df = cn_not_nan(cns_df, cn_columns, True)
    res_df = get_missing_chroms(any_nan_df, res_df, segs, assembly)

    res_df = get_covered_bases(any_nan_df, res_df, True)
    res_df = normalize_feature(res_df, "cover_any", norm_sizes)

    if len(cn_columns) == 2:
        both_nan_df = cn_not_nan(cns_df, cn_columns, False)
        res_df = get_covered_bases(both_nan_df, res_df, False)
        res_df = normalize_feature(res_df, "cover_both", norm_sizes)
    return res_df


def main_breakage(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    """
    Identifies breakpoints in CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    threshold : float, optional
        Threshold for detecting breakpoints. Default is 0.5.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.        
    segs : segments dictionary, optional
        Dictionary of segments used for selective masking. Default is None.    
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing breakpoint information.

    Notes
    -----
    This function detects breakpoints in the CNS data based on changes in copy number values.

    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    if samples_df is None:
        log_info("No samples provided, creating samples from CNS data.", suppress=not print_info)
        samples_df = samples_df_from_cns_df(cns_df)    
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    
    res_df = samples_df.copy()
    cn_columns = get_cn_cols(cns_df, cn_columns)
    if segs is not None:
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)

    # check if non of the cn_columns are NaN
    if cns_df[cn_columns].isna().any().any():
        raise RuntimeError("Cannot calculate breakage with NaN values in CN columns, infer values first.")

    total_added = False
    if len(cn_columns) == 2:
        cns_df["total_cn"] = cns_df[cn_columns].sum(axis=1)
        cn_columns.append("total_cn")
        total_added = True

    for cn_col in cn_columns:
        cns_subset_df = cns_df[["sample_id", "chrom", "start", "end", cn_col]]
        segs_df = merge_cns_df(cns_subset_df, cn_col, False)
        res_df = calc_breaks_per_sample(segs_df, res_df, cn_col, assembly)
        res_df = calc_step_per_sample(segs_df, res_df, cn_col, assembly)

    if total_added:
        cns_df.drop(columns=["total_cn"], inplace=True)

    return res_df


def main_ploidy(cns_df, samples_df=None, cn_columns=None, segs=None, assembly=hg19, print_info=False):
    """
    Calculates ploidy statistics for CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame, optional
        DataFrame containing sample information. If None, samples are created from `cns_df`.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.   
    segs : segments dictionary, optional
        Dictionary of segments used for selective masking. Default is None.    
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing ploidy statistics for each sample.

    Notes
    -----
    This function calculates ploidy metrics such as the fraction of the genome that is aneuploid.
    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    cn_columns = get_cn_cols(cns_df, cn_columns)
    cns_df["start"] = cns_df["start"].astype(np.int64)
    cns_df["end"] = cns_df["end"].astype(np.int64)

    if samples_df is None:
        log_info("No samples provided, creating samples from CNS data.", suppress=not print_info)
        samples_df = samples_df_from_cns_df(cns_df)    
    elif not isinstance(samples_df, pd.DataFrame):
        raise ValueError(f"samples_df must be a DataFrame, got {type(samples_df)}")
    res_df = samples_df.copy()
    
    if segs is not None:
        log_info("Aggregating CN data by provided segments.", suppress=not print_info)
        cns_df = aggregate_by_segments(cns_df, segs, "none", cn_columns, print_info)

    norm_sizes = get_norm_sizes(segs, assembly)

    if cns_df[cn_columns].isna().any().any():
        log_warn("NaNs are not considered in ploidy calculations, it is recommended to infer values first.", suppress=not print_info)
        cns_df = cns_df[cns_df[cn_columns].notna().all(axis=1)]

    log_info("Calculating LOH for each sample.", suppress=not print_info)
    res_df = calc_loh_bases(res_df, cns_df, cn_columns, "both", assembly)
    res_df = normalize_feature(res_df, "loh_both", norm_sizes)
    res_df = calc_loh_bases(res_df, cns_df, cn_columns, "any", assembly)
    res_df = normalize_feature(res_df, "loh_any", norm_sizes)
    log_info("Calculating aneuploidy for each sample.", suppress=not print_info)
    res_df = calc_ane_bases(res_df, cns_df, cn_columns, "any", assembly)
    res_df = normalize_feature(res_df, "ane_any", norm_sizes)
    if len(cn_columns) == 2:
        res_df = calc_ane_bases(res_df, cns_df, cn_columns, "both", assembly)
        res_df = normalize_feature(res_df, "ane_both", norm_sizes)
        log_info("Calculating imbalance for each sample.", suppress=not print_info)
        for col_i in range(2):
            res_df = calc_imb_bases(cns_df, res_df, cn_columns, col_index=col_i, assembly=assembly)
            res_df = normalize_feature(res_df, f"imb_{cn_columns[col_i]}", norm_sizes)

    log_info("Calculating ploidy for each sample.", suppress=not print_info)
    for cn_col in cn_columns:
        res_df[f"ploidy_{cn_col}"] = calc_ploidy_per_column(cns_df, cn_col)
    if len(cn_columns) == 2:
        res_df["ploidy_total_cn"] = res_df[f"ploidy_{cn_columns[0]}"] + res_df[f"ploidy_{cn_columns[1]}"]
    return res_df


def main_segment(
    select_segs=None,
    remove_segs=None,
    split_size=-1,
    merge_dist=-1,
    keep_ends=True,
    filter_size=-1,
    pad_size=0,
    align_to_assembly=False,
    assembly=hg19,
    print_info=False,
):
    """
    Creates a segmentation based on specific segments.

    Parameters
    ----------
    select_segs : 
        Segments to select for computation. By default covers the whole assembly.
    remove_segs : segment dictionary, optional
        Segments to remove from the selection. Be default nothing is removed.
    split_size : int, optional
        Size in base pairs to split segments. Default is -1 (no splitting).
    merge_dist : int, optional
        Distance in base pairs to merge nearby segments. Default is -1 (no merging). If 0, merges all touching segments.
    keep_ends : bool, optional
        If True, clustering (merge_dist > 0) will not cluster start and end breakpoint of each chromosomes.
    filter_size : int, optional
        Minimum size in base pairs to filter segments. Default is -1 (no filtering).
    align_to_assembly : bool, optional
        If True, aligns segments to the assembly. Default is False.
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    dictionary of segments
        Dictionary of segments after processing.

    """
    if select_segs is None:
        select_segs = genome_to_segments(assembly)
    elif not isinstance(select_segs, dict):
        raise ValueError(f"input_segs must a dictionary of segments, got {type(remove_segs)}")
    if filter_size > 0:
        select_segs = filter_cons_size(select_segs, filter_size)
    if remove_segs != None:
        if not isinstance(remove_segs, dict):
            raise ValueError(f"remove_segs must be None or a dictionary of segments, got {type(remove_segs)}")
        if pad_size > 0:
            remove_segs = pad_segments(remove_segs, pad_size, assembly)
            remove_segs = merge_segments(remove_segs, True)
        if filter_size > 0:
            remove_segs = filter_cons_size(remove_segs, filter_size)
        select_segs = segment_difference(select_segs, remove_segs)
        if filter_size > 0:
            select_segs = filter_cons_size(select_segs, filter_size)
    if merge_dist > 0:
        select_segs = cluster_segments(select_segs, merge_dist, keep_ends, print_info)
    if merge_dist == 0:
        select_segs = merge_segments(select_segs)
    if split_size > 0:
        select_segs = split_segments(select_segs, split_size)
    if align_to_assembly:
        select_segs = align_segs_to_assembly(select_segs, sorted=True, assembly=assembly)
    return select_segs


def main_aggregate(cns_df, segs, how="mean", cn_columns=None, print_info=False):
    """
    Aggregates CNS data over specified genomic segments.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    segs : pandas.DataFrame
        DataFrame containing segments over which to aggregate CNS data.
    how : str, optional
        Aggregation method. Options are "mean", "min", "max", "round", or "none". Default is "mean".
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with aggregated copy number values.

    Notes
    -----
    If `how` is not "none" and there are NaNs in `cns_df`, a warning is issued because NaNs are not considered in aggregation.

    Examples
    --------
    >>> aggregated_cns = main_aggregate(cns_df, segs, how="max")
    """
    if not isinstance(cns_df, pd.DataFrame):       
        raise ValueError(f"cns_df must be a DataFrame, got {type(cns_df)}") 
    cn_columns = get_cn_cols(cns_df, cn_columns)
    if how not in ["", "none"] and cns_df[cn_columns].isna().any().any():
        log_warn("NaNs found, it is recommended to infer values first.", suppress=not print_info)
    return aggregate_by_segments(cns_df, segs, how, cn_columns, print_info)


def main_seg_agg(
    cns_df,
    select_segs=None,
    remove_segs=None,
    how="mean",
    split_size=-1,
    merge_dist=-1,
    filter_size=-1,
    pad_size=0,
    cn_columns=None,
    assembly=hg19,
    print_info=False,
):
    """
    Segments CNS data and aggregates the results.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS (Copy Number Segment) data.
    select_segs : segment dictionary, optional
        Segments to select for computation. By default covers the whole assembly.
    remove_segs : segments dictionary, optional
        Segments to remove from the selection.
    how : str, optional
        Aggregation method to use. Default is "mean".
    split_size : int, optional
        Size in base pairs to split segments. Default is -1 (no splitting).
    merge_dist : int, optional
        Distance in base pairs to merge nearby segments. Default is -1 (no merging).
    filter_size : int, optional
        Minimum size in base pairs to filter segments. Default is -1 (no filtering).
    pad_size : int, optional
        Size in base pairs to pad segments on both sizes. Default is 0 (no padding).
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from `cns_df`.
    assembly : Assembly object, optional
        Genome assembly to use. Default is `hg19`.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with aggregated CNS data.
    """
    segs = main_segment(select_segs, remove_segs, split_size, merge_dist, filter_size, pad_size, assembly, print_info)
    res_df = main_aggregate(cns_df, segs, how, cn_columns, print_info)
    return res_df