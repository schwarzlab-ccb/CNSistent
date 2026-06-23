"""
Utility functions for loading and managing copy number segment (CNS) data.

This module provides functions to load, filter, and manage CNS data from different sources
(PCAWG, TRACERx, TCGA). It includes utilities for handling samples, loading binned data,
and managing file paths.
"""

import pandas as pd
from os.path import join as pjoin, abspath, dirname
from cns.utils import log_info, select_CNS_samples, load_cns, load_samples, load_segments, find_bends, z_score_filter
import matplotlib.pyplot as plt
import os


def get_root_path():
    """Get the root path of the CNSistent package.

    Returns
    -------
    str
        Absolute path to the package root directory.
    """
    return abspath(pjoin(dirname(__file__), ".."))


img_path = pjoin(get_root_path(), "img")
out_path = pjoin(get_root_path(), "out")
data_path = pjoin(get_root_path(), "data")
docs_path = pjoin(get_root_path(), "docs/source/files")


def load_cns_file(filename, print_info=False):
    """Load CNS data from the output directory.

    Parameters
    ----------
    filename : str
        Name of the file to load from the output directory.
    print_info : bool, optional
        If True, print information about the loading process.

    Returns
    -------
    pd.DataFrame
        DataFrame containing CNS data.
    """
    source_folder = data_path if "raw" in filename else out_path
    path = abspath(pjoin(source_folder, filename))
    log_info(f"Loading CNS data  from {path}", suppress=not print_info)
    cns_df = load_cns(path)
    return cns_df


def load_samples_file(filename, use_filter=False, print_info=False):
    """Load sample data from the output directory.

    Parameters
    ----------
    filename : str
        Name of the sample file to load.
    filter : bool, optional
        Whether to filter samples based on coverage and aneuploidy.
    print_info : bool, optional
        Whether to print progress information.

    Returns
    -------
    pd.DataFrame
        DataFrame containing sample data.
    """
    source_folder = data_path if "raw" in filename else out_path
    path = abspath(pjoin(source_folder, filename))
    log_info(f"Loading samples from {path}", suppress=not print_info)
    samples_df = load_samples(path)

    if use_filter:
        # calculate bend for aneuploidy
        ane_bends = find_bends(samples_df["ane_any_aut"], max_val=0.01)
        ane_min_frac = ane_bends[0][ane_bends[2]] if ane_bends[2] > 0 else 0.01

        # calculate the z-score for coverage
        cover_filtered = z_score_filter(samples_df["cover_any_aut"])
        cover_min_frac = cover_filtered.min()

        samples_df = _filter_samples(samples_df, ane_min_frac, cover_min_frac, print_info)
    else:
        if "TCGA_type" in samples_df.columns:
            samples_df["type"] = samples_df["TCGA_type"]    

    return samples_df


def load_segs_out(filename, print_info=False):
    """Load segment data from the output directory.

    Parameters
    ----------
    filename : str
        Name of the segment file to load.
    print_info : bool, optional
        Whether to print progress information.
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing segment data.
    """
    log_info(f"Loading segments from {filename}", suppress=not print_info)
    return load_segments(pjoin(out_path, filename))



def _filter_samples(samples_df, ane_min_frac=0.001, cover_min_frac=0.95, print_info=False):
    """Filter samples based on aneuploidy and coverage criteria.

    Parameters
    ----------
    samples : pd.DataFrame
        DataFrame containing sample information.
    ane_min_frac : float, optional
        Minimum fraction for aneuploidy filtering.
    cover_min_frac : float, optional
        Minimum fraction for coverage filtering.
    print_info : bool, optional
        Whether to print progress information.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing sample information.
    """
    log_info(f"Total samples: {len(samples_df)}", suppress=not print_info)
    
    cn_neutral = samples_df.query(f"ane_any_aut < {ane_min_frac}").index
    log_info(f"{len(cn_neutral)} samples are CN neutral (below {ane_min_frac:.5f})", suppress=not print_info)

    # Find samples with low coverage (below 95% in autosomes)
    low_coverage = samples_df.query(f"cover_any_aut < {cover_min_frac}").index
    log_info(f"{len(low_coverage)} samples have low coverage (below {cover_min_frac:.5f})", suppress=not print_info)

    if "TCGA_type" in samples_df.columns:
        samples_df["type"] = samples_df["TCGA_type"]    
        samples_df.drop(columns=["TCGA_type"], inplace=True)
    samples_df["type"] = samples_df["type"].replace({"LUADx2": "LUAD"}).replace({"LUADx3": "LUAD"})
    untyped = samples_df[samples_df["type"].fillna('').apply(lambda x: any(not c.isupper() for c in x))].index
    log_info(f"{len(untyped)} samples do not have exact type", suppress=not print_info)

    filtered_df = samples_df.query("(index not in @untyped) & (index not in @cn_neutral) & (index not in @low_coverage)")

    log_info(f"Filtered samples: {len(filtered_df)}", suppress=not print_info)

    return filtered_df.copy()


def main_load(segment_type: str, dataset="all", use_filter=True, concat=True, print_info=False):
    """
    Loads and filters samples from the specified dataset along with the imputed CNS data.

    Parameters
    ----------
    segment_type : str, optional
        Bin size for loading binned data. If None, no CNS data is loaded. Options include: 
        ["1MB", "2MB", "3MB", "5MB", "10MB", "250KB", "500KB", "whole", "arms", "bands", "COSMIC", "ENSEMBL"].
    dataset : str, optional
        Dataset to load. Options include: "PCAWG", "TCGA", "TRACERx", or "all". Default is "all".
    use_filter : bool, optional
        If True, filters samples based on coverage and aneuploidy. Default is True.
    concat : bool, optional
        If True, concatenates samples and CNS data from multiple datasets, otherwise a named dictionary is returned. Default is True.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    tuple
        A tuple containing two pandas DataFrames:
        - samples_df: DataFrame containing sample information and statistics.
        - cns_df: DataFrame containing the CNS data or binned data.
        If concat is False, the DataFrames are returned as a dictionary thereof.

    Examples
    --------
    >>> samples_df, cns_df = main_load("imp")
    >>> samples_df, cns_df = main_load("3MB", "PCAWG")
    """
    if dataset == "all":
        datasets = ["PCAWG", "TRACERx", "TCGA_hg19"]
    else:
        datasets = [dataset]

    samples_dict = {}
    for dataset in datasets:
        if segment_type == "raw":
            samples = load_samples_file(f"{dataset}_samples_raw.tsv", print_info)
        else:
            samples = load_samples_file(f"{dataset}_samples.tsv", use_filter, print_info)
        samples["source"] = dataset
        samples_dict[dataset] = samples
    samples_df = pd.concat(samples_dict.values()) if (concat or len(datasets) == 1) else samples_dict
    log_info(f"Total samples: {len(samples_df)}", suppress=not print_info)

    if segment_type is None:
        log_info("No segment type specified. Returning sample data only.", suppress=not print_info)
        return samples_df, None
        
    file_type = "cns" if segment_type in ["imp", "raw", "align"] else "bin"
    cns_dict = {dataset: load_cns_file(f"{dataset}_{file_type}_{segment_type}.tsv", print_info) for dataset in datasets}
    cns_dict = {k: select_CNS_samples(v, samples_dict[k]).reset_index(drop=True) for k, v in cns_dict.items()}
    cns_df = pd.concat(cns_dict.values()) if (concat or len(datasets) == 1) else cns_dict
    log_info(f"Total CNS segments: {len(cns_df)}", suppress=not print_info)
    
    return samples_df, cns_df


def save_cns_fig(fig_name, fig=None):
    """Save a figure to the images directory.

    Parameters
    ----------
    fig_name : str
        Name of the figure file (without extension).
    fig : matplotlib.figure.Figure, optional
        Figure to save. If None, uses current figure.
    """
    os.makedirs(img_path, exist_ok=True)
    if fig == None:
        fig = plt.gcf()
    fig.savefig(f"{img_path}/{fig_name}.png", bbox_inches="tight", transparent=False, dpi=300)
    fig.savefig(f"{img_path}/{fig_name}.pdf", bbox_inches="tight", transparent=True)


def save_doc_fig(fig_name, fig=None):
    """Save a figure to the documentation directory.

    Parameters
    ----------
    fig_name : str
        Name of the figure file (without extension).
    fig : matplotlib.figure.Figure, optional
        Figure to save. If None, uses current figure.
    """
    if fig == None:
        fig = plt.gcf()
    fig.savefig(f"{docs_path}/{fig_name}.png", bbox_inches="tight", transparent=False, dpi=300)
    fig.savefig(f"{docs_path}/{fig_name}.pdf", bbox_inches="tight", transparent=True)


def load_COSMIC():
    return load_segments(pjoin(data_path, "COSMIC_consensus_genes.bed"))


def load_ENSEMBL():
    return load_segments(pjoin(data_path, "ENSEMBL_coding_genes.bed"))


def samples_above_threshold(samples_df, threshold=50):
    """
    Returns the original values from samples_df where 'type' occurs at least `threshold` times and is not 'Other'.

    Parameters
    ----------
    samples_df : pd.DataFrame
        DataFrame containing a 'type' column.
    threshold : int, optional
        Minimum number of occurrences for a type to be included.

    Returns
    -------
    pd.Series
        Series of 'type' values from the original DataFrame that meet the criteria.
    """
    valid_types = samples_df["type"].value_counts()
    valid_types = valid_types[(valid_types >= threshold) & (valid_types.index != "Other")].index
    return samples_df.loc[samples_df["type"].isin(valid_types)]