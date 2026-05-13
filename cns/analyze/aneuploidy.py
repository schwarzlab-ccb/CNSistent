import numpy as np
import pandas as pd
from cns.utils.assemblies import hg19
from cns.utils.conversions import calc_lengths
from cns.utils.selection import get_chr_sets
from numba import njit
    

def _check_total(chrom, val, sex, feature, allele_spec, chr_x, chr_y):
    if feature == "ane":
        if sex == "xy":
            if chrom == chr_x or chrom == chr_y:
                return val != 1
            else:
                return val != 2
        else:
            if chrom == chr_y:
                return val != 0
            else:    
                return val != 2
    else:
        if sex == "xy":
            if chrom == chr_x or chrom == chr_y:
                return val == 0
            else:
                return val < 2 if allele_spec == "any" else val == 0
        else:
            if chrom == chr_y:
                return False
            else:    
                return val < 2 if allele_spec == "any" else val == 0
            

def _check_alleles(chrom, val_0, val_1, sex, feature, allele_spec, chr_x, chr_y):
    if val_0 < val_1:
        val_0, val_1 = val_1, val_0
    if sex == "xy":
        if chrom == chr_x or chrom == chr_y:
            exp_0 = 1
            exp_1 = 0
        else:
            exp_0 = exp_1 = 1
    else:
        if chrom == chr_y:
            exp_0 = exp_1 = 0
        else:
            exp_0 = exp_1 = 1
    
    if feature == "ane":
        return (exp_0 != val_0 or exp_1 != val_1) if allele_spec == "any" else (exp_0 != val_0 and exp_1 != val_1)
    else:  # check_type == "loh"
        return (val_0 < exp_0 or val_1 < exp_1) if allele_spec == "any" else (val_0 < exp_0 and val_1 < exp_1)
    

def _get_feature_per_seg(cns_df, samples_df, cn_columns, feature, allele_spec, assembly=hg19):
    res = []
    chr_x = assembly.chr_x
    chr_y = assembly.chr_y
    for sample_id, groupd_df in cns_df.groupby("sample_id"):
        sex = samples_df.loc[sample_id]["sex"]
        if len(cn_columns) == 1:
            res.append(groupd_df.apply(lambda row: _check_total(row["chrom"], row[cn_columns[0]], sex, feature, allele_spec, chr_x, chr_y), axis=1).values)
        else:
            res.append(groupd_df.apply(lambda row: _check_alleles(row["chrom"], row[cn_columns[0]], row[cn_columns[1]], sex, feature, allele_spec, chr_x, chr_y), axis=1).values)
    return np.concatenate(res)


def _calc_bases_per_chr_group(res, masked_cns_df, label, groups):
    for suffix, names in groups.items():
        subset = masked_cns_df.query("chrom in @names")
        length = calc_lengths(subset)
        res[f"{label}_{suffix}"] = length.groupby(subset["sample_id"]).sum()
        res[f"{label}_{suffix}"] = res[f"{label}_{suffix}"].infer_objects().fillna(0)
        res[f"{label}_{suffix}"] = res[f"{label}_{suffix}"].astype(np.int64)
    return res


def _count_bases_with_feature(res_df, cns_df, cn_columns, feature, allele_spec, assembly):
    label = feature + "_" + allele_spec
    mask = _get_feature_per_seg(cns_df, res_df, cn_columns, feature, allele_spec, assembly)
    chr_sets = get_chr_sets(cns_df, assembly)
    return _calc_bases_per_chr_group(res_df, cns_df[mask], label, chr_sets)


def calc_loh_bases(samples_df, cns_df, cn_columns, allele_spec, assembly=hg19):
    """
    Calculates the length of Loss of Heterozygosity (LOH) bases for each sample.

    Parameters
    ----------
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str
        List of column names for copy number data.
    allele_spec : str
        Allele specification, either "any" or "both".
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the length of LOH bases for each sample.
    """
    res_df = samples_df.copy()
    return _count_bases_with_feature(res_df, cns_df, cn_columns, "loh", allele_spec, assembly)


def calc_ane_bases(samples_df, cns_df, cn_columns, allele_spec, assembly=hg19):
    """
    Calculates the length of aneuploidy bases for each sample.

    Parameters
    ----------
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str
        List of column names for copy number data.
    allele_spec : str
        Allele specification, either "any" or "both".
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the length of aneuploidy bases for each sample.
    """
    res_df = samples_df.copy()
    return _count_bases_with_feature(res_df, cns_df, cn_columns, "ane", allele_spec, assembly)


def calc_imb_bases(cns_df, samples_df, cn_columns, col_index=0, assembly=hg19):
    """
    Calculates the length of imbalance bases for each sample.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    cn_columns : list of str
        List of column names for copy number data.
    col_index : int, optional
        Index of the column to use for imbalance calculation. Default is 0.
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the length of imbalance bases for each sample.
    """
    res = samples_df.copy()
    if len(cn_columns) != 2:
        raise ValueError("There must be two CN columns to calculate imbalance score")
    cn_col1 = cn_columns[col_index]
    cn_col2 = cn_columns[1 - col_index]
    mask = cns_df[cn_col1] > cns_df[cn_col2]
    label = "imb_" + cn_col1
    chr_sets = get_chr_sets(cns_df, assembly)
    res = _calc_bases_per_chr_group(res, cns_df[mask], label, chr_sets)
    return res


@njit
def _calc_ploidy_per_sample(start, end, cn_values):
    lengths = end - start
    ploidy = (cn_values * lengths).sum() / lengths.sum()
    return ploidy


def calc_ploidy_per_column(cns_df, cn_column):
    """
    Calculates the ploidy for each sample based on a specified CN column.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_column : str
        Column name for copy number data.

    Returns
    -------
    pandas.Series
        Series with the ploidy value for each sample.
    """
    grouped = cns_df.groupby('sample_id')
    res = {}
    for sample_id, group_df in grouped:
        start = group_df["start"].values
        end = group_df["end"].values
        cn_values = group_df[cn_column].values
        ploidy = _calc_ploidy_per_sample(start, end, cn_values)
        res[sample_id] = ploidy
    return pd.Series(res)


def calc_chrom_var(df, cn_column):
    """
    Calculates the variance of a specified column grouped by chromosome.
    
    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_column : str
        Column name for which to calculate the variance.

    Returns
    -------
    pandas.Series
        Series with the variance of the specified column for each chromosome.
    """
    return df.groupby("chrom")[cn_column].var()


def calc_chrom_mean(df, cn_column):
    """
    Calculates the mean of a specified column grouped by chromosome.
    
    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_column : str
        Column name for which to calculate the mean.

    Returns
    -------
    pandas.Series
        Series with the mean of the specified column for each chromosome.
    """
    return df.groupby("chrom")[cn_column].mean()

