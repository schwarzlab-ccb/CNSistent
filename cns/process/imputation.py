import numpy as np
import pandas as pd
from cns.utils.canonization import get_cn_cols
from cns.utils.logging import log_info
from cns.utils.assemblies import hg19


def add_tails(cns_df, assembly=hg19, print_info=True):
    """
    Adds tails to the CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    assembly : Assembly object
        Assembly object containing chromosome length information.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with tails added.
    """
    chr_lens = assembly.chr_lens
    grouped = cns_df.groupby(["sample_id", "chrom"]).agg({"start": "min", "end": "max"})
    grouped = grouped.rename(columns={"start": "min_start", "end": "max_end"})
    grouped = grouped.reset_index()
    missing_ranges = []
    for _, row in grouped.iterrows():
        if row.min_start > 0:
            missing_ranges.append(
                {
                    "sample_id": row.sample_id,
                    "chrom": row.chrom,
                    "start": 0,
                    "end": row.min_start
                }
            )
        if row.max_end < chr_lens[f"{row.chrom}"]:
            missing_ranges.append(
                    {
                        "sample_id": row.sample_id,
                        "chrom":  row.chrom,
                        "start": row.max_end,
                        "end": chr_lens[str(row.chrom)]
                    }
            )

    if len(missing_ranges) == 0:
        log_info(f"No missing ends found.", suppress=not print_info)
        return cns_df.copy()
    else:
        log_info(f"Adding {len(missing_ranges)} missing ends", suppress=not print_info)
        new_cns_df = pd.DataFrame(missing_ranges, columns=cns_df.columns)
        res_df = pd.concat([cns_df, new_cns_df])
        res_df.sort_values(
            by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True
        )
        return res_df


def fill_gaps(cns_df, print_info=True):
    """
    Fills gaps in the CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with gaps filled.
    """
    # Iterate over the rows
    new_rows = []
    for i in range(len(cns_df) - 1):
        # Check if the next row has the same 'sample_id' and 'chrom'
        if (
            cns_df.at[i, "sample_id"] == cns_df.at[i + 1, "sample_id"]
            and cns_df.at[i, "chrom"] == cns_df.at[i + 1, "chrom"]
        ):
            # Calculate the range
            range_start = cns_df.at[i, "end"]
            range_end = cns_df.at[i + 1, "start"]
            # If the range is greater than 1, add a new row
            if range_end > range_start:
                new_rows.append(
                    {
                        "sample_id": cns_df.at[i, "sample_id"],
                        "chrom": cns_df.at[i, "chrom"],
                        "start": range_start,
                        "end": range_end,
                    }
                )

    if len(new_rows) == 0:
        log_info(f"No gaps found.", suppress=not print_info)
        return cns_df.copy()
    else:
        # Concatenate the cns_dfs
        log_info(f"Filling {len(new_rows)} gaps.", suppress=not print_info)
        new_cns_df = pd.DataFrame(new_rows, columns=cns_df.columns)
        res_df = pd.concat([cns_df, new_cns_df])
        res_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
        return res_df


# Add fully missing chromosomes
def add_missing(cns_df, samples_df=None, assembly=hg19, print_info=True):
    """
    Adds missing chromosomes to the CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    chr_lens : dict
        Dictionary of chromosome lengths.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with missing chromosomes added.
    """
    res_df = cns_df.set_index("sample_id")

    new_entries = []
    for sample in res_df.index.unique():
        cns_sample_df = res_df.loc[sample]
        sample_chroms = cns_sample_df["chrom"].values
        for chromosome in assembly.chr_names:
            if chromosome not in sample_chroms:
                if chromosome == assembly.chr_x and samples_df is None:
                    continue
                if chromosome == assembly.chr_y and (samples_df is None or samples_df.loc[sample].sex == "xx"):
                    continue
                new_entry = {
                    "sample_id": sample,
                    "chrom": chromosome,
                    "start": 0,
                    "end": assembly.chr_lens[chromosome],
                }
                new_entries.append(new_entry)

    if len(new_entries) == 0:
        log_info(f"No missing chromosomes found.", suppress=not print_info)
        return cns_df.copy()
    else:
        log_info(f"Adding {len(new_entries)} missing chromosomes.", suppress=not print_info)
        empty_chrs_df = pd.DataFrame(new_entries)
        res_df.reset_index(inplace=True)
        res_df = pd.concat([res_df, empty_chrs_df])
        res_df.sort_values(
            by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True
        )
        return res_df


def remove_outliers(cns_df, assembly=hg19, print_info=True):
    """
    Removes outliers from the CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    assembly : Assembly object
        Assembly object containing chromosome length information.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with outliers removed.
    """
    res_df = cns_df.copy()
    chr_lens = assembly.chr_lens
    idx_to_remove = []
    for i, (index, row) in enumerate(res_df.iterrows()):
        if row.start < 0:
            if row.end < 0:
                idx_to_remove.append(i)
                log_info(f"Removed outlier:\n{row}", suppress=not print_info)
            else:
                res_df.at[index, "start"] = 0
        if row.end > chr_lens[row.chrom]:
            if row.start >= chr_lens[row.chrom]:
                idx_to_remove.append(i)
                log_info(f"Removed outlier:\n{row}", suppress=not print_info)
            else:
                res_df.at[index, "end"] = chr_lens[row.chrom]

    log_info(f"Removed outliers: {len(idx_to_remove)}", suppress=not print_info)
    # remove from cns_df where idx_to_remove is in the index
    res_df = res_df.drop(res_df.index[idx_to_remove]).sort_values(by=["sample_id", "chrom", "start"], ignore_index=True)
    return res_df

# Makes sure that the columns are of the correct type
def _are_mergeable(a, b, cn_columns):
    return (
        a.sample_id == b.sample_id
        and a.chrom == b.chrom
        and a.end == b.start
        and all([(a[col] == b[col]) or (np.isnan(a[col]) and np.isnan(b[col])) for col in cn_columns])
    )


def merge_cns_df(cns_df, cn_columns=None, print_info=True):
    """
    Merges consecutive CNS segments with the same copy number values.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with merged CNS segments.
    """
    cn_columns = get_cn_cols(cns_df, cn_columns)    
    res_df = cns_df.copy()
    idx_to_remove = []
    use_name = "name" in res_df.columns

    for i, (index, row) in enumerate(res_df.iterrows()):
        if i != 0 and _are_mergeable(prev, row, cn_columns):
            idx_to_remove.append(i - 1)
            res_df.at[index, "start"] = prev.start
            row.start = prev.start  # update the comparison copy too
            if use_name:
                row.name = str(prev.name)  + "&" + str(row.name)  # update the name for logging
        prev = row

    log_info(f"Merged entries: {len(idx_to_remove)}", suppress=not print_info)

    # remove from cns_df where idx_to_remove is in the index
    res_df = res_df.drop(res_df.index[idx_to_remove]).sort_values(by=["sample_id", "chrom", "start"], ignore_index=True)
    return res_df


def _is_same_contig(df, id, chrom, j):
    return df.at[j, "sample_id"] == id and df.at[j, "chrom"] == chrom

 
def _impute_extend(cns_df, cn_columns, print_info=True):
    """
    For each column with NaN, find the previous and next existing value and fill up to midpoint. 
    If only one exists, fill up to the next or previous value. If none exist, fill with 0.
    """
    new_entries = []
    new_vals = {}
    for i in range(len(cns_df)):
        if any(np.isnan(cns_df.at[i, col]) for col in cn_columns):
            
            id = cns_df.at[i, "sample_id"]
            chrom = cns_df.at[i, "chrom"]            
            start = cns_df.at[i, "start"]
            end = cns_df.at[i, "end"]

            for col in cn_columns:
                if not np.isnan(cns_df.at[i, col]):
                    new_vals[col] = cns_df.at[i, col]
                    continue
                prev_idx = i - 1
                next_idx = i + 1
                while prev_idx >= 0 and np.isnan(cns_df.at[prev_idx, col]) and _is_same_contig(cns_df, id, chrom, prev_idx) and cns_df.at[prev_idx, "end"] == cns_df.at[prev_idx + 1, "start"]:
                    prev_idx -= 1
                if prev_idx < 0 or np.isnan(cns_df.at[prev_idx, col]) or not _is_same_contig(cns_df, id, chrom, prev_idx) or cns_df.at[prev_idx, "end"] != cns_df.at[prev_idx + 1, "start"]:
                    prev_idx = -1
                while next_idx < len(cns_df) and np.isnan(cns_df.at[next_idx, col]) and _is_same_contig(cns_df, id, chrom, next_idx) and cns_df.at[next_idx - 1, "end"] == cns_df.at[next_idx, "start"]:
                    next_idx += 1
                if next_idx >= len(cns_df) or np.isnan(cns_df.at[next_idx, col]) or not _is_same_contig(cns_df, id, chrom, next_idx) or cns_df.at[next_idx - 1, "end"] != cns_df.at[next_idx, "start"]:
                    next_idx = -1
                if prev_idx == -1 and next_idx == -1:
                    new_vals[col] = 0
                elif prev_idx == -1 and next_idx != -1:
                    new_vals[col] = cns_df.at[next_idx, col]
                elif prev_idx != -1 and next_idx == -1:
                    new_vals[col] = cns_df.at[prev_idx, col]
                else:
                    prev_end = cns_df.at[prev_idx, "end"]
                    next_start = cns_df.at[next_idx, "start"]
                    midpoint = prev_end + (next_start - prev_end) // 2
                    if midpoint <= start:
                        new_vals[col] = cns_df.at[next_idx, col]
                    elif midpoint >= end:
                        new_vals[col] = cns_df.at[prev_idx, col]
                    else:
                        new_vals[col] = (midpoint, cns_df.at[prev_idx, col], cns_df.at[next_idx, col])
            mid_count = 0
            for k, v in new_vals.items():
                if isinstance(v, tuple):
                    mid_count += 1
                    break
            if mid_count == 0:
                new_simple = [id, chrom, start, end] + [new_vals[col] for col in cn_columns]
                new_entries.append(new_simple)
            elif mid_count == 1 or new_vals[cn_columns[0]][0] == new_vals[cn_columns[1]][0]:
                midpoint = new_vals[cn_columns[0]][0] if isinstance(new_vals[cn_columns[0]], tuple) else new_vals[cn_columns[1]][0]
                first_half = [id, chrom, start, midpoint] + [new_vals[col][1] if isinstance(new_vals[col], tuple) else new_vals[col] for col in cn_columns]
                second_half = [id, chrom, midpoint, end] + [new_vals[col][2] if isinstance(new_vals[col], tuple) else new_vals[col] for col in cn_columns]
                new_entries.append(first_half)                    
                new_entries.append(second_half)
            else:
                midpoints = sorted([new_vals[col][0] for col in cn_columns])
                first_part = [id, chrom, start, midpoints[0]] + [new_vals[col][1] for col in cn_columns]
                second_part = [id, chrom, midpoints[0], midpoints[1]] + [new_vals[col][1] if new_vals[col][0] <= midpoints[1] else new_vals[col][2] for col in cn_columns]
                third_part = [id, chrom, midpoints[1], end] + [new_vals[col][2] for col in cn_columns]
                new_entries.append(first_part)
                new_entries.append(second_part)
                new_entries.append(third_part)

    
    new_cols = ["sample_id", "chrom", "start", "end"] + cn_columns
    imputation_df = pd.DataFrame(new_entries, columns=new_cols)
    query = ' or '.join([f"{col}.isnull()" for col in cn_columns])
    idx_to_remove = cns_df.query(query).index

    log_info(f"New entries: {imputation_df.shape[0]}\nRemoved entries: {len(idx_to_remove)}", suppress=not print_info)
    # remove from cns_df where idx_to_remove is in the index
    filtered_df = cns_df.drop(idx_to_remove)
    # concat the new_table to cns_df
    res_df = filtered_df if len(imputation_df) == 0 else pd.concat([filtered_df, imputation_df])
    # sort cns_df by sample_id, chr, start
    res_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    return res_df


def _impute_diploid(cns_df, samples_df, cn_columns, print_info=True):
    if len(cn_columns) > 2 or len(cn_columns) < 1:
        raise ValueError("Diploid imputation can only be done for one (total CN) or two (major, minor CN) columns.")
    
    aut_df = cns_df.query("chrom != 'chrX' and chrom != 'chrY'")
    xx_samples = samples_df.query("sex == 'xx'").index
    xy_samples = samples_df.query("sex == 'xy'").index
    xx_cns_df = cns_df[cns_df["sample_id"].isin(xx_samples)]
    xx_x_chrom_df = xx_cns_df.query("chrom == 'chrX'")
    xx_y_chrom_df = xx_cns_df.query("chrom == 'chrY'")
    xy_cns_df = cns_df[cns_df["sample_id"].isin(xy_samples)]
    xy_x_chrom_df = xy_cns_df.query("chrom == 'chrX'")
    xy_y_chrom_df = xy_cns_df.query("chrom == 'chrY'")
    if len(cn_columns) == 2:
        aut_df[cn_columns] = aut_df[cn_columns].fillna(1)
        xx_x_chrom_df[cn_columns] = xx_x_chrom_df[cn_columns].fillna(1)      
        xx_y_chrom_df[cn_columns] = xx_y_chrom_df[cn_columns].fillna(0)
        xy_x_chrom_df[cn_columns[0]] = xy_x_chrom_df[cn_columns[0]].fillna(1)
        xy_x_chrom_df[cn_columns[1]] = xy_x_chrom_df[cn_columns[1]].fillna(0)
        xy_y_chrom_df[cn_columns[0]] = xy_y_chrom_df[cn_columns[0]].fillna(1)
        xy_y_chrom_df[cn_columns[0]] = xy_y_chrom_df[cn_columns[0]].fillna(0)
    else:
        col = cn_columns[0]
        aut_df[col] = aut_df[col].fillna(2)
        xx_x_chrom_df[col] = xx_x_chrom_df[col].fillna(2)
        xx_y_chrom_df[col] = xx_y_chrom_df[col].fillna(0)
        xy_x_chrom_df[col] = xy_x_chrom_df[col].fillna(1)
        xy_y_chrom_df[col] = xy_y_chrom_df[col].fillna(1)

        # Update the original dataframe in place
    cns_df.update(aut_df)
    cns_df.update(xx_x_chrom_df)
    cns_df.update(xx_y_chrom_df)
    cns_df.update(xy_x_chrom_df)
    cns_df.update(xy_y_chrom_df)

    return cns_df


def cns_infer(cns_df, samples_df, method='extend', cn_columns=None, print_info=True):
    """
    Infers NaN values in the CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    samples_df : pandas.DataFrame
        DataFrame containing sample information.
    method : str, optional
        Inference method to use. Options are "extend", "diploid", or "zero". Default is "extend".
        Note - if "name" is present in cns_df, it will be removed when using "extend" method.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with imputed copy number values.
    """    
    cn_columns = get_cn_cols(cns_df, cn_columns)
    if method ==  'extend':
        if 'name' in cns_df.columns:
            cns_df = cns_df.drop(columns=['name'])
        return _impute_extend(cns_df, cn_columns, print_info)
    if method == 'diploid':
        return _impute_diploid(cns_df, samples_df, cn_columns, print_info)
    if method == 'zero':
        return fill_nans_with_zeros(cns_df, cn_columns, print_info)
    else:
        msg = f"Unknown imputation method: {method}"
        raise Exception(msg)


def fill_nans_with_zeros(cns_df, cn_columns=None, print_info=True):
    """
    Fills NaN values in the CNS data with zeros.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing CNS data.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from cns_df.
    print_info : bool, optional
        If True, prints informational messages during processing. Default is False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with NaN values filled with zeros.
    """   
    cn_columns = get_cn_cols(cns_df, cn_columns)    
    res_df = cns_df.copy()
    log_info(f"Filling {res_df[cn_columns].isna().any(axis=1).sum()} NaN rows with zero", suppress=not print_info)
    # Fully missing chromosomes filled with 0
    for col in cn_columns:
        res_df[col] = res_df[col].fillna(0.0).infer_objects()
    return res_df