from os.path import abspath, exists
import numpy as np
import pandas as pd
from io import StringIO

from cns.process.segments import make_segments
from cns.process.segments import cns_df_to_segments
from cns.utils.logging import log_info
from cns.utils.canonization import canonize_cns_df, canonize_sample_id
from cns.utils.conversions import segments_to_cns_df
from cns.utils.logging import log_warn, log_info
from cns.utils.assemblies import hg19


def _get_separator(path):
    if str(path).endswith(".csv"):
        return ","
    elif str(path).endswith(".tsv"):
        return "\t"
    else:
        raise ValueError(f"Unknown file format for file {path}, cannot determine separator.")


def load_cns(path, cn_columns=None, sep=None, sort=False, change_coords=True, order_columns=False, assembly=hg19, print_info=False):
    """
    Loads a CNS file into a pandas DataFrame.
    Loading includes canonization process, where the positions column names are standardized "sample_id", "chrom", "start", "end".
    CN columns are rename to one of ["major_cn", "minor_cn"], ["hap1_cn", "hap2_cn"], ["total_cn"].
    If these are not found, other typical column names are searched. If that fails, the first 4 columns are used as position, the following 1-2 as CNs.
    Coordinates are 1-based on input by default.

    Parameters
    ----------
    path : str
        Path to the CNS file.
    cn_columns : list of str, optional
        List of column names for copy number data. If None, columns are inferred from the file.
    sep : str, optional
        Separator for the file. If None, the separator is inferred from the file extension.
    sort : bool, optional
        If True, sorts the DataFrame by sample_id, chrom, and start.
    change_coords : bool, optional
        If True, changes the coordinates to 0-based.
    order_columns : bool, optional
        If True and there are two columns, the individual will be ordered as major/minor instead of their exising order and renamed to major_cn/minor_cn
    assembly : object, optional
        Genome assembly to use. Default is hg19.
    print_info : bool, optional
        If True, prints informational messages during processing.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the CNS data.
    """
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    sep = sep if sep is not None else _get_separator(path)
    cns_df = pd.read_csv(path, sep=sep, low_memory=False)
    cns_df = canonize_cns_df(cns_df, cn_columns, order_columns, assembly, print_info)
    if change_coords:
        cns_df.loc[:, "start"] -= 1
    if sort:
        cns_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    return cns_df


def save_cns(cns_df, path, sep=None, sort=False, change_coords=True, mode="w"):
    """
    Saves a CNS DataFrame to a file. Coordinates are 1-based on output by default.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing the CNS data.
    path : str
        Path to save the file.
    sort : bool, optional
        If True, sorts the DataFrame by sample_id, chrom, and start.
    change_coords : bool, optional
        If True, changes the coordinates to 1-based before saving.
    mode : str, optional
        Mode to open the file. Default is "w" (write). For append, header is not printed.

    Returns
    -------
    None
    """
    sep = sep if sep is not None else _get_separator(path)
    if sort:
        cns_df.sort_values(by=["sample_id", "chrom", "start"], inplace=True, ignore_index=True)
    if change_coords:
        cns_df.loc[:, "start"] += 1
    cns_df.to_csv(path, sep=sep, index=False, mode=mode, header=mode=="w")
    if change_coords:
        cns_df.loc[:, "start"] -= 1


def load_samples(path, sep=None, print_info=False):
    """
    Loads a samples file into a pandas DataFrame.
    Loading includes canonization process, where the index column is set to "sample_id". 
    The column is found by exact or similar name, if not found, the first column is used.
    If sex column is not found, it is added with value "NA".

    Parameters
    ----------
    path : str
        Path to the samples file.
    sep : str, optional
        Separator for the file. If None, the separator is inferred from the file extension.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the samples data.
    """
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    sep = sep if sep is not None else _get_separator(path)
    samples_df = pd.read_csv(path, sep=sep)
    samples_df = canonize_sample_id(samples_df, print_info=print_info)
    if "sex" not in samples_df.columns:
        log_info("'sex' column not found, adding column with value 'NA'.", suppress=not print_info)
        samples_df["sex"] = "NA"
    else:
        # where samples_df["sex"] is not xy or xx, replace with NA
        unknowns = samples_df[~samples_df["sex"].isin(["xy", "xx"])]["sex"].unique()
        if len(unknowns) > 0:
            log_info(f"Found unknown sex values: {unknowns}. Replacing with 'NA'. Use ['xx', 'xy'].", suppress=not print_info)
        samples_df.loc[~samples_df["sex"].isin(["xy", "xx"]), "sex"] = "NA"

    samples_df.set_index("sample_id", inplace=True)
    return samples_df   


def save_samples(samples_df, path, mode='w'):
    """
    Saves a samples DataFrame to a file.

    Parameters
    ----------
    samples_df : pandas.DataFrame
        DataFrame containing the samples data.
    path : str
        Path to save the file.
    mode : str, optional
        Mode to open the file. Default is "w" (write). For append, header is not printed.

    Returns
    -------
    None
    """
    samples_df.to_csv(path, sep="\t", index=True, mode=mode, header=mode=="w")


def fill_sex_if_missing(cns_df, samples_df):
    """
    Fills the sex column in the samples DataFrame if missing, based on the presence of chrY in the CNS data.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing the CNS data.
    samples_df : pandas.DataFrame
        DataFrame containing the samples data.

    Returns
    -------
    pandas.DataFrame
        Updated samples DataFrame with the sex column filled if missing.
    """
    res_df = samples_df.copy()
    # Set found_sex to True for each sample if there is chrY, otherwise set it to False
    found_sex = cns_df.groupby("sample_id")["chrom"].apply(lambda chroms: "chrY" in chroms.values)
    found_sex = found_sex.map({True: "xy", False: "xx"})
    if "sex" in res_df.columns:
        res_df["found_sex"] = found_sex
        condition = (res_df["sex"] == "xx") & (res_df["found_sex"] == "xy")
        indices = res_df[condition].index
        if len(indices) > 0:
            log_warn(f"Found samples where sex is xx in data but chrY has CNs assigned: {indices.tolist()}. "\
                    "This may result in an incorrect proportions of sex-chromosome features.")          
        res_df.drop(columns=["found_sex"], inplace=True)      
    # replace values in samples["sex"] with found_sex if samples["sex"] is not xy or xx
    mask = ~res_df["sex"].isin(["xy", "xx"])
    res_df.loc[mask, "sex"] = found_sex[res_df.index[mask]].values
    return res_df


def samples_df_from_cns_df(cns_df, fill_sex=True):
    """
    Creates a samples DataFrame (sample_is, sex) from a CNS DataFrame.

    Parameters
    ----------
    cns_df : pandas.DataFrame
        DataFrame containing the CNS data.
    fill_sex : bool, optional
        If True, fills the sex column in the samples DataFrame based on the presence of chrY in the CNS data.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the samples data.
    """
    ids = cns_df["sample_id"].unique()
    samples_df = pd.DataFrame({"sample_id": ids})
    samples_df["sex"] = "NA"
    samples_df.set_index("sample_id", inplace=True)
    if fill_sex:
        samples_df = fill_sex_if_missing(cns_df, samples_df)
    return samples_df


def save_segments(segs, path):    
    """
    Saves segments (chrom, start, end, name) to a 0-indexed BED file.

    Parameters
    ----------
    segs : list of tuples
        List of segments to save.
    path : str
        Path to save the file.

    Returns
    -------
    None
    """
    is_bed = path.lower().endswith(".bed")
    if not is_bed:
        log_warn(f"Segments file {path} is not bed file, the coordinates will be 1-based.")
    seg_df = segments_to_cns_df(segs)
    if not is_bed:
        seg_df = seg_df.copy()
        seg_df.loc[:, "start"] += 1
    sel = seg_df[["chrom", "start", "end", "name"]]
    sel.to_csv(path, sep="\t", index=False, header=not is_bed)


def load_segments(path):
    """
    Loads segments (chrom, start, end, name) from a file into a list of tuples.

    Parameters
    ----------
    path : str
        Path to the segments file.

    Returns
    -------
    list of tuples
        List of segments.
    """
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    is_bed = path.lower().endswith(".bed")
    if not is_bed:
        log_warn(f"Segments file {path} is not bed file, the coordinates will be 1-based.")
    if path == "" or path is None:
        return None
    path = abspath(path)
    if not exists(path):
        raise ValueError(f"File {path} not found.")
    
    # Read file, ignore lines that do not start with 'chr'
    with open(path, "r") as f:
        lines = [line for line in f if line.lstrip().startswith("chr")]
    segs_df = pd.read_csv(StringIO("".join(lines)), sep="\t", header=(None if is_bed else 0))
    # check that columns "chrom", "start" and "end" exist, more colums may be present
    if not is_bed:
        if not all([col in segs_df.columns for col in ["chrom", "start", "end"]]):
            raise ValueError(f"File {path} must have columns 'chrom', 'start' and 'end'.")
        if "name" not in segs_df.columns:
            segs_df["name"] = np.arange(len(segs_df))
    else:
        if len(segs_df.columns) < 3:
            raise ValueError(f"File {path} must have at least 3 columns.")
        elif len(segs_df.columns) == 3:            
            segs_df["name"] = np.arange(len(segs_df))
        elif len(segs_df.columns) > 4:
            log_warn(f"File {path} has more than 4 columns. Only the first 4 columns are used.")
            segs_df = segs_df.iloc[:, :4]                
        segs_df.columns = ["chrom", "start", "end", "name"]  
    if not is_bed:
        segs_df.loc[:, "start"] -= 1
    if len(segs_df.columns) == 3:
        segs_df["name"] = np.arange(len(segs_df))

    return cns_df_to_segments(segs_df)


def obtain_segments(segs_source, in_cols = None, assembly = hg19, print_info = False):
    if segs_source[-4:] == ".bed":
        log_info(f"Loading input file {segs_source}...", suppress=not print_info)
        return load_segments(segs_source)
    elif segs_source[-4:] == ".tsv":
        log_info(f"Loading CNS input file {segs_source}...", suppress=not print_info)
        input_cns = load_cns(segs_source, cn_columns=in_cols, assembly=assembly, print_info=print_info)
        return cns_df_to_segments(input_cns, process="unify")
    else:
        log_info(f"Creating {segs_source} segments...", suppress=not print_info)
        return make_segments(segs_source, assembly)