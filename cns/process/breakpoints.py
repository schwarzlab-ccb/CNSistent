import numpy as np
from cns.utils.assemblies import hg19
from cns.utils.conversions import cytobands_to_df


def split_into_bins(reg_len, step_size, strategy="after"):
    """
    Splits a region into bins of specified size using a given strategy.

    Parameters
    ----------
    reg_len : int
        Length of the region to split.
    step_size : int
        Size of each bin.
    strategy : str, optional
        Strategy to use for splitting. Options are "scale", "pad", "after". Default is "after".

    Returns
    -------
    list of int
        List of bin boundaries.
    """
    if (step_size < 1) or (reg_len < step_size):
        return [0, reg_len]
    padding = reg_len % step_size
    if strategy=="scale":
        step_count = reg_len // step_size
        if step_count < 1:
            return [0, reg_len]
        if padding >= (step_size / 2):
            step_size -= (step_size - padding) / (step_count + 1)
        else:
            step_size += padding / step_count
        fracs = np.arange(0, reg_len + 1, step_size) ## 
        return [np.int32(np.floor(frac + .5)) for frac in fracs]
    elif strategy=="pad":
        if padding >= (step_size / 2):
            padding = (step_size - padding) / 2
            start = step_size - padding
            end = reg_len - start
            fracs = np.arange(start, end + 1, step_size)
            return [0] + [np.int32(np.floor(frac + .5)) for frac in fracs] + [reg_len]
        else:
            if reg_len // step_size <= 1:
                return [0, reg_len]
            start = step_size + padding / 2
            end = reg_len - start
            fracs = np.arange(start, end + 1, step_size)
            return [0] + [np.int32(np.floor(frac + .5)) for frac in fracs] + [reg_len]
    else:       
        if padding >= (step_size / 2):
            fracs = np.arange(0, reg_len - padding + step_size, step_size)
        else:
            fracs = np.arange(0, reg_len - padding, step_size)
        return [np.int32(np.floor(frac + .5)) for frac in fracs] + [reg_len]
        

def _calc_genome_breaks(step_size, strategy="scale", assembly=hg19):
    return { chrom: split_into_bins(length, step_size, strategy) for chrom, length in assembly.chr_lens.items() }


def _calc_arm_breaks(assembly=hg19):
    cyto_df = cytobands_to_df(assembly.cytobands)
    acen = cyto_df.query("stain == 'acen' and name.str.contains('p')", engine="python")                                
    max_ends = cyto_df.groupby("chrom")["end"].max().to_dict()
    result = { row["chrom"]: [0, row["end"], max_ends[row["chrom"]]] for _, row in acen.iterrows() }
    return result


# all the breakpoints around cytobands
def _calc_cytoband_breaks(assembly=hg19):
    cyto_df = cytobands_to_df(assembly.cytobands)
    return { chrom: [0] + [end for end in cyto_df.query(f"chrom == '{chrom}'")["end"]] 
            for chrom in cyto_df["chrom"].unique() }

# Create breakpoints
def make_breaks(break_type, strategy='scale', assembly=hg19):
    """
    Creates breakpoints based on the specified break type.

    Parameters
    ----------
    break_type : str
        Type of break to use for creating breakpoints. Options are "arms", "bands", "whole", or a step size (e.g., "1MB").
    assembly : object, optional
        Genome assembly to use. Default is hg19.

    Returns
    -------
    dict
        Dictionary with chromosome names as keys and list of breaks as values.

    Raises
    ------
    ValueError
        If the break type is not recognized.
    """
    if break_type == "arms":
        return _calc_arm_breaks(assembly)
    elif break_type == "cytobands":
        return _calc_cytoband_breaks(assembly)
    else:
        try:
            step_size = int(break_type)
        except ValueError:
            raise ValueError("break_type must be 'arms', 'cytobands' or an integer, got " + break_type)
        return _calc_genome_breaks(step_size, strategy=strategy, assembly=assembly)
