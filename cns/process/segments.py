import numpy as np
from cns.utils.assemblies import hg19
from cns.process.breakpoints import split_into_bins, make_breaks
from cns.utils.conversions import breaks_to_segments, cytobands_to_df, genome_to_segments, segments_to_breaks, tuples_to_segments


def count_segments(segs):
    """
    Counts the total number of segments in the input dictionary.

    Parameters
    ----------
    segs : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.

    Returns
    -------
    int
        Total number of segments.
    """
    return sum([len(chr_segs) for chr_segs in segs.values()])


def do_segments_overlap(segs, is_sorted=False):
    """
    Checks if any segments overlap in the input dictionary.

    Parameters
    ----------
    segs : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    is_sorted : bool, optional
        If True, assumes the segments are already sorted. Default is False.

    Returns
    -------
    bool
        True if any segments overlap, False otherwise.
    """
    for chr, chr_segs in segs.items():
        if not is_sorted:
            chr_segs.sort(key=lambda x: (x[0]))
        # Check for overlaps
        for i in range(len(chr_segs) - 1):
            current_end = chr_segs[i][1]
            next_start = chr_segs[i + 1][0]
            if current_end > next_start:
                return True
    return False


def find_overlaps(segs, is_sorted=False):
    """
    Finds overlapping segments in the input dictionary.

    Parameters
    ----------
    segs : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    is_sorted : bool, optional
        If True, assumes the segments are already sorted. Default is False.

    Returns
    -------
    dict
        Dictionary of overlapping segments with chromosome names as keys and list of overlapping segments as values.
    """
    overlaps = {}
    for chr, chr_segs in segs.items():
        if not is_sorted:
            chr_segs.sort(key=lambda x: (x[0]))
        # Iterate through all pairs of triplets to check for overlap
        n = len(chr_segs)
        for i in range(n):
            current_end = chr_segs[i][1]
            for j in range(i + 1, n):
                next_start = chr_segs[j][0]
                if current_end <= next_start:
                    break
                # Store the overlap along with the group identifiers
                if chr not in overlaps:
                    overlaps[chr] = []
                overlaps[chr].append((next_start, current_end))

    return overlaps


def merge_segments(segs, is_sorted=False):
    """
    Merges overlapping segments in the input dictionary.

    Parameters
    ----------
    segs : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    is_sorted : bool, optional
        If True, assumes the segments are already sorted. Default is False.

    Returns
    -------
    dict
        Dictionary of merged segments with chromosome names as keys and list of merged segments as values.
    """
    merged = {}
    for chr, chr_segs in segs.items():
        if len(chr_segs) == 0:
            merged[chr] = []
            continue

        if not is_sorted:
            chr_segs.sort(key=lambda x: (x[0]))

        merged[chr] = [chr_segs[0]]

        for current in chr_segs[1:]:
            last_start, last_end = merged[chr][-1][0], merged[chr][-1][1]
            last_name = merged[chr][-1][2] if len(merged[chr][-1]) > 2 else None

            # If the current segment starts at the end of the last one
            if current[0] <= last_end:
                name = current[2] if len(current) > 2 else None
                if last_name is None and name is not None:
                    last_name = name
                elif last_name is not None and name is not None:
                    last_name = f"{last_name}_{name}"                
                merged[chr][-1] = (last_start, current[1]) if last_name is None else (last_start, current[1], last_name)
            else:
                # Add the current segment as is
                merged[chr].append(current)

    return merged


def segment_union(segs_a, segs_b, merge=True):
    """
    Computes the union of two sets of segments.

    Parameters
    ----------
    segs1 : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    segs2 : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.

    Returns
    -------
    dict
        Dictionary of segments representing the union of the input segments.
    """
    # Combine and sort the segments first by group, then by start time
    keys = set(segs_a.keys()).union(set(segs_b.keys()))
    new_segs = {}
    for key in keys:
        new_segs[key] = segs_a.get(key, []) + segs_b.get(key, [])
        new_segs[key].sort(key=lambda x: (x[0]))
    if merge:
        new_segs = merge_segments(new_segs)
    return new_segs


def segment_difference(segs_a, segs_b, sorted=False):
    """
    Computes the difference between two sets of segments.

    Parameters
    ----------
    segs1 : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    segs2 : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.

    Returns
    -------
    dict
        Dictionary of segments representing the difference between the input segments.
    """
    diffs = {}
    for chr, chr_segs_a in segs_a.items():
        if chr not in segs_b:
            diffs[chr] = chr_segs_a
            continue
        chr_segs_b = segs_b[chr]

        if not sorted:
            chr_segs_a.sort(key=lambda x: (x[0]))
            chr_segs_b.sort(key=lambda x: (x[0]))

        diffs[chr] = []

        # Iterate through each segment in chr_segs_a
        for seg_a in chr_segs_a:
            new_start = seg_a[0]
            name_a = seg_a[2] if len(seg_a) > 2 else None
            subsections = []
            for seg_b in chr_segs_b:
                # Skip chr_segs_b that are in a different group or before the current segment in chr_segs_a
                if seg_b[1] < new_start:
                    continue
                # Break if the segment in chr_segs_b is beyond the current segment in chr_segs_a
                if seg_b[0] > seg_a[1]:
                    break

                # Calculate the difference if there's an overlap
                if seg_b[0] <= new_start < seg_b[1]:
                    # If chr_segs_a starts within chr_segs_b, move its start to the end of chr_segs_b
                    new_start = seg_b[1]
                elif new_start < seg_b[0] and seg_a[1] > seg_b[0]:
                    # If chr_segs_a overlaps the start of chr_segs_b, add the non-overlapping part to the difference
                    new_seg = (new_start, seg_b[0], name_a) if name_a is not None else (new_start, seg_b[0])
                    subsections.append(new_seg)
                    new_start = seg_b[1]

            # Check if there's any remaining part of chr_segs_a after processing overlaps
            if new_start < seg_a[1]:
                new_seg = (new_start, seg_a[1], name_a) if name_a is not None else (new_start, seg_a[1])
                subsections.append(new_seg)

            if len(subsections) > 1 and len(subsections[0]) > 2:
                subsections = [
                    (subsections[i][0], subsections[i][1], f"{subsections[i][2]}_{i}") for i in range(len(subsections))
                ]

            diffs[chr] += subsections

    return diffs


def pad_segments(chr_segs, pad_size, assembly=hg19):
    """
    Pads segments by a specified size on both sides.

    Parameters
    ----------
    segs : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    pad_size : int
        Size to pad segments on both sides.
    assembly : object, optional
        Genome assembly to use. Default is `hg19`.

    Returns
    -------
    dict
        Dictionary of padded segments.
    """
    res = {}
    for chrom, seg_group in chr_segs.items():
        res[chrom] = []
        for seg in seg_group:
            new_start = max(0, seg[0] - pad_size)
            new_end = min(seg[1] + pad_size, assembly.chr_lens[chrom])
            if len(seg) > 2:
                res[chrom].append((new_start, new_end, seg[2]))
            else:
                res[chrom].append((new_start, new_end))
    return res



def filter_cons_size(chr_segs, min_size):
    """
    Filters segments based on a minimum size, keeping only consecutive segments.

    Parameters
    ----------
    segs : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    min_size : int
        Minimum size of segments to keep.

    Returns
    -------
    dict
        Dictionary of filtered segments.
    """
    res = {}
    for chrom, seg_group in chr_segs.items():
        cons_groups = get_consecutive_segs(seg_group)
        res[chrom] = []
        for con_group in cons_groups:
            if con_group[-1][1] - con_group[0][0] >= min_size:
                res[chrom] += con_group
    return res        


def filter_min_size(chr_segs, min_size, merge_first=False):
    """
    Filters segments based on a minimum size.

    Parameters
    ----------
    segs : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    min_size : int
        Minimum size of segments to keep.

    Returns
    -------
    dict
        Dictionary of filtered segments.
    """
    if merge_first:
        chr_segs = merge_segments(chr_segs)
    return {chr: [seg for seg in chr_segs if seg[1] - seg[0] >= min_size] for chr, chr_segs in chr_segs.items()}


def split_segment(seg_start, seg_end, seg_name, step_size, strategy="scale"):
    """
    Splits a segment into smaller segments of specified size.

    Parameters
    ----------
    start : int
        Start position of the segment.
    end : int
        End position of the segment.
    step_size : int
        Size of each smaller segment.
    strategy : str, optional
        Strategy to use for splitting. Options are "scale", "pad", "after". Default is "after".

    Returns
    -------
    list of tuples
        List of smaller segments.
    """
    length = seg_end - seg_start
    breaks = split_into_bins(length, step_size, strategy)
    breaks = (np.array(breaks) + seg_start).tolist()
    if seg_name == None:
        res = [(breaks[i], breaks[i + 1]) for i in range(len(breaks) - 1)]
    else:
        res = [(breaks[i], breaks[i + 1], f"{seg_name}_{i}") for i in range(len(breaks) - 1)]
    return res


def split_segments(segments, step_size, strategy="scale"):
    """
    Splits segments into smaller segments of specified size.

    Parameters
    ----------
    segments : dict
        Dictionary of segments with chromosome names as keys and list of segments as values.
    step_size : int
        Size of each smaller segment.
    strategy : str, optional
        Strategy to use for splitting. Options are "scale", "pad", "after". Default is "after".

    Returns
    -------
    dict
        Dictionary of smaller segments.
    """
    res = {}
    for chr, chr_segs in segments.items():
        res[chr] = []
        for seg in chr_segs:
            seg_name = seg[2] if len(seg) > 2 else None
            res[chr] += split_segment(seg[0], seg[1], seg_name, step_size, strategy)
    return res


def _get_centromeres(assembly):
    cent_range = { chrom: (10**9, -1 ) for chrom in assembly.chr_lens.keys() }
    for band in assembly.cytobands:
        chrom = band[0]
        start = band[1]
        end = band[2]
        if "acen" in band[4]:
            cent_range[chrom] = (min(cent_range[chrom][0], start), max(cent_range[chrom][1], end))
    cents = { chrom: [(cent_range[chrom][0], cent_range[chrom][1], f"{chrom}_cen")] for chrom in cent_range.keys() }
    return cents


def make_segments(segment_type, assembly=hg19):
    """
    Selects and returns specific genomic regions of the selected type.

    Parameters
    ----------
    select : str
        The selection criteria for the regions. Options include:
        - "": Returns an empty dictionary.
        - "arms": Returns chromosome arms.
        - "bands": Returns cytogenetic bands.
        - "whole": Returns whole genome segments.
        - "gaps": Returns gap regions.
        - "centromeres": Returns centromeric regions.
        - "chrX": Returns the whole chromosome X (replace X with the chromosome number or name).
    assembly : object, optional
        The genome assembly to use. Default is hg19.

    Returns
    -------
    dict
        A dictionary with chromosome names as keys and lists of tuples representing the regions as values.
        Each tuple contains (start, end, name).

    Raises
    ------
    ValueError
        If an invalid chromosome is specified in the select parameter.
    """
    if segment_type == "":
        return {}
    if segment_type == "arms":
        arm_breaks = make_breaks("arms", assembly=assembly)
        return {
            chrom: [(breaks[0], breaks[1], f"{chrom}p"), (breaks[1], breaks[2], f"{chrom}q")]
            for chrom, breaks in arm_breaks.items()
        }
    elif segment_type == "bands":
        bands_df = cytobands_to_df(assembly.cytobands)
        return {
            chrom: [(start, end, name) for start, end, name in zip(subset["start"], subset["end"], subset["name"])]
            for chrom, subset in bands_df.groupby("chrom")
        }
    elif segment_type == "whole":
        return genome_to_segments(assembly)
    elif segment_type == "gaps":
        return tuples_to_segments(assembly.gaps)
    elif segment_type == "centromeres":
        return _get_centromeres(assembly)
    elif segment_type[0:3] == "chr":
        if segment_type not in assembly.aut_names:
            raise ValueError(f"Invalid chromosome {segment_type} for assembly {assembly.name}")
        return {segment_type: [(0, assembly.chr_lens[segment_type], segment_type)]}
    else:
        raise ValueError(f"Unable to create segments of type {segment_type}")


def get_consecutive_segs(segs):
    """
    Groups consecutive segments for each chromosome.

    Parameters
    ----------
    segs : list of tuples
        List of segments for a chromosome.

    Returns
    -------
    list of lists
        List of lists of consecutive segments.
    """
    if len(segs) == 0:
        return []
    res = []
    last_end = segs[0][1]
    res.append([segs[0]])
    for seg in segs[1:]:
        if seg[0] == last_end:
            res[-1].append(seg)
        else:
            res.append([seg])
        last_end = seg[1]
    return res


def align_segs_to_assembly(segs, sorted=False, assembly=hg19):
    """
    Aligns segments to the reference assembly by adding missing segments and adjusting boundaries.

    Parameters
    ----------
    segs : dict
        A dictionary with chromosome names as keys and lists of segments as values.
    sorted : bool, optional
        If True, sorts the segments for each chromosome. Default is False.
    assembly : object, optional
        The genome assembly to use. Default is hg19.

    Returns
    -------
    dict
        A dictionary with chromosome names as keys and lists of segments as values.
    """
    segs = segs.copy()
    for chr, end in assembly.chr_lens.items():
        if chr not in segs or len(segs[chr]) == 0:
            segs[chr] = [(0, end)]
        else:
            if not sorted:
                segs[chr].sort(key=lambda x: (x[0]))
            if segs[chr][0][0] > 0:
                segs[chr] = [(0, segs[chr][0][0], f"{chr}_start")] + segs[chr]
            if segs[chr][-1][1] < end:
                segs[chr].append((segs[chr][-1][1], end, f"{chr}_end"))
    return segs


def cns_df_to_segments(segs_df, process=None):
    """
    Converts a DataFrame of segments to a dictionary of segments by chromosome.

    Parameters
    ----------
    segs_df : pd.DataFrame
        A DataFrame containing segment information with columns "chrom", "start", "end", and optionally "name".
    process : str, optional
        The processing step to apply to the segments. Can be "merge" or "unify" or None.

    Returns
    -------
    dict
        A dictionary with chromosome names as keys and lists of segments as values.
    """
    segs = {}
    for chrom, group in segs_df.groupby("chrom"):
        names = group["name"] if "name" in group.columns else [f"{chrom}_{i}" for i in range(len(group))]
        segs[chrom] = list(zip(group["start"], group["end"], names))
    if process == "merge":
        return merge_segments(segs)
    elif process == "unify":
        input_breaks = segments_to_breaks(segs)
        return breaks_to_segments(input_breaks)
    elif process is not None:
        raise ValueError(f"Unknown segment processing method '{process}'")
    return segs


def add_arms(cns_df):    
    """
    Adds a column to the CNS DataFrame indicating the chromosome arm for each segment. (p or q)
    Parameters
    ----------
    cns_df : pd.DataFrame
        A DataFrame containing segment information with columns "chrom", "start", and "end".
    Returns
    -------
    pd.DataFrame
        The input DataFrame with an additional column "arm" indicating the chromosome arm.
    """
    arms = make_segments("arms")
    p_end = {chrom: segs[0][1] for chrom, segs in arms.items()}
    is_p = cns_df["end"] <= cns_df["chrom"].map(p_end)
    cns_df["arm"] = np.where(is_p, "p", "q")
    return cns_df