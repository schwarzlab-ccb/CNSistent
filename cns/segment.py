#!/usr/bin/env python

import argparse
import time

from cns.utils.files import obtain_segments, save_segments
from cns.utils.assemblies import get_assembly
from cns.utils.logging import log_info, set_verbose, handle_exception
from cns.utils.misc import parse_cncols, save_time
from cns.pipelines import main_segment


def _get_version():
    """Get version from pyproject.toml file."""
    try:
        with open("pyproject.toml", "r", encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("version"):
                    # Handles: version = "1.2.3"
                    return line.split("=")[1].strip().strip('"').strip("'")
    except Exception:
        return "unknown"


def _parse_args():
    """Parse command line arguments for segment action."""
    parser = argparse.ArgumentParser(
        description="Calculate segmentation regions for CNS data",
        prog="cns-segment"
    )
    
    # Add version argument
    parser.add_argument(
        "-v", "--version", 
        action="version", 
        version=f"%(prog)s {_get_version()}"
    )
    
    # Required arguments
    parser.add_argument(
        "data", 
        type=str, 
        help="Path to a TSV file with copy number segments, or a predefined segment type (e.g., 'whole', 'arms', 'bands')"
    )
    
    # Optional arguments
    parser.add_argument(
        "--out", 
        type=str, 
        help="Output file path", 
        required=False, 
        default="./segments.out.tsv"
    )
    
    parser.add_argument(
        "--assembly", 
        type=str, 
        help="Assembly to use. One of: hg19, hg38", 
        required=False, 
        default="hg19"
    )
    
    parser.add_argument(
        "--cncols",
        type=str,
        help="The name of either a single CN column or two comma separated columns. E.g. 'total_cn' or 'major_cn,minor_cn'",
        required=False,
        default=None,
    )
    
    parser.add_argument(
        "--split",
        type=int,
        help="Distance in which regions should be split, can be a positive integer or negative for no splitting (whole regions)",
        required=False,
        default=-1,
    )
    
    parser.add_argument(
        "--remove",
        type=str,
        help="Remove the regions after selection, before segmentation, can be either 'gaps', path to a BED file, or empty",
        required=False,
        default="",
    )
    
    parser.add_argument(
        "--filter",
        type=int,
        help="If set, regions smaller than the given size are excluded from selection and gaps. If negative, no filtering is done",
        required=False,
        default=-1,
    )
    
    parser.add_argument(
        "--merge",
        type=int,
        help="Maximum distance between breakpoint clusters for breakpoint merging. If negative, no breakpoints are merged",
        required=False,
        default=-1,
    )

    parser.add_argument(
        "--pad",
        type=int,
        help="Size in base pairs to pad removal segments on both sides. Default is 0 (no padding).",
        required=False,
        default=0,
    )
    
    parser.add_argument(
        "--keep-ends",
        help="If clustering (merge > 0) is enabled, do not cluster start and end breakpoints of each chromosome",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--align-to-assembly",
        help="Align segments to the assembly",
        action="store_true",
        default=False,
    )
    
    parser.add_argument(
        "--verbose", 
        help="Print progress to console", 
        action="store_true"
    )
    
    parser.add_argument(
        "--time", 
        help="Save runtime info", 
        action="store_true"
    )

    args = parser.parse_args()
    
    # Validation
    if args.filter is not None and args.filter == 0:
        raise ValueError("The --filter option must be negative (no filtering) or positive (minimum size).")
    
    if args.merge is not None and args.merge == 0:
        parser.error("The --merge option must be negative (no merging) or positive (maximum distance).")
    
    if args.split is not None and args.split == 0:
        parser.error("The --split option must be negative (no splitting) or positive (split distance).")

    return args



@handle_exception
def main():
    """Main function for segment action."""
    args = _parse_args()
    assembly = get_assembly(args.assembly)
    print_info = args.verbose
    set_verbose(print_info)  # Configure logging level based on verbose flag
    in_cols = parse_cncols(args.cncols)
    out_file = args.out

    log_info("***** segment *****")

    start = time.time()
    
    # Process segments
    log_info(f"Loading input segments from {args.data}...")
    input_segs = obtain_segments(args.data, in_cols, assembly, print_info)
    
    if args.remove:
        log_info(f"Loading regions to remove from {args.remove}...")
        remove_regs = obtain_segments(args.remove, in_cols, assembly, print_info)
    else:
        remove_regs = None
    
    log_info("Computing segmentation...")
    res_segs = main_segment(
        input_segs, 
        remove_regs, 
        args.split, 
        args.merge, 
        args.keep_ends,
        args.filter, 
        args.pad,
        args.align_to_assembly,
        assembly,
        print_info
    )
    
    runtime = time.time() - start
    
    if print_info:
        log_info(f"Finished in {runtime:.3f} seconds. Writing to {out_file}...")
        if args.time:
            save_time(out_file, runtime, start)
    
    save_segments(res_segs, out_file)
    log_info("Done.")


if __name__ == "__main__":
    main()


