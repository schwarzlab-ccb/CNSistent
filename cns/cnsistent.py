#!/usr/bin/env python

from multiprocessing import Pool
import time
import argparse
from os.path import exists

from cns.utils.misc import parse_cncols, save_time
from cns.utils.selection import dataframe_array_split
from cns.utils.logging import log_info, set_verbose, setup_mp_logging, stop_mp_logging, configure_worker_logging, handle_exception
from cns.pipelines import *


def _add_sp_args(action, parser):
    parser.add_argument("data", type=str, help="path to a TSV file with copy number segments, potentially with multiple samples.")
    parser.add_argument("--samples", type=str, help="path to the samples file", required=False, default="")
    parser.add_argument("--out", type=str, help="output file path", required=False, default="./cns.out.tsv")
    parser.add_argument("--assembly", type=str, help="assembly to use. One of: hg19, hg38.", required=False, default="hg19")
    parser.add_argument(
        "--cncols",
        type=str,
        help="The name of either a single CN column or two comma separated columns. E.g. 'total_cn' or 'major_cn,minor_cn'.",
        required=False,
        default=None,
    )
    parser.add_argument("--threads", type=int, help="number of threads to use", required=False, default=1)
    parser.add_argument(
        "--subsplit",
        type=int,
        help="will split the processing into chunks to lower memory needs",
        required=False,
        default=1
    )
    parser.add_argument("--verbose", help="print progress to console", action="store_true")
    parser.add_argument("--time", help="save runtime info", action="store_true")
    parser.add_argument(
        "--segments",
        type=str,
        help="A file with segments that create a mask over the CNS file. Preferably a .bed file.",
        required=False,
    )

    if action in ["infer", "impute"]:
        parser.add_argument(
            "--method",
            type=str,
            help='Inference method to use. Options are "extend", "diploid", or "zero". Default is "extend".',
            required=False,
            default="extend"
        )
        
    if action == "aggregate":
        parser.add_argument(
            "--how",
            type=str,
            help="The aggregation function, one of ['min', 'max', 'mean', 'none']",
            required=False,
            default="mean",
        )


def _get_version():
    try:
        with open("pyproject.toml", "r", encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("version"):
                    # Handles: version = "1.2.3"
                    return line.split("=")[1].strip().strip('"').strip("'")
    except Exception:
        return "unknown"


def _parse_args():
    # Top-level parser
    parser = argparse.ArgumentParser(description="Impute missing values in CNS data")
    subparsers = parser.add_subparsers(dest="action", help="cns action to perform. One of:")
    # Parse version from pyproject.toml

    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {_get_version()}")

    sp_dict = {}
    sp_dict["align"] = subparsers.add_parser("align", help=f"Adds Nan regions to the CNS data to match the assembly.")
    sp_dict["infer"] = subparsers.add_parser("infer", help=f"Infers values for NaNs in the CNS data.")
    sp_dict["impute"] = subparsers.add_parser("impute", help=f"Imputes missing values in the CNS data. (combines the align and infer commands)")
    sp_dict["coverage"] = subparsers.add_parser("coverage", help=f"Calculates coverage for aligned (but not imputed) CNS data." )
    sp_dict["ploidy"] = subparsers.add_parser("ploidy", help=f"Calculates the overall ploidy, nuli- and hemizygosity, and for 2 alleles also imballance (NaNs are ignored).")
    sp_dict["breakage"] = subparsers.add_parser("breakage", help=f"Conducts breakpoint analysis for CNS data (NaNs are ignored, but imputation is recommended).")
    sp_dict["aggregate"] = subparsers.add_parser("aggregate", help=f"Aggregate (bin) copy number values across samples to match provided segments.")
    for action, sp in sp_dict.items():
        _add_sp_args(action=action, parser=sp)

    args = parser.parse_args()
    if args.action is None:
        parser.print_help()
        exit(1)
        
    if args.action not in sp_dict:
        raise ValueError(f"Action {args.action} not recognized.")

    if not exists(args.data):
       raise ValueError(f"Data file {args.data} not found.")

    if args.threads <= 0:
        raise ValueError("The --threads option must be greater than 0.")

    if args.subsplit <= 0:
        raise ValueError("The --subsplit option must be greater than 0.")

    return args


def _action_to_fun(action):
    if action == "align":
        return main_align
    elif action == "infer":
        return main_infer
    elif action == "impute":
        return main_impute
    elif action == "coverage":
        return main_coverage
    elif action == "ploidy":
        return main_ploidy
    elif action == "breakage":
        return main_breakage
    elif action == "aggregate":
        return main_aggregate
    else:
        raise ValueError(f"Action {action} not recognized.")


def _get_blocks(action, input_block, samples_blocks, cn_cols, segs_block, assembly, args):
    block_count = len(input_block)
    # Apply process_block to each pair of blocks
    assembly_block = [assembly] * block_count
    ver_block = [False] * block_count
    ver_block[-1] = args.verbose
    cols_block = [cn_cols] * block_count
    if action == "infer":        
        method_block = [args.method] * block_count
        return zip(input_block, samples_blocks, cols_block, segs_block, method_block, ver_block)
    if action == "align":
        return zip(input_block, samples_blocks, cols_block, segs_block, assembly_block, ver_block)
    if action == "impute":
        method_block = [args.method] * block_count
        return zip(input_block, samples_blocks, cols_block, segs_block, method_block, assembly_block, ver_block)
    elif action in ["coverage", "ploidy", "breakage"]:
        return zip(input_block, samples_blocks, cols_block, segs_block, assembly_block, ver_block)
    elif action == "aggregate":
        if segs_block is None:
            raise ValueError("Segmentation blocks must be provided for this action.")
        fun_block = [args.how] * block_count
        return zip(input_block, segs_block, fun_block, cols_block, ver_block)
    else:
        raise ValueError(f"Unknown action {action}")


def _worker_init(queue, verbose):
    """Initialize logging in worker process."""
    configure_worker_logging(queue, verbose)


def _process(action, cns_df, samples_df, cn_cols, select_segs, assembly, args):
    main_fun = _action_to_fun(action)
    threads = abs(args.threads)
    samples_blocks = dataframe_array_split(samples_df, threads)
    cns_blocks = [cns_df.query("sample_id in @block.index").reset_index(drop=True) for block in samples_blocks]
    segs_blocks = [select_segs] * len(cns_blocks)
    zip_blocks = _get_blocks(action, cns_blocks, samples_blocks, cn_cols, segs_blocks, assembly, args)
    if threads == 1:
        return [main_fun(*list(*zip_blocks))]
    else:
        log_info(f"Multiprocessing with {threads} threads..")
        queue = setup_mp_logging()
        try:
            with Pool(threads, initializer=_worker_init, initargs=(queue, args.verbose)) as pool:
                res_blocs = pool.starmap(main_fun, zip_blocks)
                pool.close()
                pool.join()
        finally:
            stop_mp_logging()
        return res_blocs


@handle_exception
def main():
    args = _parse_args()
    action = args.action
    assembly = get_assembly(args.assembly)
    print_info = args.verbose
    set_verbose(print_info)  # Configure logging level based on verbose flag
    in_cols = parse_cncols(args.cncols)
    out_file = args.out

    log_info(f"***** cns {action} *****")
    log_info(f"Loading CNS input file {args.data}...")
    input_data = load_cns(args.data, cn_columns=in_cols, assembly=assembly, print_info=print_info)
    cn_columns = get_cn_cols(input_data, in_cols)
    if args.samples == "":
        samples_df = samples_df_from_cns_df(input_data, False)
    else:
        samples_df = load_samples(args.samples)
    samples_df = fill_sex_if_missing(input_data, samples_df)
    samples_blocks = dataframe_array_split(samples_df, args.subsplit)
    select_segs = load_segments(args.segments) if "segments" in args and args.segments is not None else None

    # Process blocks
    for i in range(args.subsplit):
        log_info(f"Processing block {i+1}/{args.subsplit}...")

        samples_block = samples_blocks[i]

        # Perform the action
        start = time.time()
        
        res_list = _process(action, input_data, samples_block, cn_columns, select_segs, assembly, args)
        runtime = time.time() - start

        if print_info:
            log_info(f"Finished in {runtime:.3f} seconds. Writing to {out_file} ...")
            if args.time:
                save_time(action, out_file, runtime, start, args.threads)

        for j in range(len(res_list)):
            mode = "w" if i == 0 and j == 0 else "a"
            res = res_list[j]
            if action in ["align", "infer", "impute", "aggregate"]:
                save_cns(res, out_file, change_coords=True, mode=mode)
            elif action in ["coverage", "ploidy", "breakage"]:
                save_samples(res, out_file, mode=mode)
            else:
                raise ValueError(f"Unknown action {action}")

    log_info("Done.")


if __name__ == "__main__":
    main()
