#!/usr/bin/env python3

import argparse
from rich_argparse import RichHelpFormatter
from pathlib import Path
import polars as pl 
import sys 

from .read_data import read_dataframe_from_bed
from .imprinting import call_imprinted_loci
from .write_data import write_dataframe_to_bed

def main():
    class RichDefaultsFormatter(RichHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description="Call imprinted loci from a BED file of haplotype-specific methylation differences for a set of samples.",
        formatter_class=RichDefaultsFormatter
    )

    # Input/Output
    parser.add_argument(
        "--delta-meth-bed", 
        required=True, 
        type=str, 
        help="Path to the input BED file containing delta methylation data."
    )
    parser.add_argument(
        "--imprinted-bed", 
        required=True, 
        type=str, 
        help="Path to save the candidate imprinted loci."
    )

    # Logic Parameters
    parser.add_argument(
        "--meth-mode", 
        type=str, 
        default="model", 
        choices=["count", "model"],
        help="The means to compute methylation."
    )
    parser.add_argument(
        "--delta-meth-threshold", 
        type=float, 
        default=0.76, 
        help="Minimum absolute difference in methylation between haplotypes."
    )
    parser.add_argument(
        "--min-valid-cpgs-per-hap", 
        type=int, 
        default=5, 
        help="Minimum number of valid CpGs required per haplotype."
    )
    parser.add_argument(
        "--min-valid-cpg-ratio", 
        type=float, 
        default=0.5, 
        help="Minimum fraction of CpGs per haplotype that are valid."
    )

    # Logic: If no arguments are passed, print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(1)

    args = parser.parse_args()

    df = read_dataframe_from_bed(args.delta_meth_bed)

    df_candidates = call_imprinted_loci(
        df=df,
        meth_mode=args.meth_mode,
        delta_meth_threshold=args.delta_meth_threshold,
        num_valid_cpgs_per_hap_threshold=args.min_valid_cpgs_per_hap,
        valid_cpg_ratio_threshold=args.min_valid_cpg_ratio
    )

    df_candidates = df_candidates.with_columns(
        pl.col("imprinted_samples").list.join(",") 
    )

    imprinted_bed = Path(args.imprinted_bed)
    imprinted_bed.parent.mkdir(parents=True, exist_ok=True)
    imprinted_bed = str(imprinted_bed)

    write_dataframe_to_bed(df_candidates, imprinted_bed, source=__file__)

if __name__ == "__main__":
    main()