#!/usr/bin/env python3

import argparse
from rich_argparse import RichHelpFormatter
import sys 

from .imprinting import compute_delta_methylation_all_samples
from .write_data import write_dataframe_to_bed

def main(): 
    class RichDefaultsFormatter(RichHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description="Compute methylation difference between haplotypes for multiple samples",
        formatter_class=RichDefaultsFormatter
    )

    parser.add_argument(
        '--sample-meth-beds', 
        required=True, 
        type=str,
        help='File containing file paths of bed files containing founder-phased DNA methylation at CpG sites'
    )
    parser.add_argument(
        '--delta-meth-bed', 
        type=str,
        required=True, 
        help='Bed file containing difference of methylation levels between haplotypes for each tile and for each sample'
    )

    parser.add_argument(
        '--reference-genome', 
        type=str, 
        default="hg38", 
        help='Reference genome'
    )
    parser.add_argument(
        '--tile-size', 
        type=int, 
        default=1000, 
        help='Size in bp of genomic tile over which to average DNA methylation'
    )
    parser.add_argument(
        '--testing', 
        action='store_true', 
        help='Testing mode, where program takes only 2 mins to run for 2 samples'
    )

    # Logic: If no arguments are passed, print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(1)

    args = parser.parse_args()

    df_delta_meth_all_samples = compute_delta_methylation_all_samples(
        args.reference_genome, args.tile_size, args.sample_meth_beds, testing=args.testing
    )
    write_dataframe_to_bed(df_delta_meth_all_samples, args.delta_meth_bed, source=__file__)
    print(f"Wrote '{args.delta_meth_bed}'")
    print(f"Done running '{__file__}'")
    
if __name__ == '__main__': 
    main()