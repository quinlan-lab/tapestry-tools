#!/usr/bin/env python3

import argparse
from rich_argparse import RichHelpFormatter
import sys 

from .methylation import compute_methylation_all_samples_at_given_loci
from .write_data import write_dataframe_to_bed

def main(): 
    class RichDefaultsFormatter(RichHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description="Compute founder-phased methylation at a set of loci for a set of samples",
        formatter_class=RichDefaultsFormatter
    )

    parser.add_argument(
        '--loci-bed', 
        type=str,
        required=True, 
        help='Bed file containing loci at which to compute methylation'
    )
    parser.add_argument(
        '--sample-meth-beds', 
        required=True, 
        type=str,
        help='File containing file paths of bed files containing founder-phased DNA methylation at CpG sites'
    )
    parser.add_argument(
        '--loci-meth-bed', 
        required=True, 
        type=str,
        help='File to store methylation at the given loci in the given samples'
    )
    parser.add_argument(
        '--testing', 
        action='store_true', 
        help='Testing mode, where program runs only 2 samples'
    )

    # Logic: If no arguments are passed, print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(1)

    args = parser.parse_args()

    df_meth_at_loci_all_samples = compute_methylation_all_samples_at_given_loci(
        args.loci_bed, args.sample_meth_beds, testing=args.testing
    )
    write_dataframe_to_bed(df_meth_at_loci_all_samples, args.loci_meth_bed, source=__file__)
    print(f"Wrote '{args.loci_meth_bed}'")
    print(f"Done running '{__file__}'")
    
if __name__ == '__main__': 
    main()