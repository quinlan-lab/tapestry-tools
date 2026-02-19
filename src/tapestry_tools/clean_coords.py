#!/usr/bin/env python3

import argparse
from rich_argparse import RichHelpFormatter
import sys


def convert_coords(coord_string):
    chrom, positions = coord_string.split(":")
    start_end = positions.replace(",", "").split("-")
    result = f"{chrom}\t{start_end[0]}\t{start_end[1]}"
    return result


def main():
    class RichDefaultsFormatter(
        RichHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
    ):
        pass

    parser = argparse.ArgumentParser(
        description="Convert genomic coordinates from colon-separated to tab-separated format.",
        formatter_class=RichDefaultsFormatter,
    )

    parser.add_argument(
        "coordinate",
        type=str,
        help="Coordinate string in format 'chr:start-end' (e.g., 'chr1:153,617,723-153,617,921')",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(1)

    args = parser.parse_args()

    result = convert_coords(args.coordinate)
    print(result)
