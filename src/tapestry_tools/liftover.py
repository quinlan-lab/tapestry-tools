#!/usr/bin/env python3

import argparse
from rich_argparse import RichHelpFormatter
import sys
from pyliftover import LiftOver


def lift_coord(lo, chrom, start, end):
    converted_start = lo.convert_coordinate(chrom, start - 1)
    converted_end = lo.convert_coordinate(chrom, end - 1) if end != start + 1 else None

    if converted_start and len(converted_start) > 0:
        new_chrom, new_start = converted_start[0][:2]
        new_start_bed = new_start + 1
        if converted_end and len(converted_end) > 0:
            new_end_bed = converted_end[0][1] + 1
            print(f"{new_chrom}\t{new_start_bed}\t{new_end_bed}")
        else:
            print(f"{new_chrom}\t{new_start_bed}\t{new_start_bed + 1}")


def parse_coordinate(coord_str):
    coord_str = coord_str.replace(",", "")  # clean IGV-style coordinate formatting
    parts = coord_str.split(":")
    if len(parts) != 2:
        sys.stderr.write(
            f"Error: Invalid coordinate format '{coord_str}'. Expected 'chr:start' or 'chr:start-end'.\n"
        )
        sys.exit(1)

    chrom = parts[0]
    range_parts = parts[1].split("-")

    if len(range_parts) == 1:
        start = int(range_parts[0])
        end = start + 1
    elif len(range_parts) == 2:
        start = int(range_parts[0])
        end = int(range_parts[1])
    else:
        sys.stderr.write(
            f"Error: Invalid coordinate format '{coord_str}'. Expected 'chr:start' or 'chr:start-end'.\n"
        )
        sys.exit(1)

    return chrom, start, end


def main():
    class RichDefaultsFormatter(
        RichHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
    ):
        pass

    parser = argparse.ArgumentParser(
        description="Lift over genomic coordinates from one reference to another.",
        formatter_class=RichDefaultsFormatter,
    )

    parser.add_argument(
        "--old-ref",
        required=True,
        type=str,
        help="Old reference genome (e.g., 'hg19', 'mm9').",
    )
    parser.add_argument(
        "--new-ref",
        required=True,
        type=str,
        help="New reference genome (e.g., 'hg38', 'mm10').",
    )
    parser.add_argument(
        "--coord",
        type=str,
        help="Single coordinate in format 'chr:start' or 'chr:start-end'.",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(1)

    args = parser.parse_args()

    lo = LiftOver(args.old_ref, args.new_ref)

    chrom, start, end = parse_coordinate(args.coord)
    
    lift_coord(lo, chrom, start, end)


if __name__ == "__main__":
    main()
