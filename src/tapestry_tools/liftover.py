#!/usr/bin/env python3

import argparse
from rich_argparse import RichHelpFormatter
from pathlib import Path
import polars as pl
import sys

from .read_data import read_dataframe_from_bed
from .write_data import write_dataframe_to_bed


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
        "--input-bed",
        required=True,
        type=str,
        help="Input BED file with chrom, start, end columns (0-based).",
    )
    parser.add_argument(
        "--output-bed",
        required=True,
        type=str,
        help="Path to save the lifted-over coordinates as a BED file.",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stdout)
        sys.exit(1)

    args = parser.parse_args()

    from pyliftover import LiftOver

    lo = LiftOver(args.old_ref, args.new_ref)

    df = read_dataframe_from_bed(args.input_bed)

    converted_rows = []
    for row in df.iter_rows(named=True):
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]

        result = lo.convert_coordinate(chrom, start)

        if result is not None and len(result) > 0:
            for r in result:
                new_chrom, new_pos, strand, original_pos = r
                converted_rows.append(
                    {
                        "chrom": new_chrom,
                        "start": new_pos,
                        "end": new_pos + 1,
                        "original_chrom": chrom,
                        "original_start": start,
                        "original_end": end,
                        "strand": strand,
                        "original_pos": original_pos,
                    }
                )
        else:
            converted_rows.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": start + 1,
                    "original_chrom": chrom,
                    "original_start": start,
                    "original_end": end,
                    "strand": ".",
                    "original_pos": start,
                }
            )

    df_converted = pl.DataFrame(converted_rows)

    output_path = Path(args.output_bed)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    write_dataframe_to_bed(df_converted, str(output_path), source=__file__)
