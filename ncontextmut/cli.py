"""Command-line interface for nContextMut."""

import argparse
import sys

from .core import run_analysis


def main():
    parser = argparse.ArgumentParser(
        prog="ncontextmut",
        description="Trinucleotide-context mutation analysis for multiple sequence alignments.",
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA file")
    parser.add_argument(
        "-o", "--output_folder", required=True, help="Path to output folder (absolute or relative)"
    )
    parser.add_argument(
        "-prefix",
        "--output_filename_prefix",
        required=True,
        help="Prefix for output files",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.6,
        help="Consensus frequency threshold (default: 0.6)",
    )

    args = parser.parse_args()

    try:
        run_analysis(
            fasta_path=args.input,
            output_folder=args.output_folder,
            output_prefix=args.output_filename_prefix,
            consensus_threshold=args.threshold,
        )
        print("Process completed successfully.")
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
