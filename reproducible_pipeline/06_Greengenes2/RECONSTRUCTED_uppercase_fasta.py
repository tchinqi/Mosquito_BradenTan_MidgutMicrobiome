#!/usr/bin/env python3

"""
RECONSTRUCTED CODE

The original command used to uppercase the final LC/MC VSEARCH97 centroid
FASTA files was not preserved as standalone source code.

This reconstruction preserves FASTA header lines exactly and converts only
sequence lines to uppercase.

Validation performed 2026-10-06:
- LC uppercase FASTA reproduced the historical output byte-for-byte.
- MC uppercase FASTA reproduced the historical output byte-for-byte.
"""

import argparse
from pathlib import Path


def uppercase_fasta(input_path: Path, output_path: Path) -> None:
    with input_path.open("r") as inp, output_path.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                out.write(line)
            else:
                out.write(line.upper())


def main():
    parser = argparse.ArgumentParser(
        description="Uppercase FASTA sequence lines while preserving headers."
    )
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_fasta", type=Path)
    args = parser.parse_args()

    uppercase_fasta(args.input_fasta, args.output_fasta)


if __name__ == "__main__":
    main()
