"""Prepare LibRPA input files from a calculation directory.

Currently supports ABACUS STRU or single-geometry running logs and writes
stru_out."""

import argparse
from pathlib import Path
import sys

from .detect import detect_source


def main(argv=None):
    parser = argparse.ArgumentParser(prog="taesc-prep", description=__doc__)
    parser.add_argument("directory", type=Path, help="calculation directory")
    parser.add_argument("-o", "--output", type=Path, metavar="DIRECTORY",
                        help="destination directory (default: input directory)")
    parser.add_argument("--no-symmetry", action="store_true",
                        help="omit spatial symmetry operations")
    args = parser.parse_args(argv)
    try:
        source = detect_source(args.directory)
        result = source.convert(args.output, include_symmetry=not args.no_symmetry)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    for note in result.notes:
        print(f"Note: {note}", file=sys.stderr)
    print(f"Wrote {result.output / 'stru_out'}: {result.n_atoms} atoms, "
          f"{result.n_symops} symmetry operations; no k-point block")
    return 0

