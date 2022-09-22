#!/usr/bin/env python3
"""
Prunes a reference strain from the provided tree.
"""
import argparse
from Bio import Phylo
import shutil
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", help="Newick tree to prune")
    parser.add_argument("--reference", nargs="?", help="Name of the reference strain to prune")
    parser.add_argument("--output", help="Output Newick tree file")

    args = parser.parse_args()

    # If reference is not provided, then just copy the input to output without modifications
    if not args.reference:
        print("WARNING: No reference was provided, copying input tree to output tree", file=sys.stdout)
        shutil.copy(args.tree, args.output)
    else:
        T = Phylo.read(args.tree, "newick")
        references = [ c for c in T.find_clades() if str(c.name) == args.reference ]
        if references:
            T.prune(references[0])

        Phylo.write(T, args.output, "newick")

