#!/usr/bin/env python3
"""
Prunes a reference strain from the provided tree.
"""
import argparse
from augur.io import read_sequences
from Bio import Phylo
import shutil
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", help="Newick tree to prune")
    parser.add_argument("--reference", nargs="?", help="FASTA file for the reference used to root the tree and to prune from the input")
    parser.add_argument("--output", help="Output Newick tree file")

    args = parser.parse_args()

    # If reference is not provided, then just copy the input to output without modifications
    if not args.reference:
        print("WARNING: No reference was provided, copying input tree to output tree", file=sys.stdout)
        shutil.copy(args.tree, args.output)
    else:
        # Open the reference sequence to get the name of the reference strain.
        reference = next(read_sequences(args.reference))
        reference_name = reference.id

        T = Phylo.read(args.tree, "newick")
        references = [c for c in T.find_clades(terminal=True) if c.name == reference_name]
        if references:
            T.prune(references[0])

        Phylo.write(T, args.output, "newick")
