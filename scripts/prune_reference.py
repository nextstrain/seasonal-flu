#!/usr/bin/env python3
"""
Prunes a reference strain from the provided tree.
"""
import argparse
from Bio import Phylo


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", help="Newick tree to prune")
    parser.add_argument("--reference", help="Name of the reference strain to prune")
    parser.add_argument("--output", help="Output Newick tree file")

    args = parser.parse_args()

    T = Phylo.read(args.tree, "newick")
    references = [ c for c in T.find_clades() if str(c.name) == args.reference ]
    if references:
        T.prune(references[0])

    Phylo.write(T, args.output, "newick")

