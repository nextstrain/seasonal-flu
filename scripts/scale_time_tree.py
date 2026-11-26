#!/usr/bin/env python3
"""
Write a time tree with its branch lengths in substitutions per site.

TreeTime writes the time tree with branch lengths in years. augur refine wrote it with
the time of each branch multiplied by the clock rate, which is what the steps after it,
such as the ancestral reconstruction and the trait inference, expect.
"""
import argparse
import json

from Bio import Phylo

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tree', required=True, help='Time tree in Newick format, in years')
    parser.add_argument('--node-data', required=True,
                        help='Node data JSON of the time tree with the clock rate')
    parser.add_argument('--output', required=True,
                        help='Time tree with branch lengths in substitutions per site')
    args = parser.parse_args()

    rate = json.load(open(args.node_data))['clock']['rate']
    tree = Phylo.read(args.tree, 'newick')
    for node in tree.find_clades():
        node.branch_length = (node.branch_length or 0.0) * rate
    Phylo.write(tree, args.output, 'newick', format_branch_length='%1.8f')
