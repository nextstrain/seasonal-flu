#!/usr/bin/env python3
import argparse
import sys

from augur.utils import read_tree, InvalidTreeError
from Bio import Phylo, AlignIO
from treetime import TreeTime, utils

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--trees", nargs="+", help="trees to sanitize by pruning leaves that do not appear in all trees.")
    parser.add_argument("--alignments", nargs="+", help="corresponding sequence alignments to remove short branches.")
    parser.add_argument("--metadata", help="metadata")
    parser.add_argument("--clock-filter", type=int, help="allowed clock deviation")
    parser.add_argument("--output-trees", nargs="+", help="sanitized trees, one for each input tree.")
    args = parser.parse_args()

    trees = []
    try:
        for tree_file in args.trees:
            tree = read_tree(tree_file)
            trees.append(tree)
    except InvalidTreeError as error:
        print(error, file=sys.stderr)
        sys.exit(1)

    alignments = [AlignIO.read(fname, 'fasta') for fname in args.alignments]
    dates = utils.parse_dates(args.metadata)

    for ti,tree in enumerate(trees):
        tt = TreeTime(tree=tree, dates=dates, aln=alignments[ti])
        tt.clock_filter(reroot='least-squares', n_iqd=args.clock_filter)

        for leaf in [l.name for l in tt.tree.get_terminals() if l.bad_branch==True]:
            tt.tree.prune(leaf)

    common_leaves = set.intersection(*[set(x.name for x in tree.find_clades(terminal=True)) for tree in trees])
    for ti,tree in enumerate(trees):
        for leaf in set([x.name for x in tree.get_terminals()]).difference(common_leaves):
            tree.prune(leaf)

        tt = TreeTime(tree=tree, dates=dates, aln=alignments[ti])
        tt.infer_ancestral_sequences(infer_gtr=True)
        tt.prune_short_branches()
        tree.ladderize()
        Phylo.write(tree, args.output_trees[ti], 'newick')
