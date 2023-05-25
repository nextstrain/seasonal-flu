#!/usr/bin/env python3
import argparse
import sys
from treetime import TreeAnc
from augur.utils import read_tree, InvalidTreeError
import Bio.Phylo


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--trees", nargs="+", help="trees to sanitize by pruning leaves that do not appear in all trees.")
    parser.add_argument("--alignments", nargs="+", help="corresponding sequence alignments to remove short branches.")
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

    common_leaves = set.intersection(*[set(x.name for x in tree.find_clades(terminal=True)) for tree in trees])
    for output_tree_file, tree, aln in zip(args.output_trees, trees, args.alignments):
        for leaf in set(x.name for x in tree.find_clades(terminal=True)).difference(common_leaves):
            tree.prune(leaf)

        tt = TreeAnc(tree=tree, aln=aln)
        tt.infer_ancestral_sequences(infer_gtr=True)
        tt.prune_short_branches()
        tt.tree.ladderize()
        Bio.Phylo.write(tt.tree, output_tree_file, 'newick')
