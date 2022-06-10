#!/usr/bin/env python3
import argparse
import sys

from augur.utils import read_tree, InvalidTreeError
import Bio.Phylo


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--trees", nargs="+", help="trees to sanitize by pruning leaves that do not appear in all trees.")
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
    for output_tree_file, tree in zip(args.output_trees, trees):
        for leaf in set(x.name for x in tree.find_clades(terminal=True)).difference(common_leaves):
            tree.prune(leaf)

        tree.root_at_midpoint()
        tree.ladderize()
        Bio.Phylo.write(tree, output_tree_file, 'newick')
