#!/usr/bin/env python3
"""
Label the lineages of a tree from the lineages of its tips.

The label of a lineage goes on the MRCA of its tips and those of its sublineages; see
tidytree.labels for how lineages without a node of their own are handled. With
--refine-calls, tips take a call from the current nomenclature where it is their lineage
or a sublineage of it, which catches up calls made with an older nomenclature.

The output has the format of augur clades.
"""
import argparse
import json

import pandas as pd
from Bio import Phylo

from tidytree.guide import GuideTree
from tidytree.labels import compare, label_lineages, refine_calls


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tree', required=True, help='Tree in Newick format with named internal nodes')
    parser.add_argument('--lineages', required=True,
                        help='Lineage of each tip, TSV with columns seq_id and lineage, as '
                             'written by tidy_tree.py --output-lineages')
    parser.add_argument('--guide-tree', required=True, help='Lineage guide tree in Newick format')
    parser.add_argument('--membership-name', default='subclade',
                        help='Node attribute holding the lineage (default: subclade)')
    parser.add_argument('--label-name', default='Subclade',
                        help='Branch label marking the root of a lineage (default: Subclade)')
    parser.add_argument('--default-lineage', default='unassigned',
                        help='Lineage of nodes outside every labelled lineage (default: unassigned)')
    parser.add_argument('--compare', help='Node data JSON of another labelling, e.g. augur clades, '
                                          'to report disagreements with')
    parser.add_argument('--refine-calls',
                        help='Node data JSON with lineage calls of the tips from the current '
                             'nomenclature, e.g. augur clades with its subclades.tsv. A call is '
                             'taken where it is the lineage of the tip or one of its sublineages, '
                             'so that tips called with an older nomenclature pick up lineages '
                             'designated or renamed since')
    parser.add_argument('--refine-calls-membership', default='subclade',
                        help='Node attribute holding the lineage in --refine-calls (default: subclade)')
    parser.add_argument('--output-node-data', required=True, help='Output node data JSON')
    return parser.parse_args()


def main():
    args = parse_arguments()
    tree = Phylo.read(args.tree, 'newick')
    guide = GuideTree(args.guide_tree, args.default_lineage)
    placed = pd.read_csv(args.lineages, sep='\t', dtype=str)
    tip_lineage = dict(zip(placed['seq_id'], placed['lineage']))

    if args.refine_calls:
        calls = {name: data.get(args.refine_calls_membership)
                 for name, data in json.load(open(args.refine_calls))['nodes'].items()}
        refined = refine_calls(tip_lineage, calls, guide)
        print(f"Refined the lineage of {sum(refined.values())} tips with {args.refine_calls}:")
        for (lineage, call), n in refined.most_common():
            print(f"  {n} tips: {lineage} -> {call}")

    labelling = label_lineages(tree, tip_lineage, guide, args.default_lineage)
    if labelling.missing:
        print(f"Warning: {len(labelling.missing)} tips have no lineage and count as "
              f"{args.default_lineage} (e.g. {', '.join(labelling.missing[:5])})")
    labelling.print()

    node_data = labelling.node_data(args.membership_name, args.label_name)
    with open(args.output_node_data, 'w') as fh:
        json.dump(node_data, fh, indent=1)

    if args.compare:
        tip_names = [tip.name for tip in tree.get_terminals()]
        disagreements = compare(node_data, json.load(open(args.compare))['nodes'], tip_names,
                                args.membership_name)
        print(f"\nTip labels agree with {args.compare} for "
              f"{len(tip_names) - sum(disagreements.values())} of {len(tip_names)} tips")
        for (ours, theirs), n in disagreements.most_common(20):
            print(f"  {n} tips: {ours} here, {theirs} there")


if __name__ == '__main__':
    main()
