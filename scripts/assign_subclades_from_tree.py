#!/usr/bin/env python3
"""
Label the lineages of a tree from the lineages of its tips.

Every tip carries its lineage as a node of the guide tree. A lineage comprises the tips
of the lineage itself and of all lineages descending from it in the guide tree. Its
label goes on the branch leading to the MRCA of these tips, and every node in the
subtree below the MRCA is assigned to it. Lineages are processed from the root of the
guide tree downwards, so a sublineage overwrites its parent within its own subtree.

If the subtree of the MRCA also contains tips of other lineages, the lineage has no
node of its own, typically because the branch leading to it had zero length and was
collapsed into a polytomy. Then its largest subtrees that contain only its own tips
are labelled instead, so that it does not take over the rest of the polytomy.

Unlike labelling by signature mutations, this does not depend on the mutations
reconstructed on the tree, only on the tips that are present, so it is robust to
pruning tips and collapsing branches after the tree was built.

The output has the format of augur clades.
"""
import argparse
import json
from collections import Counter, defaultdict

import pandas as pd
from Bio import Phylo


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
    parser.add_argument('--output-node-data', required=True, help='Output node data JSON')
    return parser.parse_args()


def guide_tree_ancestry(path):
    """Parent of each lineage in the guide tree and the lineages in preorder."""
    guide = Phylo.read(path, 'newick')
    parent = {}
    preorder = []
    for clade in guide.find_clades(order='preorder'):
        if clade.name:
            preorder.append(clade.name)
        for child in clade.clades:
            if child.name:
                parent[child.name] = clade.name
    return parent, preorder


def lineage_and_ancestors(lineage, parent):
    while lineage is not None:
        yield lineage
        lineage = parent.get(lineage)


def mrca(nodes, up, depth):
    """Most recent common ancestor of nodes, given parent pointers and depths."""
    nodes = iter(nodes)
    ancestor = next(nodes)
    for node in nodes:
        while depth[node] > depth[ancestor]:
            node = up[node]
        while depth[ancestor] > depth[node]:
            ancestor = up[ancestor]
        while node is not ancestor:
            node, ancestor = up[node], up[ancestor]
    return ancestor


def main():
    args = parse_arguments()

    tree = Phylo.read(args.tree, 'newick')
    up, depth = {}, {tree.root: 0}
    for node in tree.find_clades(order='preorder'):
        for child in node.clades:
            up[child] = node
            depth[child] = depth[node] + 1

    placed = pd.read_csv(args.lineages, sep='\t', dtype=str)
    tip_lineage = dict(zip(placed['seq_id'], placed['lineage']))
    parent, preorder = guide_tree_ancestry(args.guide_tree)

    tips = tree.get_terminals()
    missing = [tip.name for tip in tips if tip.name not in tip_lineage]
    if missing:
        print(f"Warning: {len(missing)} tips have no lineage and count as {args.default_lineage} "
              f"(e.g. {', '.join(missing[:5])})")

    # the tips of each lineage, including those of its sublineages
    members = defaultdict(list)
    for tip in tips:
        for lineage in lineage_and_ancestors(tip_lineage.get(tip.name), parent):
            members[lineage].append(tip)

    # number of tips below each node, for telling whether a subtree is pure
    n_tips = {}
    for node in tree.find_clades(order='postorder'):
        n_tips[node] = 1 if node.is_terminal() else sum(n_tips[c] for c in node.clades)

    membership = {node: args.default_lineage for node in tree.find_clades()}
    labels = []   # (node, lineage), in the order lineages are assigned
    split = {}    # lineage -> (tips of other lineages under its MRCA, number of subtrees)
    for lineage in preorder:
        if lineage == args.default_lineage or lineage not in members:
            continue
        own = members[lineage]
        root = mrca(own, up, depth)
        if n_tips[root] == len(own):
            roots = [root]
        else:
            # count the tips of the lineage below each node of the MRCA subtree
            n_own = Counter()
            for tip in own:
                node = tip
                while True:
                    n_own[node] += 1
                    if node is root:
                        break
                    node = up[node]
            roots, stack = [], [root]
            while stack:
                node = stack.pop()
                if n_own[node] == n_tips[node]:
                    roots.append(node)
                elif n_own[node]:
                    stack.extend(node.clades)
            split[lineage] = (n_tips[root] - len(own), len(roots))
        for subtree_root in roots:
            labels.append((subtree_root, lineage))
            for node in subtree_root.find_clades():
                membership[node] = lineage

    print("Lineages without a node of their own, labelled on their largest pure subtrees:")
    for lineage, (n_other, n_subtrees) in split.items():
        print(f"  {lineage}: {len(members[lineage])} tips in {n_subtrees} subtrees, "
              f"{n_other} tips of other lineages under their MRCA")
    if not split:
        print("  none")

    # a branch keeps the label of the most specific lineage starting on it
    branch_labels = {node.name: lineage for node, lineage in labels}
    node_data = {
        'nodes': {node.name: {args.membership_name: lineage} for node, lineage in membership.items()},
        'branches': {name: {'labels': {args.label_name: lineage}} for name, lineage in branch_labels.items()},
    }
    with open(args.output_node_data, 'w') as fh:
        json.dump(node_data, fh, indent=1)

    if args.compare:
        other = json.load(open(args.compare))['nodes']
        tip_names = {tip.name for tip in tips}
        disagreements = Counter()
        for name in tip_names:
            ours = node_data['nodes'][name][args.membership_name]
            theirs = other.get(name, {}).get(args.membership_name)
            if ours != theirs:
                disagreements[(ours, theirs)] += 1
        n_agree = len(tip_names) - sum(disagreements.values())
        print(f"\nTip labels agree with {args.compare} for {n_agree} of {len(tip_names)} tips")
        for (ours, theirs), n in disagreements.most_common(20):
            print(f"  {n} tips: {ours} here, {theirs} there")


if __name__ == '__main__':
    main()
