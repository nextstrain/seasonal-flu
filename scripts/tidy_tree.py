#!/usr/bin/env python3
"""
Tidy tree: build a phylogeny guided by the hierarchy of a lineage nomenclature.

A subtree is built for each lineage of the guide tree from its sequences and the
founders of the lineage and its child lineages, and the subtrees are grafted onto each
other at the founders. See the tidytree package for the details.
"""
import argparse
import sys

from Bio import Phylo, SeqIO

from tidytree.assign import (assign_founders, assign_sequences, load_assignments,
                             oldest_sequence, write_lineages)
from tidytree.build import collapse_short_branches, stitch_subtrees
from tidytree.guide import GuideTree


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--alignment', '-s', required=True,
                        help='Aligned input sequences in FASTA format')
    parser.add_argument('--founder-sequences', '-f',
                        help='Aligned lineage founder sequences in FASTA format, if they are not '
                             'in the alignment, where they are recognised by their lineage names')
    parser.add_argument('--guide-tree', '-g', required=True,
                        help='Lineage guide tree in Newick format')
    parser.add_argument('--lineage-assignments', '-a', required=True,
                        help='Table assigning sequences to lineages (TSV format)')
    parser.add_argument('--seq-id-column', default='seq_id',
                        help='Column name for sequence IDs in assignments file (default: seq_id)')
    parser.add_argument('--lineage-column', default='lineage',
                        help='Column name for lineage in assignments file (default: lineage)')
    parser.add_argument('--output-tree', '-o', required=True,
                        help='Output tree in Newick format')
    parser.add_argument('--model', default='GTR+G',
                        help='Substitution model for IQ-TREE (default: GTR+G)')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of subtrees built at the same time (default: 1)')
    parser.add_argument('--iqtree-args', default='--ninit 2 -n 2 -czb --epsilon 0.05',
                        help='Additional arguments to pass to IQ-TREE')
    parser.add_argument('--output-lineages',
                        help='Table of the lineage of each sequence in terms of the guide tree (TSV '
                             'with columns seq_id, lineage), from which the clades of the final '
                             'tree can be labelled')
    parser.add_argument('--keep-founders', action='store_true',
                        help='Keep founder sequences as leaves in the final tree')
    parser.add_argument('--root-lineage', default='auto',
                        help="Root lineage for the tree. 'auto' (default) uses the most specific "
                             "lineage of the guide tree that contains the lineages of all "
                             "sequences. With an explicit root lineage, sequences outside it are "
                             "excluded from the tree and reported.")
    parser.add_argument('--explicit-root',
                        help='Sequence on which to root the tree. If the root lineage has no founder '
                             'and this is not given, the oldest sequence of the root lineage is used, '
                             'which requires --metadata')
    parser.add_argument('--metadata',
                        help='Metadata with sampling dates, used to choose the oldest sequence of a '
                             'root lineage without founder as outgroup')
    parser.add_argument('--metadata-id-column', default='strain',
                        help='Column of --metadata holding the sequence IDs (default: strain)')
    parser.add_argument('--date-column', default='date',
                        help='Column of --metadata holding the sampling dates (default: date)')
    parser.add_argument('--default-lineage', default='unassigned',
                        help='Lineage of sequences without a lineage call (default: unassigned)')
    parser.add_argument('--ignore-missing-founders', action='store_true',
                        help='Fold lineages that lack a founder sequence into their parent '
                             '(default: raise error)')
    parser.add_argument('--unknown-lineage-policy', choices=['parent', 'root', 'drop'],
                        default='parent',
                        help='What to do with sequences whose lineage is not a node of the guide '
                             'tree: place them in the nearest enclosing lineage, in the root '
                             'lineage, or drop them (default: parent)')
    parser.add_argument('--max-dropped-fraction', type=float, default=0.05,
                        help='Fail if a larger fraction of the alignment ends up without a place in '
                             'the tree (default: 0.05)')
    parser.add_argument('--collapse-threshold', type=float, default=0.1,
                        help='Collapse internal branches shorter than this many substitutions into '
                             'polytomies, so that downstream steps need not remove branches. The '
                             'default, 0.1, is the length below which TreeTime, and thus augur '
                             'refine, collapses branches (default: 0.1)')
    return parser.parse_args()


def main():
    args = parse_arguments()

    sequences = {record.id: record for record in SeqIO.parse(args.alignment, 'fasta')}
    print(f"Loaded {len(sequences)} sequences")
    guide = GuideTree(args.guide_tree, args.default_lineage)

    if args.founder_sequences:
        founders = {record.id: record for record in SeqIO.parse(args.founder_sequences, 'fasta')}
        print(f"Loaded {len(founders)} founder sequences")
    else:
        # founders are named after their lineage; recognise them against the full guide
        # tree, not only the part below the root lineage
        founders = {seq_id: seq for seq_id, seq in sequences.items() if seq_id in guide}
        print(f"Gathered {len(founders)} founder sequences from the alignment")

    assignments = load_assignments(args.lineage_assignments, args.seq_id_column,
                                   args.lineage_column, args.default_lineage)
    print(f"Loaded {len(assignments)} assignments")

    root_lineage = args.root_lineage
    if root_lineage == 'auto':
        root_lineage = guide.automatic_root_lineage(
            guide.resolve(assignments[seq_id]) for seq_id in sequences
            if seq_id not in founders and seq_id in assignments)
        print(f"Root lineage chosen automatically as the lineage containing all sequences: "
              f"{root_lineage}")
    elif root_lineage not in guide:
        print(f"Error: root lineage '{root_lineage}' is not a lineage of the guide tree")
        return 1
    lineage_root = guide.lineage_tree(root_lineage)

    # founders of lineages outside the root lineage are not used
    below_root = set(guide.lineages_below(root_lineage))
    unused = sorted(seq_id for seq_id in founders if seq_id not in below_root)
    if unused:
        print(f"  ignoring {len(unused)} founders of lineages outside the root lineage: "
              f"{', '.join(unused)}")
        founders = {seq_id: seq for seq_id, seq in founders.items() if seq_id not in unused}
        sequences = {seq_id: seq for seq_id, seq in sequences.items() if seq_id not in unused}

    report = assign_sequences(lineage_root, sequences, assignments, guide,
                              args.unknown_lineage_policy, set(founders))
    total = len(sequences) - len(set(founders) & set(sequences))
    report.print(total)
    if not report.check(total, args.max_dropped_fraction):
        return 1

    folded = assign_founders(lineage_root, founders, args.ignore_missing_founders)
    if folded:
        print("\nLineages without a founder sequence, folded into their parent:")
        for name, parent, n_seqs in folded:
            print(f"  {name} -> {parent} ({n_seqs} sequences kept)")

    if args.output_lineages:
        n = write_lineages([seq_id for seq_id in sequences if seq_id not in founders],
                           assignments, guide, args.output_lineages)
        print(f"Wrote the lineages of {n} sequences to {args.output_lineages}")

    # the subtree of each lineage is rooted on its founder; a root lineage without one,
    # such as the context above the nomenclature, needs an outgroup instead
    explicit_root = args.explicit_root
    if lineage_root.founder_seq is None and explicit_root is None:
        if not args.metadata:
            print(f"Error: the root lineage '{lineage_root.name}' has no founder. Give an "
                  f"outgroup with --explicit-root or dates with --metadata to use the "
                  f"oldest sequence of the root lineage.")
            return 1
        explicit_root = oldest_sequence(lineage_root.sequences, args.metadata,
                                        args.metadata_id_column, args.date_column)
        if explicit_root is None:
            print(f"Error: the root lineage '{lineage_root.name}' has no founder and none "
                  f"of its sequences has a complete sampling date to root on.")
            return 1
        print(f"\nRoot lineage '{lineage_root.name}' has no founder, rooting on its oldest "
              f"sequence {explicit_root}")
    elif explicit_root is not None:
        print(f"\nRooting the root lineage '{lineage_root.name}' on {explicit_root}")

    print(f"\nBuilding subtrees using IQ-TREE, {args.threads} at a time...")
    tree = stitch_subtrees(lineage_root, {**sequences, **founders}, args.model,
                           args.iqtree_args, args.threads, args.keep_founders, explicit_root)
    if tree is None:
        print("Error: Failed to build final tree")
        return 1

    if not args.keep_founders:
        # grafting removes the founders of grafted lineages, but not that of the root
        # lineage or of lineages whose subtree is the founder alone
        for founder in [tip for tip in tree.get_terminals() if tip.name in founders]:
            tree.prune(founder)

    alignment_length = len(next(iter(sequences.values())).seq)
    threshold = args.collapse_threshold / alignment_length
    n_collapsed = collapse_short_branches(tree, threshold)
    print(f"\nCollapsed {n_collapsed} internal branches shorter than {args.collapse_threshold} "
          f"substitutions ({threshold:.2e} per site) or leading to a single child")
    tree.ladderize()

    Phylo.write(tree, args.output_tree, 'newick')
    print(f"Wrote the tree to {args.output_tree}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
