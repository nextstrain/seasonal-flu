#!/usr/bin/env python3
"""
Tidy Tree: Build phylogenetic trees guided by hierarchical lineages.

This program constructs a phylogenetic tree by:
1. Building subtrees for each lineage (with its sequences + child lineage founders)
2. Stitching subtrees together following the lineage guide tree
"""

import argparse
import io
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord


class LineageNode:
    """Represents a node in the lineage guide tree."""
    def __init__(self, name: str):
        self.tree = None
        self.name = name
        self.parent = None
        self.children = []
        self.sequences = []  # Sequence IDs belonging to this lineage
        self.founder_seq = None  # Founder sequence for this lineage
        self.subtree = None  # Phylogenetic subtree for this lineage


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Build phylogenetic tree guided by lineage hierarchy'
    )
    parser.add_argument(
        '--alignment', '-s',
        required=True,
        help='Aligned input sequences in FASTA format'
    )
    parser.add_argument(
        '--founder-sequences', '-f',
        required=False,
        help='Aligned lineage founder sequences in FASTA format (if not provided in alignment)'
    )
    parser.add_argument(
        '--guide-tree', '-g',
        required=True,
        help='Lineage guide tree in Newick format'
    )
    parser.add_argument(
        '--lineage-assignments', '-a',
        required=True,
        help='Table assigning sequences to lineages (TSV format)'
    )
    parser.add_argument(
        '--seq-id-column',
        default='seq_id',
        help='Column name for sequence IDs in assignments file (default: seq_id)'
    )
    parser.add_argument(
        '--lineage-column',
        default='lineage',
        help='Column name for lineage in assignments file (default: lineage)'
    )
    parser.add_argument(
        '--output-tree', '-o',
        required=True,
        help='Output tree in Newick format'
    )
    parser.add_argument(
        '--iqtree-path',
        default='iqtree',
        help='Path to IQ-TREE executable (default: iqtree)'
    )
    parser.add_argument(
        '--model',
        default='GTR+G',
        help='Substitution model for IQ-TREE (default: GTR+G)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads for IQ-TREE (default: 1)'
    )
    parser.add_argument(
        '--iqtree-args',
        default='--ninit 2 -n 2 --epsilon 0.05',
        help='Additional arguments to pass to IQ-TREE'
    )
    parser.add_argument(
        '--keep-founders',
        action='store_true',
        help='Keep founder sequences as leaves in the final tree'
    )
    parser.add_argument('--root-lineage', default=None, help='Root lineage for the tree')
    parser.add_argument('--explicit-root', default=None, type=str, help='Sequence from the context on which to root the tree')
    parser.add_argument('--default-lineage', default='unassigned', help='Default lineage for sequences with lineage assigned that evaluates to false (emty string, None, 0, etc.)')
    parser.add_argument(
        '--ignore-missing-founders',
        action='store_true',
        help='Ignore lineages in guide tree that lack founder sequences (default: raise error)'
    )
    return parser.parse_args()


def load_sequences(fasta_path: str) -> Dict[str, SeqRecord]:
    """Load sequences from FASTA file into a dictionary."""
    sequences = {}
    for record in SeqIO.parse(fasta_path, 'fasta'):
        sequences[record.id] = record
    return sequences


def load_assignments(
    assignment_path: str,
    seq_id_column: str = 'seq_id',
    lineage_column: str = 'lineage',
    default_lineage: str = 'unassigned'
) -> Dict[str, str]:
    """Load sequence-to-lineage assignments from TSV file using pandas."""
    # Read TSV file, handling comments
    df = pd.read_csv(assignment_path, sep='\t', comment='#')

    # Check if required columns exist
    if seq_id_column not in df.columns:
        raise ValueError(f"Column '{seq_id_column}' not found in assignments file. "
                        f"Available columns: {', '.join(df.columns)}")
    if lineage_column not in df.columns:
        raise ValueError(f"Column '{lineage_column}' not found in assignments file. "
                        f"Available columns: {', '.join(df.columns)}")

    # Fill missing or falsy lineage values with default_lineage
    df[lineage_column] = df[lineage_column].fillna(default_lineage)
    df.loc[~df[lineage_column].astype(bool), lineage_column] = default_lineage

    # Convert to dictionary
    assignments = dict(zip(df[seq_id_column], df[lineage_column]))

    return assignments


def parse_guide_tree(newick_path: str, root_lineage: str = None) -> LineageNode:
    """Parse the lineage guide tree and build LineageNode structure."""
    tree = Phylo.read(newick_path, 'newick')
    if root_lineage:
        matches = [x for x in tree.find_clades(name=root_lineage)]
        if matches:
            tree.root = matches[0]
            print(f"Re-rooted lineage tree to '{root_lineage}'")
        else:
            print(f"Warning: Specified root lineage '{root_lineage}' not found. Using original root.")

    # Convert Bio.Phylo tree to LineageNode structure
    def build_lineage_tree(phylo_node, parent=None, save_phylo_node=False):
        node = LineageNode(phylo_node.name or f"internal_{id(phylo_node)}")
        node.parent = parent
        if save_phylo_node:
            node.tree = phylo_node
        for child in phylo_node.clades:
            child_node = build_lineage_tree(child, node)
            node.children.append(child_node)

        return node

    return build_lineage_tree(tree.root, save_phylo_node=True)


def assign_sequences_to_lineages(
    root: LineageNode,
    assignments: Dict[str, str],
    sequences: Dict[str, SeqRecord]
):
    """Assign sequences to their respective lineage nodes."""
    # Create a map of lineage name to node
    lineage_map = {}

    def map_lineages(node):
        lineage_map[node.name] = node
        for child in node.children:
            map_lineages(child)

    map_lineages(root)

    # Assign sequences
    for seq_id, lineage in assignments.items():
        if lineage in lineage_map and seq_id in sequences:
            lineage_map[lineage].sequences.append(seq_id)


def assign_founders_to_lineages(
    root: LineageNode,
    founders: Dict[str, SeqRecord],
    ignore_missing_founders: bool = False
):
    """Assign founder sequences to their respective lineage nodes."""
    def assign(node):
        if node.name in founders:
            node.founder_seq = founders[node.name]
        children_to_keep = []
        # note that this allows a missing founder for the root node
        for child in node.children:
            if ignore_missing_founders and child.name not in founders:
                print(f"Warning: Founder sequence for lineage '{child.name}' not found. Skipping this lineage.")
            elif child.name not in founders:
                raise ValueError(f"Founder sequence for lineage '{child.name}' not found.")
            else:
                assign(child)
                children_to_keep.append(child)
        node.children = children_to_keep

    assign(root)


def run_iqtree(
    sequences: List[SeqRecord],
    iqtree_path: str,
    model: str,
    threads: int,
    root: str,
    extra_args: str = ''
) -> Optional[Phylo.BaseTree.Tree]:
    """
    Run IQ-TREE on a set of aligned sequences.

    Returns the resulting phylogenetic tree.
    """
    if len(sequences) < 3:
        return None

    # Create temporary directory for IQ-TREE files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        input_file = tmpdir_path / 'input.fasta'

        # Write sequences to file
        SeqIO.write(sequences, input_file, 'fasta')

        # Build IQ-TREE command
        cmd = [
            iqtree_path,
            '-s', str(input_file),
            '-m', model,
            '-nt', str(threads),
            '-quiet',  # Reduce output
            '-redo',   # Overwrite existing results
        ]

        # Add extra arguments if provided
        if extra_args:
            cmd.extend(extra_args.split())

        # Run IQ-TREE
        try:
            result = subprocess.run(
                cmd,
                cwd=tmpdir,
                capture_output=True,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"IQ-TREE error: {e.stderr}")
            return None

        # Read the resulting tree
        tree_file = input_file.with_suffix('.fasta.treefile')
        if not tree_file.exists():
            print(f"IQ-TREE did not produce expected tree file: {tree_file}")
            return None

        tree = Phylo.read(str(tree_file), 'newick')
        tree.root_with_outgroup(root)
        return tree


def build_subtree(
    node: LineageNode,
    all_sequences: Dict[str, SeqRecord],
    iqtree_path: str,
    model: str,
    threads: int,
    extra_args: str,
    explicit_root: str = None
) -> Optional[Phylo.BaseTree.Tree]:
    """
    Build a phylogenetic subtree for a lineage using IQ-TREE.

    Includes:
    - Sequences assigned to this lineage
    - Founder sequence of this lineage
    - Founder sequences of child lineages
    - Contextual sequences (for root lineage only)
    """
    # Collect sequences for this subtree
    subtree_seqs = []

    # Add sequences belonging to this lineage
    for seq_id in node.sequences:
        if seq_id in all_sequences:
            subtree_seqs.append(all_sequences[seq_id])

    # Add founder of this lineage
    if node.founder_seq:
        subtree_seqs.append(node.founder_seq)

    # Add founders of child lineages
    for child in node.children:
        if child.founder_seq:
            subtree_seqs.append(child.founder_seq)

    # Need at least 3 sequences to build a tree
    print(f"Lineage {node.name}: total sequences for subtree: {len(subtree_seqs)}")
    print([seq.id for seq in subtree_seqs])
    if len(subtree_seqs) < 3:
        if node.name.startswith('internal_'):
            tree = Phylo.read(io.StringIO(f"({','.join(seq.id for seq in subtree_seqs)});"), 'newick')
        elif len(subtree_seqs) == 2:
            tree = Phylo.read(io.StringIO(f"({','.join(seq.id for seq in subtree_seqs)}){node.name};"), 'newick')
        elif len(subtree_seqs) == 1:
            tree = Phylo.read(io.StringIO(f"{subtree_seqs[0].id};"), 'newick')
    else:
        print(f"Building subtree for lineage {node.name} with {len(subtree_seqs)} sequences")

        # Build tree using IQ-TREE
        print(f"Using explicit root: {explicit_root}, node name: {node.name}")
        tree = run_iqtree(subtree_seqs, iqtree_path, model, threads, root=explicit_root or node.name, extra_args=extra_args)

    return tree


def find_clade_by_name(tree, name: str):
    """Find a clade/node in the tree by its name."""
    def search(clade):
        if clade.name == name:
            return clade
        for child in clade.clades:
            result = search(child)
            if result:
                return result
        return None

    return search(tree.root)


def graft_subtree(
    parent_tree: Phylo.BaseTree.Tree,
    child_tree: Phylo.BaseTree.Tree,
    connection_point_id: str,
    keep_founders: bool = False
):
    """
    Graft a child subtree onto a parent tree at the connection point.

    The connection point is the child lineage's founder sequence,
    which should appear as a leaf in the parent tree and as the root
    or an internal node in the child tree.
    """
    # Find the connection point in the parent tree (should be a leaf)
    parent_node = find_clade_by_name(parent_tree, connection_point_id)
    if not parent_node:
        print(f"Warning: Could not find {connection_point_id} in parent tree")
        return

    assert parent_node.is_terminal(), "Connection point in parent tree must be a leaf"

    # Find the connection point in the child tree
    child_node = find_clade_by_name(child_tree, connection_point_id)

    if not child_node:
        print(f"Warning: Could not find {connection_point_id} in child tree")
        return

    # Replace the leaf in parent tree with the subtree from child tree
    # If child_node is a leaf in child_tree, nothing to graft
    if child_tree.root == child_node and len(child_tree.root.clades) == 0:
        # The founder is a leaf in both trees, no grafting needed
        return
    elif child_tree.root == child_node:
        if keep_founders:
            clades_to_graft = [c for c in child_tree.root.clades]
            child_node.clades = []  # Clear clades to avoid duplication
            clades_to_graft.append(child_node)  # Retain the founder as a leaf
        else:
            clades_to_graft = child_tree.root.clades
    elif child_tree.root != child_node:
        clades_to_graft = [c for c in child_tree.root.clades if (c.name != connection_point_id) or keep_founders]

    # Replace parent_node's children with child_node's children
    parent_node.clades = list(clades_to_graft)
    if parent_node.name:
        parent_node.name = parent_node.name + "_internal"


def stitch_subtrees(
    root: LineageNode,
    all_sequences: Dict[str, SeqRecord],
    iqtree_path: str,
    model: str,
    threads: int,
    extra_args: str,
    keep_founders: bool = False,
    explicit_root: str = None
) -> Optional[Phylo.BaseTree.Tree]:
    """
    Stitch subtrees together following the guide tree structure.

    Process (post-order traversal):
    1. Build subtrees for children first (bottom-up)
    2. Build subtree for current node
    3. Graft child subtrees onto current subtree at founder positions
    """
    def process_node(node, explicit_root: str = None):
        # First process all children (post-order)
        for child in node.children:
            process_node(child, explicit_root=None)

        # Build subtree for this node
        node.subtree = build_subtree(
            node, all_sequences, iqtree_path, model, threads, extra_args=extra_args, explicit_root=explicit_root
        )

        # If this node has a subtree and children with subtrees,
        # graft child subtrees onto this subtree
        if node.subtree and node.children:
            for child in node.children:
                if child.subtree and child.founder_seq:
                    print(f"Grafting {child.name} subtree onto {node.name} at {child.founder_seq.id}")
                    graft_subtree(node.subtree, child.subtree, child.founder_seq.id, keep_founders)
                else:
                    print(f"Skipping grafting for {child.name} onto {node.name} (missing subtree or founder)")

    process_node(root, explicit_root=explicit_root)

    root.subtree.ladderize()
    return root.subtree


def main():
    """Main execution function."""
    args = parse_arguments()

    # Load input data
    print("Loading sequences...")
    sequences = load_sequences(args.alignment)
    print(f"  Loaded {len(sequences)} sequences")

    print("Parsing lineage guide tree...")
    lineage_root = parse_guide_tree(args.guide_tree, args.root_lineage)

    if args.founder_sequences:
        print("Loading founder sequences...")
        founders = load_sequences(args.founder_sequences)
        print(f"  Loaded {len(founders)} founder sequences")
    else:
        founders = {}
        lineages_in_guide_tree = [clade.name for clade in lineage_root.tree.find_clades() if clade.name]
        for seq_id, seq in sequences.items():
            if seq_id in lineages_in_guide_tree:
                founders[seq_id] = seq
        print(f"  gathered {len(founders)} founder sequences from alignment")

    print("Loading sequence-to-lineage assignments...")
    assignments = load_assignments(
        args.lineage_assignments,
        args.seq_id_column,
        args.lineage_column,
        args.default_lineage
    )
    print(f"  Loaded {len(assignments)} assignments")

    # Prepare lineage structure
    print("\nAssigning sequences to lineages...")
    assign_sequences_to_lineages(lineage_root, assignments, sequences)
    assign_founders_to_lineages(lineage_root, founders, args.ignore_missing_founders)

    # Build and stitch trees
    print("\nBuilding subtrees using IQ-TREE...")
    all_sequences = {**sequences, **founders}
    final_tree = stitch_subtrees(
        lineage_root,
        all_sequences,
        args.iqtree_path,
        args.model,
        args.threads,
        args.iqtree_args,
        args.keep_founders,
        args.explicit_root,
        )

    if final_tree is None:
        print("Error: Failed to build final tree")
        return 1

    # Write output
    print(f"\nWriting final tree to {args.output_tree}...")
    Phylo.write(final_tree, args.output_tree, 'newick')

    print("Done!")
    return 0


if __name__ == '__main__':
    exit(main())
