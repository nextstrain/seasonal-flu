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
import sys
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord

sys.setrecursionlimit(10000)

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
        default='--ninit 2 -n 2 -czb --epsilon 0.05',
        help='Additional arguments to pass to IQ-TREE'
    )
    parser.add_argument(
        '--output-lineages',
        help='Table of the lineage of each sequence in terms of the guide tree (TSV with '
             'columns seq_id, lineage), from which the clades of the final tree can be labelled'
    )
    parser.add_argument(
        '--keep-founders',
        action='store_true',
        help='Keep founder sequences as leaves in the final tree'
    )
    parser.add_argument('--root-lineage', default='auto',
                        help="Root lineage for the tree. 'auto' (default) uses the most specific "
                             "lineage of the guide tree that contains the lineages of all "
                             "sequences. With an explicit root lineage, sequences outside it are "
                             "excluded from the tree and reported.")
    parser.add_argument('--explicit-root', default=None, type=str,
                        help='Sequence on which to root the tree. If the root lineage has no founder '
                             'and this is not given, the oldest sequence of the root lineage is used, '
                             'which requires --metadata')
    parser.add_argument('--metadata', default=None,
                        help='Metadata with sampling dates, used to choose the oldest sequence of a '
                             'root lineage without founder as outgroup')
    parser.add_argument('--metadata-id-column', default='strain',
                        help='Column of --metadata holding the sequence IDs (default: strain)')
    parser.add_argument('--date-column', default='date',
                        help='Column of --metadata holding the sampling dates (default: date)')
    parser.add_argument('--default-lineage', default='unassigned', help='Default lineage for sequences with lineage assigned that evaluates to false (emty string, None, 0, etc.)')
    parser.add_argument(
        '--ignore-missing-founders',
        action='store_true',
        help='Fold lineages that lack a founder sequence into their parent (default: raise error)'
    )
    parser.add_argument(
        '--unknown-lineage-policy',
        choices=['parent', 'root', 'drop'],
        default='parent',
        help='What to do with sequences whose lineage is not a node of the guide tree: place them '
             'in the nearest enclosing lineage, in the root lineage, or drop them (default: parent)'
    )
    parser.add_argument(
        '--max-dropped-fraction',
        type=float,
        default=0.05,
        help='Fail if a larger fraction of the alignment ends up without a place in the tree '
             '(default: 0.05)'
    )
    parser.add_argument(
        '--collapse-threshold',
        type=float,
        default=1e-5,
        help='Collapse internal branches shorter than this into polytomies, so that downstream '
             'steps need not remove branches (default: 1e-5, well below one substitution)'
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


def parse_guide_tree(newick_path: str, root_lineage: str = None,
                     default_lineage: str = 'unassigned') -> LineageNode:
    """Parse the lineage guide tree and build LineageNode structure."""
    tree = Phylo.read(newick_path, 'newick')
    # as in read_guide_tree_parents, an unnamed root is the default lineage
    tree.root.name = tree.root.name or default_lineage
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


def oldest_sequence(
    seq_ids: List[str],
    metadata_path: str,
    id_column: str = 'strain',
    date_column: str = 'date',
) -> Optional[str]:
    """
    The sequence with the earliest fully specified sampling date.

    Used as outgroup for a root lineage without a founder, typically the context of
    sequences that predate the nomenclature. Ambiguous dates such as 2015-XX-XX are
    ignored, as they cannot be ordered reliably.
    """
    from augur.dates import get_numerical_date_from_value

    df = pd.read_csv(metadata_path, sep='\t', usecols=[id_column, date_column], dtype=str)
    dates = dict(zip(df[id_column], df[date_column]))
    dated = []
    for seq_id in seq_ids:
        value = dates.get(seq_id)
        if not isinstance(value, str):
            continue
        numdate = get_numerical_date_from_value(value, '%Y-%m-%d')
        if isinstance(numdate, float):
            dated.append((numdate, seq_id))
    if not dated:
        return None
    return min(dated)[1]


def read_guide_tree_parents(newick_path: str, default_lineage: str = 'unassigned'):
    """
    Parent of each lineage of the full guide tree, and the name of its root.

    An unnamed root takes the default lineage, which is what sequences without a
    lineage call are assigned to.
    """
    guide = Phylo.read(newick_path, 'newick')
    guide_root = guide.root.name or default_lineage
    parent = {guide_root: None}
    for clade in guide.find_clades(order='preorder'):
        name = guide_root if clade is guide.root else clade.name
        for child in clade.clades:
            if child.name:
                parent[child.name] = name
    return parent, guide_root


def resolve_lineage(lineage: Optional[str], parent: Dict[str, Optional[str]], guide_root: str) -> str:
    """
    The lineage of the full guide tree a lineage call corresponds to: the call itself,
    else its nearest enclosing lineage, else the root of the guide tree.
    """
    if lineage in parent:
        return lineage
    if isinstance(lineage, str):
        return enclosing_lineage(lineage, parent) or guide_root
    return guide_root


def automatic_root_lineage(lineages, parent: Dict[str, Optional[str]]) -> str:
    """The most specific lineage of the guide tree containing all given lineages."""
    def path_to_root(lineage):
        path = []
        while lineage is not None:
            path.append(lineage)
            lineage = parent[lineage]
        return path[::-1]

    common = None
    for lineage in set(lineages):
        path = path_to_root(lineage)
        if common is None:
            common = path
        else:
            n = 0
            while n < min(len(common), len(path)) and common[n] == path[n]:
                n += 1
            common = common[:n]
    return common[-1]


def enclosing_lineage(lineage: str, lineage_map: Dict[str, LineageNode]) -> Optional[str]:
    """
    Nearest ancestor of `lineage` that is present in the guide tree.

    Lineages that are added to the nomenclature after the guide tree was built,
    or that are more finely resolved than the guide tree, are placed in the
    lineage they descend from, e.g. 'K.1.2' is placed in 'K.1' or else in 'K'.
    """
    parts = lineage.split('.')
    for i in range(len(parts) - 1, 0, -1):
        candidate = '.'.join(parts[:i])
        if candidate in lineage_map:
            return candidate
    return None


def assign_sequences_to_lineages(
    root: LineageNode,
    assignments: Dict[str, str],
    sequences: Dict[str, SeqRecord],
    guide_parent: Dict[str, Optional[str]],
    guide_root: str,
    unknown_lineage_policy: str = 'parent',
    founder_ids: Optional[set] = None,
) -> Dict[str, object]:
    """
    Assign sequences to their respective lineage nodes.

    Every sequence of the alignment ends up in exactly one of the reported
    categories, so that sequences are never dropped without notice.
    """
    lineage_map = {}

    def map_lineages(node):
        lineage_map[node.name] = node
        for child in node.children:
            map_lineages(child)

    map_lineages(root)

    report = {
        'assigned': 0,
        'founders': 0,
        'placed_in_enclosing_lineage': {},   # lineage not in guide tree -> count
        'placed_in_root': {},                # no enclosing lineage in guide tree -> count
        'dropped_unknown_lineage': {},       # only when the policy says to drop
        'outside_root': {},                  # lineage outside the root lineage -> count
        'without_assignment': [],            # in the alignment but not in the metadata
    }

    founder_ids = founder_ids or set()
    for seq_id in sequences:
        if seq_id in founder_ids:
            # founders are placed by assign_founders_to_lineages, not by their lineage call
            report['founders'] += 1
            continue
        lineage = assignments.get(seq_id)
        if lineage is None:
            report['without_assignment'].append(seq_id)
            continue
        if lineage in lineage_map:
            lineage_map[lineage].sequences.append(seq_id)
            report['assigned'] += 1
            continue

        # a lineage of the guide tree outside the root lineage, or a call that is not
        # a node of the guide tree and resolves to one outside the root lineage: these
        # are excluded on purpose by the choice of the root lineage
        resolved = resolve_lineage(lineage, guide_parent, guide_root)
        if resolved not in lineage_map:
            report['outside_root'].setdefault(lineage, 0)
            report['outside_root'][lineage] += 1
            continue

        # the lineage is not a node of the guide tree but descends from the root lineage
        if unknown_lineage_policy == 'drop':
            report['dropped_unknown_lineage'].setdefault(lineage, 0)
            report['dropped_unknown_lineage'][lineage] += 1
        elif unknown_lineage_policy == 'parent' and resolved != root.name:
            lineage_map[resolved].sequences.append(seq_id)
            report['placed_in_enclosing_lineage'].setdefault(f"{lineage} -> {resolved}", 0)
            report['placed_in_enclosing_lineage'][f"{lineage} -> {resolved}"] += 1
        else:
            root.sequences.append(seq_id)
            report['placed_in_root'].setdefault(f"{lineage} -> {root.name}", 0)
            report['placed_in_root'][f"{lineage} -> {root.name}"] += 1

    return report


def report_assignments(report: Dict[str, object], total: int, max_dropped_fraction: float) -> int:
    """Print what happened to every sequence and fail if too many were lost."""
    placed = (report['assigned']
              + sum(report['placed_in_enclosing_lineage'].values())
              + sum(report['placed_in_root'].values()))
    dropped = sum(report['dropped_unknown_lineage'].values()) + len(report['without_assignment'])
    excluded = sum(report['outside_root'].values())

    print("\nSequence assignment:")
    print(f"  {report['assigned']} assigned to their own lineage")
    for label, count in sorted(report['placed_in_enclosing_lineage'].items(), key=lambda x: -x[1]):
        print(f"  {count} placed in enclosing lineage: {label}")
    for label, count in sorted(report['placed_in_root'].items(), key=lambda x: -x[1]):
        print(f"  {count} placed in the root lineage: {label}")
    for label, count in sorted(report['dropped_unknown_lineage'].items(), key=lambda x: -x[1]):
        print(f"  {count} DROPPED, lineage not in guide tree: {label}")
    if report['without_assignment']:
        examples = ', '.join(report['without_assignment'][:5])
        print(f"  {len(report['without_assignment'])} DROPPED, no lineage assignment "
              f"(e.g. {examples})")
    for label, count in sorted(report['outside_root'].items(), key=lambda x: -x[1]):
        print(f"  {count} excluded, outside the root lineage: {label}")
    if report['founders']:
        print(f"  ({report['founders']} founder sequences in the alignment, placed as founders)")
    print(f"  {placed} of {total} sequences placed, {dropped} dropped, {excluded} excluded "
          f"as outside the root lineage")
    # sequences outside the root lineage are excluded by choice, not lost
    total -= excluded

    # sequences that land in the root lineage are not lost, but a build in which most
    # sequences do is usually reading the wrong column or using a stale guide tree
    if total and report['assigned'] / total < 0.5:
        print(f"Warning: only {report['assigned'] / total:.1%} of sequences were assigned to "
              f"their own lineage. Check that --lineage-column names the column holding the "
              f"lineages of the guide tree and that the guide tree is up to date.")

    if total and dropped / total > max_dropped_fraction:
        print(f"Error: dropped {dropped / total:.1%} of sequences, more than the allowed "
              f"{max_dropped_fraction:.1%}. Check that --lineage-column names the column "
              f"holding the lineages of the guide tree.", file=sys.stderr)
        return 1
    return 0


def assign_founders_to_lineages(
    root: LineageNode,
    founders: Dict[str, SeqRecord],
    ignore_missing_founders: bool = False
):
    """
    Assign founder sequences to their respective lineage nodes.

    A lineage without a founder cannot be the root of its own subtree. Rather than
    discarding it with its sequences and descendants, it is folded into its parent:
    its sequences are built as part of the parent subtree and its children are
    re-attached to the parent. This keeps newly designated lineages, for which no
    founder has been defined yet, in the tree.
    """
    folded = []

    def assign(node):
        if node.name in founders:
            node.founder_seq = founders[node.name]
        children_to_keep = []
        # note that this allows a missing founder for the root node
        for child in node.children:
            if child.name in founders:
                assign(child)
                children_to_keep.append(child)
            elif ignore_missing_founders:
                assign(child)
                folded.append((child.name, node.name, len(child.sequences)))
                node.sequences.extend(child.sequences)
                children_to_keep.extend(child.children)
                for grandchild in child.children:
                    grandchild.parent = node
            else:
                raise ValueError(f"Founder sequence for lineage '{child.name}' not found.")
        node.children = children_to_keep

    assign(root)

    if folded:
        print("\nLineages without a founder sequence, folded into their parent:")
        for name, parent, n_seqs in folded:
            print(f"  {name} -> {parent} ({n_seqs} sequences kept)")
    return folded


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
        print(f"Rooting tree with {root}")
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
    if len(subtree_seqs) == 0:
        # a lineage can end up empty when its sequences were folded into a parent
        print(f"Lineage {node.name}: no sequences, no subtree built")
        return None
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


def write_guide_tree_lineages(
    seq_ids,
    assignments: Dict[str, str],
    parent: Dict[str, Optional[str]],
    guide_root: str,
    path: str,
):
    """
    Write the lineage of each sequence as a node of the full guide tree.

    This is the lineage call itself if it is a node of the guide tree, otherwise the
    nearest enclosing lineage, otherwise the root of the guide tree. Unlike the lineage
    a sequence was built into, it does not depend on lineages being folded for lack of
    a founder.
    """
    rows = [(seq_id, resolve_lineage(assignments.get(seq_id), parent, guide_root))
            for seq_id in seq_ids]
    pd.DataFrame(rows, columns=['seq_id', 'lineage']).to_csv(path, sep='\t', index=False)
    print(f"Wrote the lineages of {len(rows)} sequences to {path}")


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


def collapse_short_branches(tree, threshold: float) -> int:
    """
    Collapse internal branches shorter than threshold into polytomies, in place.

    A lineage whose founder sits at zero distance from the node where tips of its parent
    lineage branch off thereby gets a polytomy that contains those tips, rather than a
    zero length stem that a later step would collapse in an arbitrary way.
    """
    n_collapsed = 0
    for node in tree.find_clades(order='postorder'):
        children = []
        for child in node.clades:
            if child.clades and (child.branch_length or 0) < threshold:
                for grandchild in child.clades:
                    grandchild.branch_length = (grandchild.branch_length or 0) + (child.branch_length or 0)
                children.extend(child.clades)
                n_collapsed += 1
            else:
                children.append(child)
        node.clades = children
    return n_collapsed


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

    guide_parent, guide_root = read_guide_tree_parents(args.guide_tree, args.default_lineage)

    if args.founder_sequences:
        print("Loading founder sequences...")
        founders = load_sequences(args.founder_sequences)
        print(f"  Loaded {len(founders)} founder sequences")
    else:
        # founders are named after their lineage; recognise them against the full guide
        # tree, not only the part below the root lineage
        founders = {seq_id: seq for seq_id, seq in sequences.items() if seq_id in guide_parent}
        print(f"  gathered {len(founders)} founder sequences from alignment")

    print("Loading sequence-to-lineage assignments...")
    assignments = load_assignments(
        args.lineage_assignments,
        args.seq_id_column,
        args.lineage_column,
        args.default_lineage
    )
    print(f"  Loaded {len(assignments)} assignments")

    root_lineage = args.root_lineage
    if root_lineage == 'auto':
        root_lineage = automatic_root_lineage(
            [resolve_lineage(assignments.get(seq_id), guide_parent, guide_root)
             for seq_id in sequences if seq_id not in founders and seq_id in assignments],
            guide_parent)
        print(f"Root lineage chosen automatically as the lineage containing all sequences: {root_lineage}")
    elif root_lineage not in guide_parent:
        print(f"Error: root lineage '{root_lineage}' is not a lineage of the guide tree")
        return 1

    print("Parsing lineage guide tree...")
    lineage_root = parse_guide_tree(args.guide_tree, root_lineage, args.default_lineage)

    # founders of lineages outside the root lineage are not used
    lineages_below_root = {clade.name for clade in lineage_root.tree.find_clades() if clade.name}
    lineages_below_root.add(lineage_root.name)
    unused = sorted(seq_id for seq_id in founders if seq_id not in lineages_below_root)
    if unused:
        print(f"  ignoring {len(unused)} founders of lineages outside the root lineage: "
              f"{', '.join(unused)}")
        founders = {seq_id: seq for seq_id, seq in founders.items() if seq_id not in unused}
        sequences = {seq_id: seq for seq_id, seq in sequences.items() if seq_id not in unused}

    # Prepare lineage structure
    print("\nAssigning sequences to lineages...")
    report = assign_sequences_to_lineages(lineage_root, assignments, sequences,
                                          guide_parent, guide_root,
                                          args.unknown_lineage_policy, set(founders))
    if report_assignments(report, len(sequences) - len(set(founders) & set(sequences)),
                          args.max_dropped_fraction):
        return 1
    assign_founders_to_lineages(lineage_root, founders, args.ignore_missing_founders)

    if args.output_lineages:
        write_guide_tree_lineages([seq_id for seq_id in sequences if seq_id not in founders],
                                  assignments, guide_parent, guide_root, args.output_lineages)

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
        explicit_root,
        )

    if final_tree is None:
        print("Error: Failed to build final tree")
        return 1

    if not args.keep_founders:
        # grafting removes the founders of grafted lineages, but not that of the root
        # lineage or of lineages whose subtree is the founder alone
        for founder in [tip for tip in final_tree.get_terminals() if tip.name in founders]:
            final_tree.prune(founder)

    n_collapsed = collapse_short_branches(final_tree, args.collapse_threshold)
    print(f"\nCollapsed {n_collapsed} internal branches shorter than {args.collapse_threshold}")

    # Write output
    print(f"\nWriting final tree to {args.output_tree}...")
    Phylo.write(final_tree, args.output_tree, 'newick')

    print("Done!")
    return 0


if __name__ == '__main__':
    exit(main())
