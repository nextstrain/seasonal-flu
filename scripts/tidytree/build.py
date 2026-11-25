"""Building the subtree of each lineage and stitching them along the guide tree."""
import io
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional

from Bio import Phylo, SeqIO
from Bio.SeqRecord import SeqRecord

from .guide import LineageNode


def p_distance(a: SeqRecord, b: SeqRecord) -> float:
    """Fraction of sites with an unambiguous nucleotide in both sequences at which they differ."""
    compared = differ = 0
    for x, y in zip(str(a.seq).upper(), str(b.seq).upper()):
        if x in 'ACGT' and y in 'ACGT':
            compared += 1
            differ += x != y
    return differ / compared if compared else 0.0


def run_iqtree(sequences: List[SeqRecord], model: str, root: str, iqtree_args: str = '',
               threads: int = 1) -> Optional[Phylo.BaseTree.Tree]:
    """
    Tree of the aligned sequences by IQ-TREE, rooted on the sequence named root.

    Uses the IQ-TREE runner of augur tree, which escapes the characters IQ-TREE alters in
    sequence names.
    """
    from augur.tree import build_iqtree

    with tempfile.TemporaryDirectory() as tmpdir:
        alignment = Path(tmpdir) / 'alignment.fasta'
        SeqIO.write(sequences, alignment, 'fasta')
        tree, _ = build_iqtree(str(alignment), str(Path(tmpdir) / 'tree.nwk'),
                               substitution_model=model, nthreads=threads,
                               tree_builder_args=iqtree_args)
    if tree is not None:
        tree.root_with_outgroup(root)
    return tree


def subtree_sequences(node: LineageNode, sequences: Dict[str, SeqRecord]) -> List[SeqRecord]:
    """The sequences of a lineage, its founder and the founders of its children."""
    seqs = [sequences[seq_id] for seq_id in node.sequences if seq_id in sequences]
    if node.founder_seq is not None:
        seqs.append(node.founder_seq)
    seqs.extend(child.founder_seq for child in node.children if child.founder_seq is not None)
    return seqs


def build_subtree(node: LineageNode, sequences: Dict[str, SeqRecord], model: str,
                  iqtree_args: str = '', explicit_root: Optional[str] = None,
                  threads: int = 1) -> Optional[Phylo.BaseTree.Tree]:
    """
    Subtree of a lineage, rooted on its founder or on explicit_root.

    IQ-TREE needs at least three sequences. A founder with a single other sequence gets
    the distance between the two as branch length; a missing one would count as zero
    and the branch be collapsed.
    """
    seqs = subtree_sequences(node, sequences)
    if not seqs:
        # a lineage can end up empty when its sequences were folded into a parent
        return None
    if len(seqs) >= 3:
        return run_iqtree(seqs, model, explicit_root or node.name, iqtree_args, threads)
    if len(seqs) == 1:
        newick = f"{seqs[0].id};"
    elif node.founder_seq is not None:
        other = next(seq for seq in seqs if seq.id != node.founder_seq.id)
        newick = (f"({node.founder_seq.id}:0,{other.id}:{p_distance(node.founder_seq, other)})"
                  f"{node.name};")
    elif node.name.startswith('internal_'):
        newick = f"({','.join(seq.id for seq in seqs)});"
    else:
        newick = f"({','.join(seq.id for seq in seqs)}){node.name};"
    return Phylo.read(io.StringIO(newick), 'newick')


def find_clade(tree, name: str):
    """The first clade of the tree with the given name, or None."""
    return next(tree.find_clades(name=name), None)


def graft_subtree(parent_tree, child_tree, founder: str, keep_founders: bool = False) -> bool:
    """
    Graft the subtree of a lineage onto the subtree of its parent at its founder.

    The founder is a tip of the parent subtree and the root of the child subtree, or a
    child of its root. The tip is replaced by the child subtree and renamed
    '<founder>_internal'; the founder itself is kept as a tip only with keep_founders.
    Returns whether grafting took place.
    """
    parent_node = find_clade(parent_tree, founder)
    child_node = find_clade(child_tree, founder)
    if parent_node is None or child_node is None:
        print(f"Warning: Could not find {founder} in the "
              f"{'parent' if parent_node is None else 'child'} tree")
        return False
    assert parent_node.is_terminal(), "Connection point in parent tree must be a leaf"

    root = child_tree.root
    if child_node is root:
        if not root.clades:
            # the founder alone, which is the tip of the parent tree already
            return False
        clades = list(root.clades)
        if keep_founders:
            child_node.clades = []
            clades.append(child_node)
    else:
        clades = [c for c in root.clades if c.name != founder or keep_founders]
    parent_node.clades = clades
    if parent_node.name:
        parent_node.name += "_internal"
    return True


def stitch_subtrees(root: LineageNode, sequences: Dict[str, SeqRecord], model: str,
                    iqtree_args: str = '', threads: int = 1, keep_founders: bool = False,
                    explicit_root: Optional[str] = None) -> Optional[Phylo.BaseTree.Tree]:
    """
    Build the subtree of every lineage and graft each onto its parent at its founder.

    The subtrees are independent, so up to `threads` IQ-TREE runs take place at the
    same time; the grafting then proceeds from the leaves of the guide tree upwards.
    """
    lineages = root.postorder()
    with ThreadPoolExecutor(max_workers=max(threads, 1)) as pool:
        futures = {
            node: pool.submit(build_subtree, node, sequences, model, iqtree_args,
                              explicit_root if node is root else None)
            for node in lineages
        }
        for node in lineages:
            node.subtree = futures[node].result()
            n = len(subtree_sequences(node, sequences))
            print(f"Lineage {node.name}: subtree of {n} sequences"
                  + ("" if node.subtree is not None or n == 0 else ", FAILED"))

    for node in lineages:
        if node.subtree is None:
            continue
        for child in node.children:
            if child.subtree is not None and child.founder_seq is not None:
                print(f"Grafting {child.name} subtree onto {node.name} at {child.founder_seq.id}")
                graft_subtree(node.subtree, child.subtree, child.founder_seq.id, keep_founders)
            else:
                print(f"Skipping grafting for {child.name} onto {node.name} (missing subtree or founder)")
    return root.subtree


def collapse_short_branches(tree, threshold: float) -> int:
    """
    Collapse internal branches shorter than threshold into polytomies, in place, and
    merge nodes with a single child into it.

    A lineage whose founder sits at zero distance from the node where tips of its parent
    lineage branch off thereby gets a polytomy that contains those tips, rather than a
    zero length stem that a later step would collapse in an arbitrary way. Single child
    nodes are left where grafting removed the founder of a lineage. A branch without a
    length is kept, as its length is unknown.
    """
    n_collapsed = 0
    for node in tree.find_clades(order='postorder'):
        children = []
        for child in node.clades:
            short = child.branch_length is not None and child.branch_length < threshold
            if child.clades and (short or len(child.clades) == 1):
                for grandchild in child.clades:
                    grandchild.branch_length = (grandchild.branch_length or 0) + (child.branch_length or 0)
                children.extend(child.clades)
                n_collapsed += 1
            else:
                children.append(child)
        node.clades = children
    while len(tree.root.clades) == 1 and tree.root.clades[0].clades:
        tree.root = tree.root.clades[0]
        n_collapsed += 1
    return n_collapsed
