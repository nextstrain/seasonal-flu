"""
Labelling the lineages of a tree from the lineages of its tips.

Every tip carries its lineage as a node of the guide tree. A lineage comprises the tips
of the lineage itself and of all lineages descending from it in the guide tree. Its
label goes on the branch leading to the MRCA of these tips, and every node in the
subtree below the MRCA is assigned to it. Lineages are processed from the root of the
guide tree downwards, so a sublineage overwrites its parent within its own subtree.

Founder sequences kept as tips are named after their lineage and belong to it, but do
not count as tips of their own. Three cases need care:

- A lineage without tips of its own, only of sublineages, is not labelled. In tidy tree
  it is represented only by its founder and has no node of its own.
- A lineage whose stem has zero length was collapsed into a polytomy, which may also
  hold subtrees made only of tips of ancestral lineages: they branch off at the founder.
  The lineage is labelled on the polytomy all the same and those subtrees keep their own
  label, unless an ancestral lineage is labelled on the same node, as it would then lose
  its label.
- If the subtree of the MRCA contains any other tips of other lineages, the lineage has
  no node of its own. Then its largest subtrees that contain only its own tips are
  labelled instead, so that it does not take over the rest of the subtree.

Unlike labelling by signature mutations, this does not depend on the mutations
reconstructed on the tree, only on the tips that are present, so it is robust to
pruning tips and collapsing branches after the tree was built.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .guide import GuideTree


def refine_calls(tip_lineage: Dict[str, str], calls: Dict[str, Optional[str]],
                 guide: GuideTree) -> Counter:
    """
    Take the call of a tip from the current nomenclature, in place, where it is a
    sublineage of the lineage of the tip, which catches up calls made with an older
    nomenclature. Returns the number of tips per (lineage, call).
    """
    refined = Counter()
    for name, lineage in tip_lineage.items():
        call = calls.get(name)
        if call and call != lineage and lineage in guide.ancestors(guide.parent.get(call)):
            tip_lineage[name] = call
            refined[(lineage, call)] += 1
    return refined


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


@dataclass
class Labelling:
    membership: Dict = field(default_factory=dict)          # node -> lineage
    labels: List[Tuple] = field(default_factory=list)       # (node, lineage) in labelling order
    split: Dict[str, Tuple[int, int, int]] = field(default_factory=dict)
    # lineage -> (tips, tips of other lineages under its MRCA, number of subtrees labelled)
    basal: Dict[str, Tuple[int, int]] = field(default_factory=dict)
    # lineage -> (tips, tips of ancestral lineages at its MRCA)
    missing: List[str] = field(default_factory=list)        # tips without a lineage

    def node_data(self, membership_name: str = 'subclade', label_name: str = 'Subclade') -> Dict:
        """Node data in the format of augur clades; a branch keeps the most specific label."""
        branch_labels = {node.name: lineage for node, lineage in self.labels}
        return {
            'nodes': {node.name: {membership_name: lineage}
                      for node, lineage in self.membership.items()},
            'branches': {name: {'labels': {label_name: lineage}}
                         for name, lineage in branch_labels.items()},
        }

    def print(self):
        print("Lineages without a node of their own, labelled on their largest pure subtrees:")
        for lineage, (n, n_other, n_subtrees) in self.split.items():
            print(f"  {lineage}: {n} tips in {n_subtrees} subtrees, "
                  f"{n_other} tips of other lineages under their MRCA")
        if not self.split:
            print("  none")
        print("Lineages labelled on a polytomy that holds tips of ancestral lineages:")
        for lineage, (n, n_basal) in self.basal.items():
            print(f"  {lineage}: {n} tips, {n_basal} tips of ancestral lineages")
        if not self.basal:
            print("  none")


def label_lineages(tree, tip_lineage: Dict[str, str], guide: GuideTree,
                   default_lineage: str = 'unassigned') -> Labelling:
    """Label the lineages of the tree from the lineages of its tips."""
    up, depth = {}, {tree.root: 0}
    for node in tree.find_clades(order='preorder'):
        for child in node.clades:
            up[child] = node
            depth[child] = depth[node] + 1

    tips = tree.get_terminals()
    tip_lineage = dict(tip_lineage)
    # founder sequences kept as tips are named after their lineage and belong to it
    founders = {tip.name for tip in tips if tip.name not in tip_lineage and tip.name in guide}
    tip_lineage.update({name: name for name in founders})
    result = Labelling(missing=[tip.name for tip in tips if tip.name not in tip_lineage])
    lineage_of = {tip: tip_lineage.get(tip.name, default_lineage) for tip in tips}

    # the tips of each lineage, including those of its sublineages
    members = defaultdict(list)
    for tip in tips:
        for lineage in guide.ancestors(lineage_of[tip]):
            members[lineage].append(tip)

    # number of tips below each node, for telling whether a subtree is pure
    n_tips = {}
    for node in tree.find_clades(order='postorder'):
        n_tips[node] = 1 if node.is_terminal() else sum(n_tips[c] for c in node.clades)

    result.membership = {node: default_lineage for node in tree.find_clades()}
    labelled = set()
    # lineages called for at least one tip; one represented only by its founder in tidy
    # tree has no node of its own to label
    called = {lineage_of[tip] for tip in tips if tip.name not in founders}
    for lineage in guide.preorder:
        if lineage == default_lineage or lineage not in called:
            continue
        own = members[lineage]
        root = mrca(own, up, depth)
        ancestors = set(guide.ancestors(guide.parent.get(lineage)))
        # an ancestral lineage labelled on the same node would lose its label
        at_founder = set() if root in labelled else {
            child for child in root.clades
            if all(lineage_of[tip] in ancestors for tip in child.get_terminals())}
        n_at_founder = sum(n_tips[child] for child in at_founder)
        if n_tips[root] == len(own) + n_at_founder:
            roots = [root]
            if at_founder:
                result.basal[lineage] = (len(own), n_at_founder)
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
            result.split[lineage] = (len(own), n_tips[root] - len(own), len(roots))
        for subtree_root in roots:
            result.labels.append((subtree_root, lineage))
            labelled.add(subtree_root)
            for child in subtree_root.clades:
                if child not in at_founder:
                    for node in child.find_clades():
                        result.membership[node] = lineage
            result.membership[subtree_root] = lineage
    return result


def compare(node_data: Dict, other: Dict, tip_names, membership_name: str) -> Counter:
    """Tips whose lineage differs from another labelling, by (ours, theirs)."""
    disagreements = Counter()
    for name in tip_names:
        ours = node_data['nodes'][name][membership_name]
        theirs = other.get(name, {}).get(membership_name)
        if ours != theirs:
            disagreements[(ours, theirs)] += 1
    return disagreements
