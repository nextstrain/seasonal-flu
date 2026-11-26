"""The lineage guide tree: the hierarchy of a lineage nomenclature as a tree."""
from typing import Dict, Iterable, List, Optional

from Bio import Phylo


class LineageNode:
    """A lineage of the guide tree, with what tidy tree builds for it."""

    def __init__(self, name: str, parent: Optional['LineageNode'] = None):
        self.name = name
        self.parent = parent
        self.children: List['LineageNode'] = []
        self.sequences: List[str] = []  # IDs of the sequences built into this lineage
        self.founder_seq = None         # founder sequence, the root of its subtree
        self.subtree = None             # phylogeny of the lineage

    def preorder(self) -> List['LineageNode']:
        """This lineage and all below it, parents before children."""
        nodes, stack = [], [self]
        while stack:
            node = stack.pop()
            nodes.append(node)
            stack.extend(reversed(node.children))
        return nodes

    def postorder(self) -> List['LineageNode']:
        """This lineage and all below it, children before parents."""
        return self.preorder()[::-1]


class GuideTree:
    """
    Parent of each lineage of the full guide tree and the name of its root.

    An unnamed root takes the default lineage, which is what sequences without a
    lineage call are assigned to.
    """

    def __init__(self, path: str, default_lineage: str = 'unassigned'):
        self.tree = Phylo.read(path, 'newick')
        self.tree.root.name = self.tree.root.name or default_lineage
        self.root = self.tree.root.name
        self.parent: Dict[str, Optional[str]] = {self.root: None}
        self.preorder: List[str] = []  # lineages, parents before children
        for clade in self.tree.find_clades(order='preorder'):
            if clade.name:
                self.preorder.append(clade.name)
            for child in clade.clades:
                if child.name:
                    self.parent[child.name] = clade.name

    def __contains__(self, lineage) -> bool:
        return lineage in self.parent

    def ancestors(self, lineage: Optional[str]) -> List[str]:
        """
        The lineage and its ancestors, from the lineage up to the root. A lineage that is
        not in the guide tree has no ancestors.
        """
        path = []
        while lineage is not None:
            path.append(lineage)
            lineage = self.parent.get(lineage)
        return path

    def enclosing_lineage(self, lineage: str) -> Optional[str]:
        """
        Nearest ancestor of a lineage that is not in the guide tree, by its name.

        Lineages added to the nomenclature after the guide tree was made, or more finely
        resolved than the guide tree, are placed in the lineage they descend from, e.g.
        'K.1.2' in 'K.1' or else in 'K'.
        """
        parts = lineage.split('.')
        for i in range(len(parts) - 1, 0, -1):
            candidate = '.'.join(parts[:i])
            if candidate in self.parent:
                return candidate
        return None

    def resolve(self, lineage: Optional[str]) -> str:
        """
        The lineage of the guide tree a lineage call corresponds to: the call itself,
        else its nearest enclosing lineage, else the root of the guide tree.
        """
        if lineage in self.parent:
            return lineage
        if isinstance(lineage, str):
            return self.enclosing_lineage(lineage) or self.root
        return self.root

    def automatic_root_lineage(self, lineages: Iterable[str]) -> str:
        """The most specific lineage of the guide tree containing all given lineages."""
        common = None
        for lineage in set(lineages):
            path = self.ancestors(lineage)[::-1]
            if common is None:
                common = path
            else:
                n = 0
                while n < min(len(common), len(path)) and common[n] == path[n]:
                    n += 1
                common = common[:n]
        return common[-1] if common else self.root

    def lineages_below(self, lineage: str) -> List[str]:
        """The lineage and all lineages below it."""
        start = next(self.tree.find_clades(name=lineage))
        return [clade.name for clade in start.find_clades() if clade.name]

    def lineage_tree(self, root_lineage: str) -> LineageNode:
        """The part of the guide tree below root_lineage as LineageNodes."""
        start = next(self.tree.find_clades(name=root_lineage))
        root = LineageNode(start.name)
        stack = [(start, root)]
        while stack:
            clade, node = stack.pop()
            for child in clade.clades:
                child_node = LineageNode(child.name or f"internal_{id(child)}", node)
                node.children.append(child_node)
                stack.append((child, child_node))
        return root
