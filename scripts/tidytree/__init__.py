"""
Tidy tree: a phylogeny built along the hierarchy of a lineage nomenclature.

A subtree is built for each lineage of the guide tree, from the sequences called for it
and the founder sequences of the lineage and its child lineages, and rooted on its own
founder. The subtrees are then grafted onto each other at the founders. This gives more
stable trees than a single global tree, especially with little historical context.

- guide: the lineage guide tree and the resolution of lineage calls
- assign: assigning sequences and founders to lineages
- build: building, grafting and tidying the subtrees
- labels: labelling the lineages of a tree from the lineages of its tips
"""
