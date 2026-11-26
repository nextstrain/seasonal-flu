"""Assigning sequences and founders to the lineages of the guide tree."""
import sys
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set

import pandas as pd

from .guide import GuideTree, LineageNode


def load_assignments(path: str, seq_id_column: str = 'seq_id', lineage_column: str = 'lineage',
                     default_lineage: str = 'unassigned') -> Dict[str, str]:
    """Lineage call of each sequence from a TSV table; missing calls get the default lineage."""
    df = pd.read_csv(path, sep='\t', dtype=str)
    for column in (seq_id_column, lineage_column):
        if column not in df.columns:
            raise ValueError(f"Column '{column}' not found in assignments file. "
                             f"Available columns: {', '.join(df.columns)}")
    lineages = df[lineage_column].fillna(default_lineage)
    lineages[~lineages.astype(bool)] = default_lineage
    return dict(zip(df[seq_id_column], lineages))


def oldest_sequence(seq_ids: Iterable[str], metadata_path: str, id_column: str = 'strain',
                    date_column: str = 'date') -> Optional[str]:
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
    return min(dated)[1] if dated else None


@dataclass
class AssignmentReport:
    """What happened to each sequence of the alignment, which ends up in exactly one category."""
    assigned: int = 0
    founders: int = 0
    placed_in_enclosing_lineage: Counter = field(default_factory=Counter)  # "call -> lineage"
    placed_in_root: Counter = field(default_factory=Counter)               # "call -> root"
    dropped_unknown_lineage: Counter = field(default_factory=Counter)      # only with policy drop
    outside_root: Counter = field(default_factory=Counter)                 # excluded on purpose
    without_assignment: List[str] = field(default_factory=list)            # no lineage call

    @property
    def placed(self) -> int:
        return (self.assigned + sum(self.placed_in_enclosing_lineage.values())
                + sum(self.placed_in_root.values()))

    @property
    def dropped(self) -> int:
        return sum(self.dropped_unknown_lineage.values()) + len(self.without_assignment)

    @property
    def excluded(self) -> int:
        return sum(self.outside_root.values())

    def print(self, total: int):
        print("\nSequence assignment:")
        print(f"  {self.assigned} assigned to their own lineage")
        for label, count in self.placed_in_enclosing_lineage.most_common():
            print(f"  {count} placed in enclosing lineage: {label}")
        for label, count in self.placed_in_root.most_common():
            print(f"  {count} placed in the root lineage: {label}")
        for label, count in self.dropped_unknown_lineage.most_common():
            print(f"  {count} DROPPED, lineage not in guide tree: {label}")
        if self.without_assignment:
            print(f"  {len(self.without_assignment)} DROPPED, no lineage assignment "
                  f"(e.g. {', '.join(self.without_assignment[:5])})")
        for label, count in self.outside_root.most_common():
            print(f"  {count} excluded, outside the root lineage: {label}")
        if self.founders:
            print(f"  ({self.founders} founder sequences in the alignment, placed as founders)")
        print(f"  {self.placed} of {total} sequences placed, {self.dropped} dropped, "
              f"{self.excluded} excluded as outside the root lineage")

    def check(self, total: int, max_dropped_fraction: float) -> bool:
        """Whether few enough sequences were dropped; warns if few found their own lineage."""
        # sequences outside the root lineage are excluded by choice, not lost
        total -= self.excluded
        if not total:
            return True
        # sequences that land in the root lineage are not lost, but a build in which most
        # sequences do is usually reading the wrong column or using a stale guide tree
        if self.assigned / total < 0.5:
            print(f"Warning: only {self.assigned / total:.1%} of sequences were assigned to "
                  f"their own lineage. Check that --lineage-column names the column holding the "
                  f"lineages of the guide tree and that the guide tree is up to date.")
        if self.dropped / total > max_dropped_fraction:
            print(f"Error: dropped {self.dropped / total:.1%} of sequences, more than the allowed "
                  f"{max_dropped_fraction:.1%}. Check that --lineage-column names the column "
                  f"holding the lineages of the guide tree.", file=sys.stderr)
            return False
        return True


def assign_sequences(root: LineageNode, seq_ids: Iterable[str], assignments: Dict[str, str],
                     guide: GuideTree, unknown_lineage_policy: str = 'parent',
                     founder_ids: Optional[Set[str]] = None) -> AssignmentReport:
    """
    Assign sequences to the lineages below root by their lineage calls.

    A call that is not a lineage of the guide tree is placed in its nearest enclosing
    lineage, in the root lineage or dropped, according to unknown_lineage_policy. Calls
    that resolve to a lineage outside root are excluded.
    """
    lineages = {node.name: node for node in root.preorder()}
    founder_ids = founder_ids or set()
    report = AssignmentReport()
    for seq_id in seq_ids:
        if seq_id in founder_ids:
            # founders are placed by assign_founders, not by their lineage call
            report.founders += 1
            continue
        lineage = assignments.get(seq_id)
        if lineage is None:
            report.without_assignment.append(seq_id)
            continue
        if lineage in lineages:
            lineages[lineage].sequences.append(seq_id)
            report.assigned += 1
            continue
        resolved = guide.resolve(lineage)
        if resolved not in lineages:
            report.outside_root[lineage] += 1
        elif unknown_lineage_policy == 'drop':
            report.dropped_unknown_lineage[lineage] += 1
        elif unknown_lineage_policy == 'parent' and resolved != root.name:
            lineages[resolved].sequences.append(seq_id)
            report.placed_in_enclosing_lineage[f"{lineage} -> {resolved}"] += 1
        else:
            root.sequences.append(seq_id)
            report.placed_in_root[f"{lineage} -> {root.name}"] += 1
    return report


def assign_founders(root: LineageNode, founders: Dict, ignore_missing_founders: bool = False):
    """
    Give each lineage its founder sequence.

    A lineage without a founder cannot be the root of its own subtree. Rather than
    discarding it with its sequences and descendants, it is folded into its parent: its
    sequences are built as part of the parent subtree and its children are re-attached
    to the parent. The root lineage may lack a founder.

    Returns (lineage, parent, number of sequences) for each folded lineage.
    """
    folded = []
    for node in root.postorder():
        node.founder_seq = founders.get(node.name)
        children = []
        for child in node.children:
            if child.founder_seq is not None:
                children.append(child)
            elif ignore_missing_founders:
                folded.append((child.name, node.name, len(child.sequences)))
                node.sequences.extend(child.sequences)
                for grandchild in child.children:
                    grandchild.parent = node
                children.extend(child.children)
            else:
                raise ValueError(f"Founder sequence for lineage '{child.name}' not found.")
        node.children = children
    return folded


def write_lineages(seq_ids: Iterable[str], assignments: Dict[str, str], guide: GuideTree,
                   path: str) -> int:
    """
    Write the lineage of each sequence as a node of the full guide tree.

    This is the lineage call itself if it is a node of the guide tree, otherwise the
    nearest enclosing lineage, otherwise the root of the guide tree. Unlike the lineage a
    sequence was built into, it does not depend on lineages being folded for lack of a
    founder.
    """
    rows = [(seq_id, guide.resolve(assignments.get(seq_id))) for seq_id in seq_ids]
    pd.DataFrame(rows, columns=['seq_id', 'lineage']).to_csv(path, sep='\t', index=False)
    return len(rows)
