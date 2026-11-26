from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tidytree.assign import AssignmentReport, assign_founders, assign_sequences


def test_resolve(guide):
    assert guide.resolve('A.1') == 'A.1'
    # finer than the guide tree or designated since: the enclosing lineage by name
    assert guide.resolve('A.1.1.3') == 'A.1.1'
    assert guide.resolve('A.7') == 'A'
    # nothing enclosing, or no call at all: the root of the guide tree
    assert guide.resolve('Z.1') == 'unassigned'
    assert guide.resolve(None) == 'unassigned'


def test_automatic_root_lineage(guide):
    assert guide.automatic_root_lineage(['A.1.1', 'A.2']) == 'A'
    assert guide.automatic_root_lineage(['A.1.1']) == 'A.1.1'
    assert guide.automatic_root_lineage(['A.1', 'B']) == 'unassigned'


def test_lineage_tree(guide):
    root = guide.lineage_tree('A')
    assert [node.name for node in root.preorder()] == ['A', 'A.1', 'A.1.1', 'A.2']
    assert [node.name for node in root.postorder()][-1] == 'A'
    assert set(guide.lineages_below('A')) == {'A', 'A.1', 'A.1.1', 'A.2'}


def assignment(guide, policy, calls):
    root = guide.lineage_tree('A')
    report = assign_sequences(root, list(calls), calls, guide, policy, founder_ids={'A.2'})
    return {node.name: node.sequences for node in root.preorder()}, report


CALLS = {'s1': 'A.1', 's2': 'A.1.1.4', 's3': 'A.9', 's4': 'B', 'A.2': 'A.2'}


def test_assign_policy_parent(guide):
    built, report = assignment(guide, 'parent', CALLS)
    assert built == {'A': ['s3'], 'A.1': ['s1'], 'A.1.1': ['s2'], 'A.2': []}
    assert report.assigned == 1 and report.founders == 1
    assert report.placed_in_enclosing_lineage == {'A.1.1.4 -> A.1.1': 1}
    # A.9 resolves to the root lineage itself, which counts as placed in the root
    assert report.placed_in_root == {'A.9 -> A': 1}
    assert report.outside_root == {'B': 1}
    assert (report.placed, report.dropped, report.excluded) == (3, 0, 1)


def test_assign_policy_root_and_drop(guide):
    built, _ = assignment(guide, 'root', CALLS)
    assert built['A'] == ['s2', 's3']
    built, report = assignment(guide, 'drop', CALLS)
    assert built['A.1.1'] == [] and report.dropped == 2


def test_report_check():
    report = AssignmentReport(assigned=90, without_assignment=['x'] * 10)
    assert report.check(100, max_dropped_fraction=0.2)
    assert not report.check(100, max_dropped_fraction=0.05)


def seq(name):
    return SeqRecord(Seq('ACGT'), id=name)


def test_assign_founders_folds_lineages_without_founder(guide):
    root = guide.lineage_tree('A')
    nodes = {node.name: node for node in root.preorder()}
    nodes['A.1'].sequences = ['s1']
    nodes['A.1.1'].sequences = ['s2']
    folded = assign_founders(root, {name: seq(name) for name in ['A', 'A.1.1', 'A.2']},
                             ignore_missing_founders=True)
    # A.1 has no founder: its sequences go to A and A.1.1 becomes a child of A
    assert folded == [('A.1', 'A', 1)]
    assert [child.name for child in root.children] == ['A.1.1', 'A.2']
    assert root.sequences == ['s1'] and nodes['A.1.1'].parent is root
