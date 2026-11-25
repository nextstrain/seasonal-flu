from conftest import tree
from tidytree.labels import label_lineages, refine_calls


def labels_of(newick, tip_lineage, guide):
    t = tree(newick)
    result = label_lineages(t, tip_lineage, guide)
    node_data = result.node_data()
    branch = {name: data['labels']['Subclade'] for name, data in node_data['branches'].items()}
    member = {name: data['subclade'] for name, data in node_data['nodes'].items()}
    return branch, member, result


def test_nested_lineages(guide):
    branch, member, result = labels_of(
        '((a1,a2,(b1,b2)nB)nA,(c1,c2)nC)root;',
        {'a1': 'A', 'a2': 'A', 'b1': 'A.1', 'b2': 'A.1', 'c1': 'B', 'c2': 'B'}, guide)
    assert branch == {'nA': 'A', 'nB': 'A.1', 'nC': 'B'}
    assert member['a1'] == 'A' and member['b1'] == member['nB'] == 'A.1'
    assert member['root'] == 'unassigned'
    assert not result.split and not result.basal


def test_single_tip_lineage(guide):
    branch, member, _ = labels_of(
        '((a1,a2,b1)nA,c1)root;', {'a1': 'A', 'a2': 'A', 'b1': 'A.2', 'c1': 'B'}, guide)
    # the label of a single tip lineage goes on the tip itself
    assert branch['b1'] == 'A.2' and branch['c1'] == 'B' and branch['nA'] == 'A'


def test_misplaced_tip_splits_lineage(guide):
    # a B tip within A: A has no node of its own and is labelled on its pure subtrees
    branch, member, result = labels_of(
        '(((a1,a2)x,(a3,c1)y)nA,c2)root;',
        {'a1': 'A', 'a2': 'A', 'a3': 'A', 'c1': 'B', 'c2': 'B'}, guide)
    assert result.split['A'] == (3, 1, 2)
    assert branch['x'] == 'A' and branch['a3'] == 'A'
    assert member['c1'] == 'B' and member['nA'] == 'unassigned'


def test_collapsed_stem_keeps_ancestral_tips_at_founder(guide):
    # the stem of A.1 was collapsed: tips of A branch off in the polytomy of A.1
    branch, member, result = labels_of(
        '((a1,(a2,(p1,p2,(p3,p4)q,a3)poly)z)nA,c1)root;',
        {'a1': 'A', 'a2': 'A', 'a3': 'A', 'p1': 'A.1', 'p2': 'A.1', 'p3': 'A.1', 'p4': 'A.1',
         'c1': 'B'}, guide)
    assert branch['poly'] == 'A.1'
    assert result.basal['A.1'] == (4, 1)
    # the ancestral tip keeps its own lineage
    assert member['a3'] == 'A' and member['q'] == 'A.1'


def test_no_label_on_a_node_an_ancestor_is_labelled_on(guide):
    # A and A.1 share their MRCA, which holds a tip of A: A keeps the node, A.1 is
    # labelled on its pure subtrees
    branch, member, result = labels_of(
        '((a1,(p1,p2)u,(p3,p4)v)nA,c1)root;',
        {'a1': 'A', 'p1': 'A.1', 'p2': 'A.1', 'p3': 'A.1', 'p4': 'A.1', 'c1': 'B'}, guide)
    assert branch['nA'] == 'A' and branch['u'] == branch['v'] == 'A.1'
    assert member['a1'] == 'A'


def test_lineage_represented_only_by_its_founder_is_not_labelled(guide):
    # A.1 has no tips of its own, only of its sublineage A.1.1, and its founder
    branch, member, _ = labels_of(
        '((a1,a2,(A.1,(d1,d2)n11)n1)nA,c1)root;',
        {'a1': 'A', 'a2': 'A', 'd1': 'A.1.1', 'd2': 'A.1.1', 'c1': 'B'}, guide)
    assert 'A.1' not in branch.values()
    assert branch['n11'] == 'A.1.1'
    # the founder tip belongs to its lineage
    assert member['A.1'] == 'A'


def test_founder_tips_belong_to_their_lineage(guide):
    branch, member, result = labels_of(
        '((a1,a2,(A.1,b1,b2)n1)nA,c1)root;',
        {'a1': 'A', 'a2': 'A', 'b1': 'A.1', 'b2': 'A.1', 'c1': 'B'}, guide)
    assert branch['n1'] == 'A.1' and member['A.1'] == 'A.1'
    assert not result.missing


def test_refine_calls(guide):
    tip_lineage = {'s1': 'A', 's2': 'A.1', 's3': 'A', 's4': 'unassigned'}
    calls = {'s1': 'A.1.1', 's2': 'A.2', 's3': 'A', 's4': 'B'}
    refined = refine_calls(tip_lineage, calls, guide)
    # a sublineage is taken, a sibling is not
    assert tip_lineage == {'s1': 'A.1.1', 's2': 'A.1', 's3': 'A', 's4': 'B'}
    assert refined == {('A', 'A.1.1'): 1, ('unassigned', 'B'): 1}
