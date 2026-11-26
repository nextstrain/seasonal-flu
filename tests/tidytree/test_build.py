import shutil
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from conftest import SCRIPTS, tree
from tidytree.build import build_subtree, collapse_short_branches, graft_subtree, p_distance
from tidytree.guide import LineageNode

needs_iqtree = pytest.mark.skipif(
    not any(shutil.which(x) for x in ['iqtree3', 'iqtree2', 'iqtree']), reason='needs IQ-TREE')


def test_p_distance_ignores_ambiguous_sites():
    a, b = SeqRecord(Seq('ACGTN-AC')), SeqRecord(Seq('ACGAAAAT'))
    assert p_distance(a, b) == pytest.approx(2 / 6)


def test_collapse_short_branches():
    t = tree('(((a:0.01,b:0.02)u:0.003)G:0.002,(c:0.01,d:0.01)n,(e:0.1,f:0.2)z:0.000001)r;')
    assert collapse_short_branches(t, 1e-4) == 2
    # the single child node G merges into u, the short branch z is collapsed, the
    # branch n without a length is kept
    children = {c.name or c.get_terminals()[0].name: c for c in t.root.clades}
    assert set(children) == {'u', 'n', 'e', 'f'}
    assert children['u'].branch_length == pytest.approx(0.005)
    assert children['n'].branch_length is None
    assert children['e'].branch_length == pytest.approx(0.100001)


def test_collapse_moves_single_child_root_down():
    t = tree('((a:0.01,b:0.02)x:0.005)r;')
    collapse_short_branches(t, 1e-4)
    assert t.root.name == 'x' and len(t.root.clades) == 2


@pytest.mark.parametrize('keep_founders', [False, True])
def test_graft_at_founder(keep_founders):
    parent = tree('((p1:0.1,A.1:0.2):0.1,p2:0.3);')
    # a subtree as IQ-TREE gives it, rooted on the founder
    child = tree('(A.1:0.0,(c1:0.1,c2:0.1):0.05);')
    assert graft_subtree(parent, child, 'A.1', keep_founders)
    names = {tip.name for tip in parent.get_terminals()}
    assert {'c1', 'c2'} <= names
    assert ('A.1' in names) == keep_founders
    assert next(parent.find_clades(name='A.1_internal')) is not None


def test_graft_founder_alone_leaves_the_tip():
    parent = tree('(p1:0.1,A.1:0.2);')
    assert not graft_subtree(parent, tree('A.1;'), 'A.1')
    assert {tip.name for tip in parent.get_terminals()} == {'p1', 'A.1'}


def test_small_subtree_gets_branch_lengths():
    node = LineageNode('A')
    node.founder_seq = SeqRecord(Seq('AAAAAAAAAA'), id='A')
    node.sequences = ['s1']
    seqs = {'s1': SeqRecord(Seq('AAAAAAAAAT'), id='s1')}
    subtree = build_subtree(node, seqs, 'GTR+G')
    assert subtree.root.name == 'A'
    lengths = {tip.name: tip.branch_length for tip in subtree.get_terminals()}
    assert lengths == {'A': 0.0, 's1': pytest.approx(0.1)}


@needs_iqtree
@pytest.mark.parametrize('threads', [1, 3])
def test_tidy_tree_end_to_end(synthetic, threads):
    out = synthetic / 'tree.nwk'
    result = subprocess.run([
        sys.executable, str(SCRIPTS / 'tidy_tree.py'),
        '--alignment', str(synthetic / 'alignment.fasta'),
        '--guide-tree', str(synthetic / 'guide.nwk'),
        '--lineage-assignments', str(synthetic / 'lineages.tsv'),
        '--seq-id-column', 'strain', '--lineage-column', 'clade',
        '--metadata', str(synthetic / 'metadata.tsv'),
        '--output-lineages', str(synthetic / 'tidy_lineages.tsv'),
        '--iqtree-args=-n 2 -czb -seed 1', '--threads', str(threads),
        '--output-tree', str(out),
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'Root lineage chosen automatically as the lineage containing all sequences: unassigned' in result.stdout
    assert 'rooting on its oldest sequence old/1968' in result.stdout
    assert '21 of 21 sequences placed, 0 dropped' in result.stdout

    t = Phylo.read(out, 'newick')
    tips = {tip.name for tip in t.get_terminals()}
    # every sequence is in the tree, the founders are not
    assert len(tips) == 21 and not tips & {'A', 'B', 'A.1', 'A.2', 'A.1.1'}
    assert 'old/1968' in [c.name for c in t.root.clades]
    # each lineage is a clade of the tree
    for lineage in ['A', 'A.1', 'A.1.1', 'A.2', 'B']:
        members = [x for x in tips if x.split('/')[0] == lineage
                   or x.split('/')[0].startswith(lineage + '.')]
        assert t.common_ancestor(members).count_terminals() == len(members), lineage
    # no single child nodes, tree ladderized
    for node in t.get_nonterminals():
        assert len(node.clades) > 1
        sizes = [c.count_terminals() for c in node.clades]
        assert sizes == sorted(sizes)

    # and the subclades are labelled on these clades
    labels = synthetic / 'labels.json'
    result = subprocess.run([
        sys.executable, str(SCRIPTS / 'assign_subclades_from_tree.py'),
        '--tree', str(out), '--lineages', str(synthetic / 'tidy_lineages.tsv'),
        '--guide-tree', str(synthetic / 'guide.nwk'), '--output-node-data', str(labels),
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    import json
    node_data = json.load(open(labels))
    assert sorted(b['labels']['Subclade'] for b in node_data['branches'].values()) == \
        ['A', 'A.1', 'A.1.1', 'A.2', 'B']
    assert node_data['nodes']['A.1.1/seq0']['subclade'] == 'A.1.1'
    assert node_data['nodes']['old/1968']['subclade'] == 'unassigned'


@needs_iqtree
def test_root_lineage_without_sequences_roots_on_a_child_founder(synthetic):
    # without the old sequence, the root lineage holds only the founders of A and B
    from Bio import SeqIO
    records = [r for r in SeqIO.parse(synthetic / 'alignment.fasta', 'fasta') if r.id != 'old/1968']
    SeqIO.write(records, synthetic / 'alignment.fasta', 'fasta')
    out = synthetic / 'tree.nwk'
    result = subprocess.run([
        sys.executable, str(SCRIPTS / 'tidy_tree.py'),
        '--alignment', str(synthetic / 'alignment.fasta'),
        '--guide-tree', str(synthetic / 'guide.nwk'),
        '--lineage-assignments', str(synthetic / 'lineages.tsv'),
        '--seq-id-column', 'strain', '--lineage-column', 'clade',
        '--metadata', str(synthetic / 'metadata.tsv'),
        '--iqtree-args=-n 2 -czb -seed 1', '--output-tree', str(out),
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    # A holds the oldest sequences of the synthetic data
    assert "no sequences older than those of its sublineages, rooting on the founder of A" in result.stdout
    t = Phylo.read(out, 'newick')
    assert len(t.get_terminals()) == 20
    # the two halves of the root are A and B, with lengths on both
    sides = [{tip.name.split('/')[0][0] for tip in c.get_terminals()} for c in t.root.clades]
    assert sorted(map(sorted, sides)) == [['A'], ['B']]
    assert all(c.branch_length for c in t.root.clades)


@needs_iqtree
def test_young_sequence_of_the_root_lineage_does_not_root_the_tree(synthetic):
    # the only sequence without a lineage call is younger than those of the lineages,
    # as an anomalous sequence would be: the tree is rooted on a founder instead
    meta = (synthetic / 'metadata.tsv').read_text().replace('old/1968\t1968-06-01', 'old/1968\t2030-06-01')
    (synthetic / 'metadata.tsv').write_text(meta)
    result = subprocess.run([
        sys.executable, str(SCRIPTS / 'tidy_tree.py'),
        '--alignment', str(synthetic / 'alignment.fasta'),
        '--guide-tree', str(synthetic / 'guide.nwk'),
        '--lineage-assignments', str(synthetic / 'lineages.tsv'),
        '--seq-id-column', 'strain', '--lineage-column', 'clade',
        '--metadata', str(synthetic / 'metadata.tsv'),
        '--iqtree-args=-n 2 -czb -seed 1', '--output-tree', str(synthetic / 'tree.nwk'),
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "rooting on the founder of A, the lineage of the oldest sequence A/seq0" in result.stdout
