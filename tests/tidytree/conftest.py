import io
import random
import sys
from pathlib import Path

import pytest
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SCRIPTS = Path(__file__).resolve().parents[2] / 'scripts'
sys.path.insert(0, str(SCRIPTS))

from tidytree.guide import GuideTree  # noqa: E402

# A.1.1 is a sublineage of A.1; B has no sublineages
GUIDE_NEWICK = "(((A.1.1)A.1,A.2)A,B)unassigned;"


def tree(newick):
    return Phylo.read(io.StringIO(newick), 'newick')


@pytest.fixture
def guide(tmp_path):
    path = tmp_path / 'guide.nwk'
    path.write_text(GUIDE_NEWICK)
    return GuideTree(str(path))


def mutate(seq, n, rng):
    seq = list(seq)
    for pos in rng.sample(range(len(seq)), n):
        seq[pos] = rng.choice([b for b in 'ACGT' if b != seq[pos]])
    return ''.join(seq)


@pytest.fixture
def synthetic(tmp_path):
    """
    Founders evolving along the guide tree and a few sequences per lineage, written as
    an alignment with the founders, a lineage table, metadata and the guide tree.
    """
    rng = random.Random(1)
    length = 600
    root = ''.join(rng.choice('ACGT') for _ in range(length))
    founders = {'A': mutate(root, 12, rng), 'B': mutate(root, 12, rng)}
    founders['A.1'] = mutate(founders['A'], 8, rng)
    founders['A.2'] = mutate(founders['A'], 8, rng)
    founders['A.1.1'] = mutate(founders['A.1'], 6, rng)

    records, rows, meta = [], [], []
    year = 2000
    for lineage, founder in founders.items():
        for i in range(4):
            name = f"{lineage}/seq{i}"
            records.append(SeqRecord(Seq(mutate(founder, rng.randint(1, 3), rng)), id=name, description=''))
            rows.append((name, lineage))
            meta.append((name, f"{year + i}-01-15"))
        year += 5
    # an old sequence without a lineage call, at the root of the nomenclature
    records.append(SeqRecord(Seq(mutate(root, 2, rng)), id='old/1968', description=''))
    rows.append(('old/1968', ''))
    meta.append(('old/1968', '1968-06-01'))
    records.extend(SeqRecord(Seq(seq), id=name, description='') for name, seq in founders.items())

    from Bio import SeqIO
    SeqIO.write(records, tmp_path / 'alignment.fasta', 'fasta')
    (tmp_path / 'lineages.tsv').write_text(
        'strain\tclade\n' + ''.join(f"{name}\t{lineage}\n" for name, lineage in rows))
    (tmp_path / 'metadata.tsv').write_text(
        'strain\tdate\n' + ''.join(f"{name}\t{date}\n" for name, date in meta))
    (tmp_path / 'guide.nwk').write_text(GUIDE_NEWICK)
    return tmp_path
