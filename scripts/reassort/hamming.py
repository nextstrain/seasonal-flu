import argparse, sys, os, glob, json
import collections
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(
        description="Calculate hamming distance between strains",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--fasta', required=True, help="fasta file")
    parser.add_argument('--output', required=True, help="output json file containing strain-strain hamming distances")
    args=parser.parse_args()
    return args


def hamming_distance(seqObj1, seqObj2):
    atcg = ["A", "a", "T", "t", "C", "c", "G", "g"]
    assert len(seqObj1.seq) == len(seqObj2.seq)
    return sum(x != y for x, y in zip(seqObj1.seq, seqObj2.seq) if x in atcg and y in atcg)


if __name__ == '__main__':
    """
    Compute the hamming distance between all strains in the fasta provided.
    Currently writes out [x][y] and [y][x] despite them being equal.
    """
    args = get_args()
    records = list(SeqIO.parse(args.fasta, "fasta"))
    d = collections.defaultdict(dict)
    for a in records:
        for b in records:
            d[a.name][b.name] = hamming_distance(a, b)

    # double check symmetrical
    for a in records:
        for b in records:
            assert(d[a.name][b.name] == d[b.name][a.name])

    with open(args.output, 'w') as fh:
        json.dump(d, fh, indent=2, sort_keys = True)
