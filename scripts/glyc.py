# this code export the sequence json needed for the old deprecated auspice
import argparse, json
from random import sample
import numpy as np
from Bio import Phylo, AlignIO
import re

def glycosylation_count(total_aa_seq, glyc_mask=None):
    if glyc_mask is None:
        glyc_mask = np.ones(len(total_aa_seq), dtype=bool)

    # TODO: need to restrict to surface residues.
    total_aa_seq_masked = "".join([aa if mask else 'X'
                                   for (mask, aa) in zip(glyc_mask, total_aa_seq)])

    return len(re.findall('N[^P][ST][^P]', total_aa_seq_masked))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True,
                        help="newick file with the tree")
    parser.add_argument('--alignment',  help="fasta file with ancestral translations")
    parser.add_argument('--output', type=str, help="names of files to write selected strains to, one for each gene")

    args = parser.parse_args()

    T = Phylo.read(args.tree, 'newick')

    glyc_json = {}
    aln = {s.name:str(s.seq) for s in AlignIO.read(args.alignment, 'fasta')}
    root_seq = aln[T.root.name]
    root_glyc = glycosylation_count(root_seq)
    for n in T.find_clades(order='preorder'):
        glyc_json[n.name] = {'glyc':glycosylation_count(aln[n.name]) - root_glyc}

    with open(args.output, 'wt') as fh:
        json.dump({'nodes':glyc_json, 'comment':"glycosylation motif count in HA1 relative to root sequence."}, fh)

