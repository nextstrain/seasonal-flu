# this code export the sequence json needed for the old deprecated auspice
import argparse, json
from random import sample
import numpy as np
from Bio import Phylo, AlignIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True,
                        help="newick file with the tree")
    parser.add_argument('--alignment', type=str, help="json file with ancestral reconstructions. assumes full sequence")
    parser.add_argument('--translations', nargs='+', help="fasta files with ancestral translations")
    parser.add_argument('--genes', nargs='+', help="names of the genes corresponding to the translations")
    parser.add_argument('--output', type=str, help="names of files to write selected strains to, one for each gene")

    args = parser.parse_args()

    with open(args.alignment) as fh:
        nuc = json.load(fh)["nodes"]

    T = Phylo.read(args.tree, 'newick')
    root_seq=nuc[T.root.name]['sequence']
    sequence_json = {'root':{'nuc':root_seq}}
    for n in T.find_clades(order='preorder'):
        sequence_json[n.name]={'nuc':{p:d for p,a,d in zip(range(len(root_seq)), root_seq, nuc[n.name]['sequence']) if a!=d}}

    for gene, fname in zip(args.genes, args.translations):
        aln = {s.name:str(s.seq) for s in AlignIO.read(fname, 'fasta')}
        root_seq = aln[T.root.name]
        sequence_json['root'][gene]=root_seq
        for n in T.find_clades(order='preorder'):
            sequence_json[n.name][gene] = {p:d for p,a,d in zip(range(len(root_seq)), root_seq, aln[n.name]) if a!=d}

    with open(args.output, 'wt') as fh:
        json.dump(sequence_json, fh)

