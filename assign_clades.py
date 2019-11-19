#!/usr/bin/env python3
"""Extract sequences from a given FASTA file that match the given list of sample names.
"""
import numpy as np
import argparse, sys
from Bio import AlignIO, SeqIO, Seq, SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
from augur.translate import safe_translate
from augur.clades import read_in_clade_definitions, is_node_in_clade
from augur.utils import load_features
from scripts.codon_align import align_pairwise, get_cds, codon_align

class tmpNode(object):
    def __init__(self):
        self.sequences = {}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign clades to sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file of HA sequences")
    parser.add_argument("--lineage", required=True, help="lineage of the sequences supplied")
    args = parser.parse_args()

    refname = f"config/reference_{args.lineage}_ha.gb"
    seqs = SeqIO.parse(args.sequences, 'fasta')
    ref = SeqIO.read(refname, 'genbank')
    features = load_features(refname)
    clade_designations = read_in_clade_definitions(f"config/clades_{args.lineage}_ha.tsv")

    # get sequence as string, CDS seq, amino acid sequence, and start/end pos
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(ref)

    alignment = []
    for seq in seqs:
        seq_container = tmpNode()
        seq_aln = codon_align(seq,  refstr, refAA, cds_start, cds_end)
        if seq_aln is None:
            print(f"{seq.id}\tnot translatable", file=sys.stdout)
            continue

        seq_container.sequences['nuc'] = {i:c for i,c in enumerate(seq_aln)}
        for fname, feat in features.items():
            if feat.type != 'source':
                seq_container.sequences[fname] = {i:c for i,c in enumerate(safe_translate(feat.extract(seq_aln)))}

        matches = []
        for clade_name, clade_alleles in clade_designations.items():
            if is_node_in_clade(clade_alleles, seq_container, ref):
                matches.append(clade_name)
        print(f"{seq.description}\t{', '.join(matches)}", file=sys.stdout)
