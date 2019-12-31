import argparse, sys, os, glob
import numpy as np
from datetime import datetime, timedelta, date
from random import sample
from collections import defaultdict
from Bio import SeqIO, AlignIO, SeqRecord, Seq
from treetime.utils import numeric_date
from augur.utils import read_metadata, get_numerical_dates, load_features
from augur import align
from augur.translate import safe_translate
from select_strains import read_strain_list, determine_time_interval, parse_metadata
from codon_align import codon_align, get_cds

class pseudo_args(object):
    def __init__(self, **kwargs):
        for k,v in kwargs.items():
            self.__setattr__(k,v)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Separate strains be region and align specific genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True,
                        help="file with metadata associated with viral sequences")
    parser.add_argument('--sequences', type=str, required=True,
                        help="FASTA file of virus sequences from fauna (e.g., zika.fasta)")
    parser.add_argument('--output', nargs='+', help="names of files to write selected strains to, one for each gene")
    parser.add_argument('-l', '--lineage', choices=['h3n2', 'h1n1pdm', 'vic', 'yam'],
                        default='h3n2', type=str,
                        help="single lineage to include (default: h3n2)")
    parser.add_argument('-r', '--resolution',default='3y', type = str,
                        help = "single resolution to include (default: 3y)")
    parser.add_argument('--time-interval', nargs=2, help="explicit time interval to use -- overrides resolutions"
                                                        " expects YYYY-MM-DD YYYY-MM-DD")
    parser.add_argument('--reference-sequence', required=True,
                        help='GenBank file containing the annotation')
    parser.add_argument('--exclude', help="a text file containing strains (one per line) that will be excluded")
    parser.add_argument('--genes', nargs='+', help="genes to be translated")
    parser.add_argument('--region', type=str, help="region to draw sequences from")

    args = parser.parse_args()

    region=args.region
    time_interval = sorted([numeric_date(x)
            for x in determine_time_interval(args.time_interval, args.resolution)])
    # read strains to exclude
    excluded_strains = read_strain_list(args.exclude) if args.exclude else []

    # read in meta data, parse numeric dates, and exclude outlier strains
    metadata = {k:val for k,val in parse_metadata(['segment'], [args.metadata]).items()
                if k not in excluded_strains}['segment']

    sequences = []
    for seq in SeqIO.parse(args.sequences, 'fasta'):
        if seq.name in metadata:
            if metadata[seq.name]["num_date"]>=time_interval[0] and \
               metadata[seq.name]["num_date"]<time_interval[1] and \
               metadata[seq.name]["region"]==region:
                sequences.append(seq)

    tmp_str = "".join(sample('ABCDEFGHILKLMOPQRSTUVWXYZ', 20))

    ref = SeqIO.read(args.reference_sequence, 'genbank')

    # get sequence as string, CDS seq, amino acid sequence, and start/end pos
    features_to_translate = load_features(args.reference_sequence, args.genes)
    refstr, refCDS, refAA, cds_start, cds_end = get_cds(ref)

    alignment = []
    for seq in sequences:
        seq_aln = codon_align(seq,  refstr, refAA, cds_start, cds_end)
        if seq_aln:
            seq.seq=Seq.Seq(seq_aln)
            alignment.append(seq)

    print("selected %d for region %s and date interval %f-%f"%(len(sequences), region, time_interval[0], time_interval[1]))

    for gene, fname in zip(args.genes, args.output):
        if gene not in features_to_translate:
            continue
        seqs = []
        feature = features_to_translate[gene]
        for seq in alignment:
            try:
                translation =  SeqRecord.SeqRecord(seq=Seq.Seq(safe_translate(str(feature.extract(seq.seq)))),
                                                   id=seq.name, name=seq.name, description='')
                seqs.append(translation)
            except:
                print("WARN:",seq.name,"did not translate")
        SeqIO.write(seqs, fname, 'fasta')
