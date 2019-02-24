import argparse, sys, os, glob
import numpy as np
from datetime import datetime, timedelta, date
from random import sample
from collections import defaultdict
from Bio import SeqIO, AlignIO
from treetime.utils import numeric_date
from augur.utils import read_metadata, get_numerical_dates, load_features
from augur import align
from select_strains import read_strain_list, regions, determine_time_interval, parse_metadata


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
    if not os.path.isdir('tmp'):
        os.mkdir('tmp')

    print("selected %d for region %s and date interval %f-%f"%(len(sequences), region, time_interval[0], time_interval[1]))
    features_to_translate = load_features(args.reference_sequence, args.genes)
    tmp_file = "tmp/sequence_file_%s_%s.fasta"%(region, tmp_str)
    tmp_file_out = "tmp/sequence_file_%s_%s_aln.fasta"%(region, tmp_str)
    SeqIO.write(sequences, tmp_file, 'fasta')
    fail = align.run(pseudo_args(sequences=tmp_file, reference_sequence=args.reference_sequence,
                      output = tmp_file_out, reference_name=None, remove_reference=True,
                      method='mafft', nthreads=1, fill_gaps=False))
    if fail:
        sys.exit(fail)

    aln = AlignIO.read(tmp_file_out, 'fasta')
    for seq in aln:
        seq_str = str(seq.seq)
        seq_str = seq_str.lstrip('-')
        seq_str = 'X'*(len(seq)-len(seq_str)) + seq_str

        seq_str = seq_str.rstrip('-')
        seq_str = seq_str + 'X'*(len(seq)-len(seq_str))

    for gene, fname in zip(args.genes, args.output):
        if gene not in features_to_translate:
            continue
        seqs = []
        feature = features_to_translate[gene]
        for seq in aln:
            try:
                translation = feature.extract(seq).translate(gap='-')
                translation.id, translation.name, translation.description = seq.name, seq.name, ''
                seqs.append(translation)
            except:
                print("WARN:",seq.name,"did not translate")
        SeqIO.write(seqs, fname, 'fasta')

    for fname in glob.glob("tmp/*%s_%s*"%(region, tmp_str)):
        os.remove(fname)
