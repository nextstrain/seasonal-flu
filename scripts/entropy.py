# this code export the entropy json needed for the old deprecated auspice
import argparse, json
from random import sample
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from treetime.seq_utils import alphabets
from augur.utils import load_features

def calc_SNV_frequencies(aln, alphabet='ACGT-'):
    aln_array = np.array(aln)
    f = np.zeros((len(alphabet), aln_array.shape[1]))
    for ni,nuc in enumerate(alphabet):
        f[ni] = np.mean(aln_array==nuc, axis=0)
    f /= np.sum(f,axis=0)

    return f

def calc_entropy(af, aa=False, start=0):
    res = {}
    res['val'] = [round(x,4) for x in np.sum(-af*np.log(af+1e-15), axis=0)]
    res['pos'] = [int(x) for x in start + (3 if aa else 1)*np.arange(af.shape[1])]
    if aa:
        res['codon'] = [int(x) for x in np.arange(af.shape[1])]
    else:
        res['codon'] = [int(x) for x in np.arange(af.shape[1]//3).repeat(3)]
    return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--alignment', type=str, required=True,
                        help="FASTA file of virus sequences from fauna (e.g., zika.fasta)")
    parser.add_argument('--output', type=str, help="names of files to write selected strains to, one for each gene")
    parser.add_argument('--reference-sequence', required=True,
                        help='GenBank file containing the annotation')
    parser.add_argument('--genes', nargs='+', help="genes to be translated")


    args = parser.parse_args()
    alignment = AlignIO.read(args.alignment, 'fasta')
    features_to_translate = load_features(args.reference_sequence, args.genes)

    entropy = {}
    entropy['nuc'] = calc_entropy(calc_SNV_frequencies(alignment, alphabets['nuc']))

    for gene in args.genes:
        if gene not in features_to_translate:
            continue
        seqs = []
        feature = features_to_translate[gene]
        for seq in alignment:
            try:
                translation = feature.extract(seq).translate(gap='-')
                translation.id, translation.name, translation.description = seq.name, seq.name, ''
                seqs.append(translation)
            except:
                print("WARN:",seq.name,"did not translate")
        entropy[gene] = calc_entropy(calc_SNV_frequencies(MultipleSeqAlignment(seqs), alphabets['aa']),
                                    start=feature.location.start, aa=True)

    with open(args.output, 'wt') as fh:
        json.dump(entropy, fh)

