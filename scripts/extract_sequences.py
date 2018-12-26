"""Extract sequences from a given FASTA file that match the given list of sample names.
"""
import argparse
import Bio
import Bio.SeqIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract sample sequences by name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file of all sample sequences")
    parser.add_argument("--samples", required=True, help="text file of samples names with one name per line")
    parser.add_argument("--output", required=True, help="FASTA file of extracted sample sequences")
    args = parser.parse_args()

    with open(args.samples) as infile:
        samples = set([line.strip() for line in infile])

    with open(args.output, 'w') as outfile:
        for seq in Bio.SeqIO.parse(args.sequences, 'fasta'):
            if seq.name in samples:
                Bio.SeqIO.write(seq, outfile, 'fasta')
