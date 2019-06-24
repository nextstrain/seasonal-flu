"""Extract sequences from a given FASTA file that match the given list of sample names.
"""
import numpy as np
import argparse
from Bio import AlignIO, SeqIO, Seq, SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
from augur.translate import safe_translate
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract sample sequences by name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of aligned sequences")
    parser.add_argument("--reference", required=True, help="annotated genbank file")
    parser.add_argument("--output", required=True, help="FASTA file of extracted sample sequences")
    args = parser.parse_args()

    aln = AlignIO.read(args.alignment, 'fasta')
    ref = SeqIO.read(args.reference, 'genbank')

    cds_start, cds_end = np.inf, 0

    for feature in ref.features:
        if feature.type=='CDS':
            if feature.location.start<cds_start:
                cds_start=feature.location.start
            if feature.location.end>cds_end:
                cds_end=feature.location.end

    UTR5p = {s.id:s for s in aln[:,:cds_start]} if cds_start else None
    UTR3p = {s.id:s for s in aln[:,cds_end:]} if cds_end<aln.get_alignment_length() else None

    cds = aln[:,cds_start:cds_end]
    ungapped_aa = []
    ungapped = {}
    for seq in cds:
        str_seq = str(seq.seq)
        left_gaps = len(str_seq) - len(str_seq.lstrip('-'))
        right_gaps = len(str_seq) - len(str_seq.rstrip('-'))
        ungapped[seq.id] = '-'*left_gaps + str(seq.seq.ungap('-')) + '-'*right_gaps
        aa = safe_translate(ungapped[seq.id])
        if aa.count('X') + aa.count('*')<5:
            ungapped_aa.append(SeqRecord.SeqRecord(seq=Seq.Seq(aa), id=seq.id, name=seq.id, description=''))
        else:
            print(seq.id, "didn't translate properly")

    tmp_outfile = args.output+'.tmp.fasta'
    SeqIO.write(ungapped_aa, tmp_outfile, 'fasta')

    tmp_aln_file = args.output+'.tmp_aln.fasta'
    os.system("mafft --auto %s > %s"%(tmp_outfile, tmp_aln_file))


    aa_aln = {seq.id: seq for seq in AlignIO.read(tmp_aln_file, 'fasta')}

    new_cds_aln = []
    for seq in cds:
        pos=0
        nuc_seq_aln = [str(UTR5p[seq.id].seq) if UTR5p else '']
        nuc_seq = ungapped[seq.id]
        if seq.id in aa_aln:
            for aa in aa_aln[seq.id].seq:
                if aa=='-':
                    nuc_seq_aln.append('---')
                    if nuc_seq[pos:pos+3]=='---':
                        pos+=3
                else:
                    if len(nuc_seq)>=pos+3:
                        nuc_seq_aln.append(nuc_seq[pos:pos+3])
                    else:
                        nuc_seq_aln.append('---')
                    pos+=3

            seq.seq=Seq.Seq(''.join(nuc_seq_aln)+(str(UTR3p[seq.id].seq) if UTR3p else ''))
            new_cds_aln.append(seq)
        else:
            print(seq.id, "didn't translate properly")

    AlignIO.write(MultipleSeqAlignment(new_cds_aln), args.output, 'fasta')
    os.remove(tmp_outfile)
    os.remove(tmp_aln_file)
