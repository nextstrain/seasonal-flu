"""Extract sequences from a given FASTA file that match the given list of sample names.
"""
import numpy as np
import argparse
from Bio import AlignIO, SeqIO, Seq, SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
from augur.translate import safe_translate
from seqanpy import align_overlap

scoring_params = {"score_match":3, "score_mismatch":-1, "score_gapext":-1, "score_gapopen":-10}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract sample sequences by name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file of aligned sequences")
    parser.add_argument("--reference", required=True, help="annotated genbank file")
    parser.add_argument("--output", required=True, help="FASTA file of extracted sample sequences")
    args = parser.parse_args()

    aln = SeqIO.parse(args.sequences, 'fasta')
    ref = SeqIO.read(args.reference, 'genbank')

    # assuming there is one contiguous coding region which might be
    # split into multiple sub-proteins like HA1 and HA2.
    # loop over all features, pull out min and max of their union
    cds_start, cds_end = np.inf, 0
    for feature in ref.features:
        if feature.type=='CDS':
            if feature.location.start<cds_start:
                cds_start=feature.location.start
            if feature.location.end>cds_end:
                cds_end=feature.location.end

    refstr = str(ref.seq).upper()
    refCDS = refstr[cds_start:cds_end]
    refAA = safe_translate(refstr[cds_start:cds_end])

    alignment = []
    for seq in aln:
        seqstr = str(seq.seq).upper()
        score, refaln, seqaln = align_overlap(refstr, seqstr, **scoring_params)
        if score<0: # did not align
            continue
        ref_aln_array = np.array(list(refaln))
        seq_aln_array = np.array(list(seqaln))

        # stip gaps
        ungapped = ref_aln_array!='-'
        ref_aln_array_ungapped = ref_aln_array[ungapped]
        seq_aln_array_ungapped = seq_aln_array[ungapped]

        seq5pUTR = "".join(seq_aln_array_ungapped[:cds_start])
        seq3pUTR = "".join(seq_aln_array_ungapped[cds_end:])
        seqCDS = "".join(seq_aln_array_ungapped[cds_start:cds_end])
        seqCDS_ungapped = seqCDS.replace('-', '')
        seqAA = safe_translate(seqCDS_ungapped)

        scoreAA, refalnAA, seqalnAA = align_overlap(refAA, seqAA,**scoring_params)
        if scoreAA<0 or sum(seqAA.count(x) for x in ['*', 'X'])>5 or refalnAA.count('-')>5:
            print(seq.id, "didn't translate properly")
            continue

        seqCDS_aln = seq5pUTR
        pos = 0
        for aa_ref, aa_seq in zip(refalnAA, seqalnAA):
            if aa_seq=='-':
                seqCDS_aln += '---'
                # if the nucleotide sequence is gapped
                # (i.e. because of missing data at the 5p and 3p end, advance pos)
                if seqCDS_ungapped[pos:pos+3]=='---':
                    pos += 3
            else:
                if len(seqCDS_ungapped)>=pos+3:
                    seqCDS_aln += seqCDS_ungapped[pos:pos+3]
                else:
                    seqCDS_aln += '---'
                pos += 3
        if len(seqCDS_aln)>len(refstr):
            import ipdb; ipdb.set_trace()
        seq.seq=Seq.Seq(''.join(seqCDS_aln)+seq3pUTR)
        alignment.append(seq)

    # output
    AlignIO.write(MultipleSeqAlignment(alignment), args.output, 'fasta')
