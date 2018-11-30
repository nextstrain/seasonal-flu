import numpy as np
from .graph_frequencies import load_frequencies

def rising_mutations(freq_file, dn=5, offset=0, baseline = 0.01, fname='tmp.txt', n_out=-1):
    '''
    safe a file containing all mutations and summary of their recent frequency trajectories.
    mutations are sorted by their log derivative over the past dn month
    '''
    frequencies = load_frequencies(freq_file)
    dx = {}
    counts_by_gene = {}
    for label, f in frequencies.items():
        if 'counts' in label:
            gene, mut = label.split(':')
            counts_by_gene[gene] = np.array(f)

    for label, f in frequencies.items():
        if 'counts' in label or 'pivots' in label:
            continue
        gene, mut = label.split(':')
        ind = np.arange(len(f))[-(dn+offset):len(f)-offset]
        tmp_freq = np.array(f)
        c = np.sum(counts_by_gene[gene][ind]*tmp_freq[ind])
        tmp_x = np.mean(tmp_freq[ind])
        tmp_dx = tmp_freq[-1-offset] - tmp_freq[-dn-offset]
        dx[(gene, mut)] = (tmp_x, tmp_dx, tmp_dx/(tmp_x+baseline), c)

    with open(fname, 'w') as ofile:
        ofile.write("#Frequency change over the last %d month\n"%dn)
        for k,v in sorted(dx.items(), key=lambda x:x[1][2], reverse=True)[:n_out]:
            ofile.write("%s:%s\t%1.3f\t%1.3f\t%1.3f\t%1.1f\n"%(k[0], k[1], v[0], v[1], v[2], v[3]))
        ofile.write('\n')


def recurring_mutations(aa_mutations, fname_by_position='tmp.txt', fname_by_mutation='tmp.txt', n_out=-1):
    '''
    count the number of times that each position has mutated on the tree and save to file.
    in additition, count the number of mutations that resulted in a specific substitution
    '''
    from collections import defaultdict
    by_mutation = defaultdict(int)
    by_position = defaultdict(int)
    aa_muts = load_frequencies(aa_mutations)

    for name, n in aa_muts["nodes"].items():
        for gene in n["aa_muts"]:
            for mut in n["aa_muts"][gene]:
                a, pos, d = mut[0], mut[1:-1], mut[-1]
                by_position[(gene, pos)]+=1
                by_mutation[(gene, mut[1:])]+=1

    with open(fname_by_position, 'w') as ofile:
        ofile.write("#Number of independent mutation events in the tree observed at a specific position\n")
        for mut,val in sorted(by_position.items(), key=lambda x:x[1], reverse=True)[:n_out]:
            ofile.write("%s:%s\t%d\n"%(mut[0], mut[1], val))

    with open(fname_by_mutation, 'w') as ofile:
        ofile.write("#Number of sepcific independent mutation events in the tree observed\n")
        for mut,val in sorted(by_mutation.items(), key=lambda x:x[1], reverse=True)[:n_out]:
            ofile.write("%s:%s\t%d\n"%(mut[0], mut[1], val))

