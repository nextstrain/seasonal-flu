from treetime import TreeTime
from treetime.utils import parse_dates
import argparse
import numpy as np

def calc_node_timings(T, sigma_sq, mu, eps=0.2):
    for n in T.find_clades(order='postorder'):
        if not n.keep: continue
        if n.is_terminal():
            prefactor = (n.observations/sigma_sq + mu**2/(n.nmuts+eps))
            n.a = (n.avg_date/sigma_sq + mu*n.nmuts/(n.nmuts+eps))/prefactor
        else:
            if n==T.root:
                tmp_children_1 = mu*np.sum([(mu*c.a-c.nmuts)/(eps+c.nmuts) for c in n if c.keep])
                tmp_children_2 = mu**2*np.sum([(1-c.b)/(eps+c.nmuts) for c in n if c.keep])
                prefactor = (n.observations/sigma_sq + tmp_children_2)
                n.a = (n.observations*n.avg_date/sigma_sq + tmp_children_1)/prefactor
            else:
                tmp_children_1 = mu*np.sum([(mu*c.a-c.nmuts)/(eps+c.nmuts) for c in n if c.keep])
                tmp_children_2 = mu**2*np.sum([(1-c.b)/(eps+c.nmuts) for c in n if c.keep])
                prefactor = (n.observations/sigma_sq + mu**2/(n.nmuts+eps) + tmp_children_2)
                n.a = (n.observations*n.avg_date/sigma_sq + mu*n.nmuts/(n.nmuts+eps)+tmp_children_1)/prefactor
        n.b = mu**2/(n.nmuts+eps)/prefactor
    T.root.tau = T.root.a

    for n in T.get_nonterminals(order='preorder'):
        for c in n:
            if not c.keep: c.tau=None
            else: c.tau = c.a + c.b*n.tau

def calc_scores_to_optimize(x, T):
    mu, sigma = x
    return calc_scores(T, sigma=sigma, mu=mu)['cost']

def calc_scores(T, sigma=None, mu=None):
    sigma_sq=sigma**2
    z_distribution = []
    calc_node_timings(T, sigma_sq=sigma_sq, mu=mu)
    cost = 0
    n_tips = 0
    for n in T.find_clades():
        if not n.keep: continue
        for x in n.tips.values():
            x['z'] = (x['date']-n.tau)/sigma
            z_distribution.append(x['z'])
            cost += x['z']**2
        for c in n:
            if c.keep:
                cost += (mu*(c.tau-n.tau) - c.nmuts)**2/(c.nmuts+1)
        n_tips += n.observations

    res = 0.5*cost + np.log(2*np.pi*(sigma_sq+0.1))*n_tips*0.5
    return {'cost':res, 'z_stddev':np.std(z_distribution)}

def prepare_tree(T):
    for n in T.get_nonterminals(order='postorder'):
        n.dates = []
        n.tips = {}
        for c in n:
            if c.is_terminal():
                if len(c.mutations)==0:
                    c.keep=False
                    if c.raw_date_constraint is not None:
                        n.dates.append(c.raw_date_constraint)
                        n.tips[c.name]={'date':c.raw_date_constraint}
                else:
                    c.keep=True
                    c.dates = [c.raw_date_constraint]
                    c.observations = 1
                    c.avg_date = c.raw_date_constraint
                    c.tips = {c.name:{'date':c.raw_date_constraint}}
                    c.nmuts = len([m for m in c.mutations if m[-1] in 'ACGT'])

        n.keep=any([c.keep for c in n])

        n.nmuts = len(n.mutations)
        n.observations = len(n.dates)
        n.avg_date = np.mean(n.dates) if n.observations else 0

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Flag outliers in a tree')
    parser.add_argument('--tree', type=str, help='tree file in newick format')
    parser.add_argument('--aln', type=str, help='alignment file in fasta format')
    parser.add_argument('--cutoff', type=float, default=4.0, help="z-score used to flag outliers")
    parser.add_argument('--optimize', action="store_true", help="optimize sigma and mu")
    parser.add_argument('--dates', type=str, help='csv/tsv file with dates for each sequence')
    parser.add_argument('--output-outliers', type=str, help='file for outliers')
    parser.add_argument('--output-tree', type=str, help='file for pruned tree')

    args = parser.parse_args()
    dates = parse_dates(args.dates)
    tt = TreeTime(gtr='JC69', tree=args.tree, aln=args.aln, verbose=1, dates=dates)
    tt.clock_filter(n_iqd=4, reroot='least-squares')
    if args.aln:
        tt.infer_ancestral_sequences(prune_short=True, marginal=True)

    pruned_tips = prepare_tree(tt.tree)

    mu = tt.clock_model['slope']*tt.data.full_length
    # magic number: allowing for slack in timing equivalent to 3 mutations
    # this is rescaled later such that the empirical z-score distributions as variance 1
    sigma = 3/mu
    if args.optimize:
        from scipy.optimize import minimize
        x0=(mu, sigma)
        sol = minimize(calc_scores_to_optimize, x0=x0, args=(tt.tree,), method='Nelder-Mead')
        mu = sol['x'][0]
        sigma = sol['x'][1]

    print(f"Calculating node timings using {mu=:1.3e}/year and {sigma=:1.3e}years")
    res = calc_scores(tt.tree, mu=mu, sigma=sigma)


    outliers = []
    for n in tt.tree.find_clades():
        if not n.keep: continue
        for tip, s in n.tips.items():
            if np.abs(s['z'])>args.cutoff*res['z_stddev']:
                outliers.append({"sequence": tip,
                                 "z_score": s['z']/res['z_stddev'],
                                 "expected_date": n.tau,
                                 "given_date": s['date']})

    import pandas as pd
    df = pd.DataFrame(outliers)
    print(df)
    if args.output_outliers:
        df.to_csv(args.output_outliers, index=False, sep='\t')

    if args.output_tree:
        from Bio import Phylo
        T = tt.tree
        for r, row in df.iterrows():
            T.prune(row['sequence'])
        Phylo.write(T, args.output_tree, 'newick')
