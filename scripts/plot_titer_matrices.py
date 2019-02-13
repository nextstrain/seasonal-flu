import argparse, sys, os, glob, json
from collections import defaultdict
import numpy as np
import matplotlib
# important to use a non-interactive backend, otherwise will crash on cluster
matplotlib.use('agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from select_strains import read_strain_list, regions, determine_time_interval, parse_metadata


h3n2_clades = ['A1', 'A1a', 'A1b', 'A1b/135K', 'A1b/135N', 'A1b/131K','A2', 'A2/re', 'A3', 'A4', '3c3.A']
h1n1_clades = ["6b1.A", "6b1.A/183P-1", "6b1.A/183P-2", "6b1.A/183P-3", "6b1.A/183P-5", "6b1.A/183P-6", "6b1.A/183P-7"]
vic_clades =  ["V1A", "V1A.1", "V1A/165N"]
yam_clades = ["172Q", "3"]

def load_json(fname):
    with open(fname) as fh:
        return json.load(fh)

def reduce_to_sera(titers, sera):
    sub_titers = dict()
    for refstrain in titers:
        tmp = defaultdict(dict)
        for testvirus in titers[refstrain]:
            for serum in sera:
                if serum in titers[refstrain][testvirus]:
                    tmp[testvirus][serum] = titers[refstrain][testvirus][serum]
        if len(tmp):
            sub_titers[refstrain]=tmp

    return sub_titers


def get_autologous_titers(titers):
    autologous_titers = defaultdict(dict)
    for refstrain in titers:
        for teststrain in titers[refstrain]:
            if refstrain==teststrain:
                for serum, val in titers[refstrain][refstrain].items():
                    autologous_titers[refstrain][serum] = val
    return autologous_titers

def get_viruses_by_clade(clades):
    viruses_by_clade = defaultdict(list)
    for v,d in clades.items():
        if (not v.startswith('NODE_')) and 'clade_membership' in d:
            viruses_by_clade[d['clade_membership']].append(v)
    return viruses_by_clade


def get_average_titer_by_clade(titers, clades, normalized=False,
                               median=False, geometric=False, serum=None):
    ti = 0 if normalized else 1
    def mean_func(d):
        if median:
            return np.median(d)
        elif geometric:
            return np.exp(np.mean(np.log(d)))
        else:
            return np.mean(d)

    titers_by_clade = defaultdict(list)
    for teststrain in titers:
        if teststrain not in clades:
            continue
        if serum is None:
            titers_by_clade[clades[teststrain]['clade_membership']].append(mean_func([x[ti] for x in titers[teststrain].values()]))
        elif serum in titers[teststrain]:
            titers_by_clade[clades[teststrain]['clade_membership']].append(titers[teststrain][serum][ti])

    return {k:mean_func(d) for k,d in titers_by_clade.items()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Use json summarizing titers to plot average titers by clades",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--titers', help="json with titer info")
    parser.add_argument('--model', help="json with titer model used to subtract serum potencies")
    parser.add_argument('--clades', help="json with clade info")
    parser.add_argument('--metadata', type=str, help="metadata table")
    parser.add_argument('--clades-to-plot', nargs='+', help="clades to include in matrix")
    parser.add_argument('--antigens', nargs="+", help="antigens to summarize titers for")
    parser.add_argument('--combine-sera', action='store_true', help="average values for different sera")
    parser.add_argument('--reassortants', nargs="+", help="list of sera raised against reassortants")

    parser.add_argument('--output', help="file name to save figure to")


    args=parser.parse_args()
    # metadata = {k:val for k,val in parse_metadata(['segment'], [args.metadata]).items()}['segment']

    date_cutoff = 2018

    titers = load_json(args.titers)
    potency = load_json(args.model)["potency"] if args.model else defaultdict(dict)

    if args.reassortants:
        print("reducing to ", args.reassortants)
        titers = reduce_to_sera(titers, args.reassortants)

    autologous_titers = get_autologous_titers(titers)
    # prune old measurements
    to_pop = []
    for antigen in titers:
        titers[antigen] = {k:v for k,v in titers[antigen].items()
                            if any([x in k for x in ['/2018', '/2019']])}
#                           if metadata[k]["num_date"]>date_cutoff}
        if len(titers[antigen])<5 or any([x in antigen for x in ["/2014", "/2015"]]):
            to_pop.append(antigen)

    for a in to_pop:
        print("dropping reference strain ",a)
        titers.pop(a)

    clades = load_json(args.clades)["nodes"]

    viruses_by_clade = get_viruses_by_clade(clades)

    # if no antigens are specified, take the top 15
    if args.antigens:
        antigens = args.antigens
    else:
        antigens_by_clades = defaultdict(list)
        for serum in titers:
            clade = clades[serum]['clade_membership']
            antigens_by_clades[clade].append(serum)

        selected_antigens = set()
        for clade in antigens_by_clades:
            selected_antigens.update(sorted(antigens_by_clades[clade], key=lambda x:len(titers[x]), reverse=True)[:2])

        antigens = selected_antigens.union(sorted(titers.keys(), key=lambda x:len(titers[x]), reverse=True)[:15])

    # summarize titers for each antigen
    average_titers = {}
    for antigen in antigens:
        if antigen in clades:
            c = clades[antigen]['clade_membership']
        else:
            continue
        average_titers[(c, antigen)] = get_average_titer_by_clade(titers[antigen], clades,
                                        normalized=True, geometric=False)
        for clade in average_titers[(c,antigen)]:
            average_titers[(c,antigen)][clade] -= potency[antigen].get("mean_potency",0)



    df = pd.DataFrame(average_titers).T
    # sort columns
    df = df[[x for x in h3n2_clades+h1n1_clades+vic_clades+yam_clades
             if x in df.columns]]
    rows = list(df.iterrows())
    rows.sort(key=lambda x:('Z',x[0][1]) if x[0][0]=='3c3.A' else x[0])
    df = pd.DataFrame({x[0]:x[1] for x in rows}).T

    try:
        df.pop('unassigned')
    except:
        pass

    plt.figure()
    cmap = sns.cubehelix_palette(start=2.6, rot=.1, as_cmap=True)
    sns.heatmap(df, vmin=0, vmax=4, cmap=cmap, square=True, cbar_kws={"shrink": min(1.0,(df.shape[0]+.1)/(df.shape[1]+.1))})
    plt.ylabel('')
    plt.xlabel('')
    tick_labels = []
    for x in df.T:
        median_auto = np.median([y[1] for y in autologous_titers[x[1]].values()])
        if np.isnan(median_auto):
            tick_labels.append("%s - %s  NaN"%(x[1], x[0]))
        else:
            tick_labels.append("%s - %s % 4d"%(x[1], x[0], median_auto))
    plt.yticks(0.5+np.arange(len(df)), tick_labels)
    plt.xticks(1.0+np.arange(len(df.columns)), [x for x in df.columns], rotation=60, horizontalalignment='right')
    plt.tight_layout()

    if args.output:
        plt.savefig(args.output)
