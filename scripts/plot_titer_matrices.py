import argparse, sys, os, glob, json
from collections import defaultdict
import numpy as np
import matplotlib
import pandas as pd
# important to use a non-interactive backend, otherwise will crash on cluster
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

def load_json(fname):
    with open(fname) as fh:
        return json.load(fh)

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
    parser.add_argument('--clades', help="json with clade info")
    parser.add_argument('--clades-to-plot', nargs='+', help="clades to include in matrix")
    parser.add_argument('--antigens', nargs="+", help="antigens to summarize titers for")
    parser.add_argument('--combine-sera', action='store_true', help="average values for different sera")

    parser.add_argument('--output', help="file name to save figure to")

    args=parser.parse_args()

    titers = load_json(args.titers)
    clades = load_json(args.clades)["nodes"]

    autologous_titers = get_autologous_titers(titers)
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

    df = pd.DataFrame(average_titers).T
    df.sort_index(inplace=True)
    try:
        df.pop('unassigned')
    except:
        pass

    plt.figure()
    sns.heatmap(df, vmin=-1, vmax=5, cmap='YlOrRd')
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
