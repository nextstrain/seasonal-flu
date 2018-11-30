import argparse, sys, os, glob, json
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter

region_label = {'NA': 'N America', 'AS': 'Asia', 'EU': 'Europe', 'OC': 'Oceania', 'global': 'Global'}
cols = [ '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
fs=12
ymax=300
years = YearLocator()
months = MonthLocator(range(1, 13), bymonthday=1, interval=2)
yearsFmt = DateFormatter('%Y')
monthsFmt = DateFormatter("%b")


def load_frequencies(fname):
    with open(fname) as fh:
        return json.load(fh)


def plot_mutations_by_region(frequencies, mutations, fname,show_errorbars=True, n_std_dev=2, n_smooth=3):
    regions = sorted(frequencies.keys())
    smoothed_count_by_region = {}
    for region, freqs in frequencies.items():
        for n, f in freqs.items():
            if 'count' in n:
                gene = n.split(':')[0]
                smoothed_count_by_region[(gene, region)] = np.convolve(np.ones(n_smooth, dtype=float)/n_smooth, f, mode='same')

    fig, axs = plt.subplots(len(mutations), 1, sharex=True, figsize=(8,3+3*len(mutations)))

    for mi,(mut, ax) in enumerate(zip(mutations, axs)):
        gene = mut.split(':')[0]
        for ri,region in enumerate(regions):
            col = 'C%d'%(ri+1)
            pivots = frequencies[region]["pivots"]
            if mut in frequencies[region]:
                tmp_freq = np.array(frequencies[region][mut])
                ax.plot(pivots, tmp_freq, label=region, lw=2, c=col)
                std_dev = np.sqrt(tmp_freq*(1-tmp_freq)/(smoothed_count_by_region[(gene, region)]+1))
                if show_errorbars:
                    ax.fill_between(pivots, tmp_freq-n_std_dev*std_dev, tmp_freq+n_std_dev*std_dev, facecolor=col, linewidth=0, alpha=0.1)
            else:
                print("Mutation %s not calculated in region %s"%(mut, region))
            if mi==0:
                ax.legend()
            if mi==len(mutations)-1:
                ax.set_xlabel('time')
            ax.set_ylabel(mut)
            ax.set_ylim(0,1)

    plt.savefig(fname)


def sample_count_by_region(frequencies, fname):
    regions = sorted(frequencies.keys())
    counts = {}
    for region in frequencies:
        date_bins = frequencies[region]["pivots"]
        tmp = [frequencies[region][x] for x in frequencies[region] if 'counts' in x]
        if len(tmp):
            counts[region] = tmp[0]

    if 'global' not in counts:
        if len(counts)>1:
            counts["global"] = np.sum([x for x in counts.values()], axis=0)
        else:
            counts["global"] = list(counts.values())[0]
    plot_counts(counts, date_bins, fname, drop=3)


def tree_sample_counts(tree_frequencies, fname):
    date_bins = tree_frequencies["pivots"]
    counts = {x.split(':')[0]:tree_frequencies[x] for x in tree_frequencies
              if 'counts' in x}
    if 'global' not in counts:
        counts["global"] = np.sum([x for x in counts.values()], axis=1)

    plot_counts(counts, date_bins, fname, drop=3)


def plot_counts(counts, date_bins, fname, drop=3):
    fig, ax = plt.subplots(figsize=(8, 3))
    regions = sorted(counts.keys())
    tmpcounts = np.zeros(len(date_bins)-drop)
    plt.bar(date_bins[drop:], counts['global'][drop:], width=1.0/13, linewidth=0, label="Other", color="#bbbbbb", clip_on=False)
    for c,region in zip(cols, regions):
        if region!='global':
            plt.bar(date_bins[drop:], counts[region][drop:], bottom=tmpcounts, width=1.0/13, linewidth=0,
                    label=region, color=c, clip_on=False)
            tmpcounts += np.array(counts[region][drop:])
    ax.set_xlim([date_bins[drop-1], date_bins[-1]])
    ax.set_ylim(0,min(max(counts['global']), ymax))
    ax.tick_params(axis='x', which='major', labelsize=fs, pad=20)
    ax.tick_params(axis='x', which='minor', pad=7)
    ax.set_ylabel('Sample count', fontsize=fs*1.1)
    ax.legend(loc=3, ncol=1, bbox_to_anchor=(1.02, 0.53))
    plt.subplots_adjust(left=0.1, right=0.82, top=0.94, bottom=0.22)
    if fname:
        plt.savefig(fname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Separate strains by region and align specific genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--mutation-frequencies', nargs='+', required=True,
                        help="json files containing frequencies in different regions")
    parser.add_argument('--tree-frequencies', type=str, required=True,
                        help="json files containing frequencies of clades in the tree")
    parser.add_argument('--mutations',nargs='+', required=True, help="mutations to graph")
    parser.add_argument('--regions',nargs='+', required=True, help="regions corresponding to alignment files")
    parser.add_argument('--output-mutations', help="file name to save figure to")
    parser.add_argument('--output-total-counts', help="file name to save figure to")
    parser.add_argument('--output-tree-counts', help="file name to save figure to")

    args=parser.parse_args()

    if args.mutation_frequencies:
        frequencies = {}
        for region, fname in zip(args.regions, args.mutation_frequencies):
            frequencies[region] = load_frequencies(fname)
        plot_mutations_by_region(frequencies, args.mutations, args.output_mutations)
        sample_count_by_region(frequencies, args.output_total_counts)


    if args.tree_frequencies:
        tree_frequencies = load_frequencies(args.tree_frequencies)
        tree_sample_counts(tree_frequencies, args.output_tree_counts)
