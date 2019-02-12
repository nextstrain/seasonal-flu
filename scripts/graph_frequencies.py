'''
Plot count distributions as well as global frequencies
'''

import argparse
import json
from datetime import datetime
import numpy as np
import matplotlib
# important to use a non-interactive backend, otherwise will crash on cluster
# this needs to be right after the matplotlib import!
matplotlib.use('agg')

import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import seaborn as sns

default_regions = ['north_america', 'china', 'japan_korea', 'oceania', 'europe', 'southeast_asia']

region_label = {
    'global': 'Global',
    'africa': 'Africa',
    'china': 'China',
    'europe': 'Europe',
    'japan_korea': 'Japan/Korea',
    'north_america': 'N America',
    'oceania': 'Oceania',
    'south_america': 'S America',
    'south_asia': 'S Asia',
    'southeast_asia': 'SE Asia',
    'west_asia': 'W Asia'
}

region_colors = {
    'global': '#111111',
    'africa': '#A0CCA5',
    'china': '#A76BB1',
    'europe': '#658447',
    'japan_korea': '#2A4786',
    'north_america': '#D6C568',
    'oceania': '#8E1616',
    'south_america': '#926224',
    'south_asia': '#EBA85F',
    'southeast_asia': '#8FBDD0',
    'west_asia': '#76104B'
}

fs = 12
ymax = 800
years = YearLocator()
months = MonthLocator(range(1, 13), bymonthday=1, interval=2)
yearsFmt = DateFormatter('%Y')
monthsFmt = DateFormatter("%b")
sns.set_style('ticks')

def load_frequencies(fname):
    with open(fname) as fh:
        return json.load(fh)

def plot_mutations_by_region(frequencies, mutations, fname, show_errorbars=True,
                             n_std_dev=1, n_smooth=3, drop=3, regions=None):

    if regions is None:
        regions = default_regions

    smoothed_count_by_region = {}
    # generate a temporally smoothed sample count vector for each region.
    for region, freqs in frequencies.items():
        if region=='global':
            continue
        for n, f in freqs.items():
            if 'count' in n:
                gene = n.split(':')[0]
                smoothed_count_by_region[(gene, region)] = np.convolve(np.ones(n_smooth, dtype=float)/n_smooth, f, mode='same')

    # set up a figure and plot each mutation in a different panel
    rows = int(np.ceil(len(mutations)/2.))
    fig, axs = plt.subplots(rows, 2, sharex=True, figsize=(10, 1.0+1.6*rows))

    if len(mutations) % 2 != 0:
        axs[-1, -1].axis('off')

    axs = list(np.array(axs).flatten())

    region_labeled = set()
    for mi,(mut, ax) in enumerate(zip(mutations, axs)):
        gene = mut.split(':')[0]
        for region in regions:
            pivots = frequencies[region]["pivots"]
            offset = datetime(2000,1,1).toordinal()
            pivots = [offset+(x-2000)*365.25 for x in pivots]
            if mut in frequencies[region]:
                tmp_freq = np.array(frequencies[region][mut])
                ax.plot(pivots[:-drop], tmp_freq[:-drop], '-o',
                        ms=7 if region=='global' else 4, lw=3 if region=='global' else 1,
                        label=region_label.get(region, region) if region not in region_labeled else '',
                        c=region_colors[region])
                region_labeled.add(region)
                if show_errorbars and region!="global":
                    std_dev = np.sqrt(tmp_freq*(1-tmp_freq)/(smoothed_count_by_region[(gene, region)]+1))
                    ax.fill_between(pivots[:-drop], (tmp_freq-n_std_dev*std_dev)[:-drop],
                                    (tmp_freq+n_std_dev*std_dev)[:-drop],
                                    facecolor=region_colors[region], linewidth=0, alpha=0.1)
            else:
                print("Mutation %s not calculated in region %s"%(mut, region))
            # if mi==0:
            #     ax.legend(ncol=1, bbox_to_anchor=(1.02, 0.2))
            ax.set_ylabel(mut, fontsize=fs)
            ax.set_ylim(0, 1)
            ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]])
            ax.tick_params(axis='x', which='major', labelsize=fs, pad=20)
            ax.tick_params(axis='x', which='minor', pad=7)
            ax.xaxis.set_major_locator(years)
            ax.xaxis.set_major_formatter(yearsFmt)
            ax.xaxis.set_minor_locator(months)
            ax.xaxis.set_minor_formatter(monthsFmt)

    fig.legend(loc=1, ncol=2)
    plt.tight_layout()
    sns.despine()
    plt.savefig(fname)


def plot_clades_by_region(frequencies, clades, clade_to_node, fname, show_errorbars=True,
                          regions=None, n_std_dev=1, n_smooth=3, drop=3):

    if regions is None:
        regions = default_regions

    smoothed_count_by_region = {}
    total_count_by_region = {}
    for region in frequencies['counts']:
        smoothed_count_by_region[region] = np.convolve(np.ones(n_smooth, dtype=float)/n_smooth,
                                                       frequencies['counts'][region], mode='same')
        total_count_by_region[region] = np.sum(frequencies['counts'][region])

    rows = int(np.ceil(len(clades)/2.))
    fig, axs = plt.subplots(rows, 2, sharex=True, figsize=(10, 1.0+1.6*rows))

    if len(clades) % 2 != 0:
        axs[-1, -1].axis('off')

    axs = list(np.array(axs).flatten())

    for mi,(clade, ax) in enumerate(zip(clades, axs)):
        if clade not in clade_to_node:
            print("Clade %s is not annotated"%clade)
            continue

        node = clade_to_node[clade]
        for region in regions:
            pivots = frequencies["pivots"]
            offset = datetime(2000,1,1).toordinal()
            pivots = [offset+(x-2000)*365.25 for x in pivots]
            if node in frequencies and region in frequencies[node]:
                tmp_freq = np.array(frequencies[node][region])
                ax.plot(pivots[:-drop], tmp_freq[:-drop], '-o',
                        ms=7 if region=='global' else 4, lw=3 if region=='global' else 1,
                        label=region_label.get(region, region), c=region_colors[region])
                std_dev = np.sqrt(tmp_freq*(1-tmp_freq)/(smoothed_count_by_region[region]+1))
                if show_errorbars:
                    ax.fill_between(pivots[:-drop], (tmp_freq-n_std_dev*std_dev)[:-drop],
                                    (tmp_freq+n_std_dev*std_dev)[:-drop],
                                    facecolor=region_colors[region], linewidth=0, alpha=0.1)
            else:
                print("region %s not present in node %s"%(region, node))
            ax.set_ylabel(clade, fontsize=fs)
            ax.set_ylim(0, 1)
            ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]])
            ax.tick_params(axis='x', which='major', labelsize=fs, pad=20)
            ax.tick_params(axis='x', which='minor', pad=7)
            ax.xaxis.set_major_locator(years)
            ax.xaxis.set_major_formatter(yearsFmt)
            ax.xaxis.set_minor_locator(months)
            ax.xaxis.set_minor_formatter(monthsFmt)

    fig.legend(regions, loc=1, ncol=2)
    plt.tight_layout()
    sns.despine()
    plt.savefig(fname)


def sample_count_by_region(frequencies, fname, regions=None):

    if regions is None:
        regions = default_regions

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
    plot_counts(counts, date_bins, fname, drop=3, regions=regions)


def tree_sample_counts(tree_frequencies, fname, regions=None):

    if regions is None:
        regions = default_regions

    date_bins = tree_frequencies["pivots"]
    counts = tree_frequencies["counts"]
    if 'global' not in counts:
        counts["global"] = np.sum([x for x in counts.values()], axis=1)

    plot_counts(counts, date_bins, fname, drop=3, regions=regions)


def plot_counts(counts, date_bins, fname, drop=3, regions=None):

    if regions is None:
        regions = default_regions

    offset = datetime(2000,1,1).toordinal()
    date_bins = [offset+(x-2000)*365.25 for x in date_bins]

    fig, ax = plt.subplots(figsize=(8, 3))
    tmpcounts = np.zeros(len(date_bins))
    width = 0.75*(date_bins[1] - date_bins[0])
    plt.bar(date_bins, counts['global'], width=width, linewidth=0,
            label="Other", color="#bbbbbb", clip_on=False)

    for region in regions:
        if region!='global':
            plt.bar(date_bins, counts[region], bottom=tmpcounts, width=width, linewidth=0,
                    label=region_label.get(region, region), color=region_colors[region],
                    clip_on=False, alpha=0.8)
            tmpcounts += np.array(counts[region])
    ax.set_xlim([date_bins[0]-width*0.5, date_bins[-1]])
    ax.set_ylim(0,min(max(counts['global']), ymax))
    ax.tick_params(axis='x', which='major', labelsize=fs, pad=20)
    ax.tick_params(axis='x', which='minor', pad=7)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.xaxis.set_minor_locator(months)
    ax.xaxis.set_minor_formatter(monthsFmt)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('Sample count', fontsize=fs*1.1)
    ax.legend(loc=3, ncol=1, bbox_to_anchor=(1.01, 0.35))
    plt.subplots_adjust(left=0.08, right=0.81, top=0.9, bottom=0.22)
    if fname:
        plt.savefig(fname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Separate strains by region and align specific genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--mutation-frequencies', type=str,
                        help="json files containing frequencies in different regions")
    parser.add_argument('--tree-frequencies', type=str,
                        help="json files containing frequencies of clades in the tree")
    parser.add_argument('--clade-annotations', type=str,
                        help="json files containing clade annotations to map internal nodes to clades")
    parser.add_argument('--mutations',nargs='+', help="mutations to graph")
    parser.add_argument('--clades',nargs='+', help="clades to graph")
    parser.add_argument('--regions',nargs='+', required=True, help="regions corresponding to alignment files")
    parser.add_argument('--output-mutations', help="file name to save figure to")
    parser.add_argument('--output-total-counts', help="file name to save figure to")
    parser.add_argument('--output-tree-counts', help="file name to save figure to")
    parser.add_argument('--output-clades', help="file name to save figure to")

    args=parser.parse_args()

    if args.mutation_frequencies:
        frequencies = load_frequencies(args.mutation_frequencies)
        plot_mutations_by_region(frequencies, args.mutations, args.output_mutations, drop=1, regions=args.regions)
        sample_count_by_region(frequencies, args.output_total_counts, regions=args.regions)


    if args.tree_frequencies:
        tree_frequencies = load_frequencies(args.tree_frequencies)
        tree_sample_counts(tree_frequencies, args.output_tree_counts, regions=args.regions)

        if args.clade_annotations:
            clade_annotations = load_frequencies(args.clade_annotations)
            clade_to_node = {node["clade_annotation"]:node_name for node_name, node in clade_annotations['nodes'].items()
                             if "clade_annotation" in node}

            plot_clades_by_region(tree_frequencies, args.clades, clade_to_node,
                                  args.output_clades, regions=args.regions, drop=1)
