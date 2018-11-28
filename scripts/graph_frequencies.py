import argparse, sys, os, glob, json
import numpy as np
from matplotlib import pyplot as plt

def load_frequencies(fname):
    with open(fname) as fh:
        return json.load(fh)

def plot_mutations_by_region(frequencies, mutations):
    regions = sorted(frequencies.keys())

    fig, axs = plt.subplots(len(mutations), 1, sharex=True, figsize=(5+3*len(mutations), 6))

    for mut, ax in zip(mutations, axs):
        for region in regions:
            #ax.plot(frequencies[region][mut])
            ax.plot(frequencies[region]["pivots"], frequencies[region][mut])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Separate strains be region and align specific genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--frequencies', nargs='+', required=True,
                        help="json files containing frequencies in different regions")
    parser.add_argument('--output', help="file name to save figure to")
    parser.add_argument('--mutations',nargs='+', required=True, help="mutations to graph")

    args=parser.parse_args()

    frequencies = {}
    for fname in args.frequencies:
        region = fname.split('_')[-4]
        frequencies[region] = load_frequencies(fname)


    plot_mutations_by_region(frequencies, args.mutations)
