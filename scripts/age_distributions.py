import argparse, sys, os, glob
import numpy as np
from datetime import datetime, timedelta, date
from collections import defaultdict
from Bio import SeqIO, AlignIO
from treetime.utils import numeric_date
from augur.utils import read_metadata, get_numerical_dates
from select_strains import read_strain_list, determine_time_interval, parse_metadata
from flu_regions import region_names, region_properties

def age_distribution(metadata, fname, title=None):
    import matplotlib
    # important to use a non-interactive backend, otherwise will crash on cluster
    matplotlib.use('agg')
    import matplotlib.pyplot as plt

    fs=16
    bins = np.arange(0,100,10)
    bc = 0.5*(bins[1:]+bins[:-1])
    plt.figure()

    for region in region_names:
        props = region_properties[region]
        y,x = np.histogram([m['age'] for m in metadata
                            if m['age']!='unknown' and (m['region']==region or region=='global')], bins=bins)
        total = np.sum(y)
        y = np.array(y, dtype=float)/total
        plt.plot(bc, y, label=props.get('label', region),
                 lw=4 if region=='global' else 2, c=props['color'])

    plt.legend(fontsize=fs*0.8, ncol=2)
    if title:
        plt.title(title, fontsize=1.5*fs)
    plt.ylabel('fraction in age bin', fontsize=fs)
    plt.xlabel('age', fontsize=fs)
    plt.ylim([0,0.7])
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()

def get_title(fname):
    if 'h3n2' in fname:
        return "A(H3N2)"
    elif 'h1n1pdm' in fname:
        return "A(H1N1pdm)"
    elif 'vic' in fname:
        return "B(Vic)"
    elif 'yam' in fname:
        return "B(Yam)"
    else:
        return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Separate strains be region and align specific genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True,
                        help="file with metadata associated with viral sequences")
    parser.add_argument('--output', type=str,  help="names of file to save age_distribution histogram to ")
    parser.add_argument('-r', '--resolution',default='3y', type = str,
                        help = "single resolution to include (default: 3y)")
    parser.add_argument('--time-interval', nargs=2, help="explicit time interval to use -- overrides resolutions"
                                                        " expects YYYY-MM-DD YYYY-MM-DD")
    parser.add_argument('--region', type=str, help="region to draw sequences from")
    parser.add_argument('--exclude', help="a text file containing strains (one per line) that will be excluded")

    args = parser.parse_args()

    region=args.region
    time_interval = sorted([numeric_date(x)
            for x in determine_time_interval(args.time_interval, args.resolution)])
    # read strains to exclude
    excluded_strains = read_strain_list(args.exclude) if args.exclude else []

    # read in meta data, parse numeric dates, and exclude outlier strains
    metadata = {k:val for k,val in parse_metadata(['segment'], [args.metadata]).items()
                if k not in excluded_strains}['segment']

    sequences = []
    print(time_interval)
    for seq in metadata:
        if not "num_date" in metadata[seq]:
            continue
        if metadata[seq]["num_date"]>=time_interval[0] and \
           metadata[seq]["num_date"]<time_interval[1]:
            sequences.append(metadata[seq])

    age_distribution(sequences, args.output, title = get_title(args.output))
