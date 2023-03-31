'''
Plot count distributions as well as global frequencies
'''

import argparse
import pandas as pd
import datetime
import numpy as np
from flu_regions import region_properties, region_names
import matplotlib.pyplot as plt

def parse_date(date_str):
    from augur.utils import ambiguous_date_to_date_range
    try:
        date_range = ambiguous_date_to_date_range(date_str, fmt="%Y-%m-%d")
        return date_range[0]
    except:
        return None

def plot_totals(data, date_bins):
    regions = [x for x in sorted(data.loc[~pd.isna(data['region']), 'region'].unique()) if x!='?']
    time_points  = [datetime.datetime(year=x[0], month=x[1], day=15) for x in date_bins]
    fig = plt.figure()
    width= datetime.timedelta(days=25)
    tmpcounts = np.zeros(len(date_bins))
    for i,region in enumerate(regions):
        by_region = data.loc[data.region==region, "year-month"].value_counts()
        counts = np.array([by_region.get(date,0) for date in date_bins])
        plt.bar(time_points, counts,
                bottom=tmpcounts, width=width, linewidth=0,
                label=region, color=f'C{i}',
                clip_on=False, alpha=0.8)
        tmpcounts += np.array(counts)
    plt.legend()
    plt.ylabel("sample count")
    fig.autofmt_xdate()
    return fig

def clade_frequencies_by_region(data, date_bins):
    regions = [x for x in sorted(data.loc[~pd.isna(data['region']), 'region'].unique()) if x!='?']
    time_points  = [datetime.datetime(year=x[0], month=x[1], day=15) for x in date_bins]
    fig, axs = plt.subplots(len(regions)//2, 2, figsize=(6,12))
    clades = sorted([x for x in data.clade.unique() if not pd.isna(x)])

    for ai,region in enumerate(regions):
        ax = axs[ai//2,ai%2]
        ax.set_title(region)
        totals = data.loc[data.region==region, ["year-month"]].value_counts()
        total_counts_by_month = np.array([float(totals.get(date,1)) for date in date_bins])
        print(total_counts_by_month)
        for ci, clade in enumerate(clades):
            clade_counts_tmp = data.loc[np.logical_and(data.region==region, data.clade==clade), ["year-month"]].value_counts()
            clade_counts_by_month = np.array([float(clade_counts_tmp.get(date,0)) for date in date_bins])
            freqs = clade_counts_by_month/total_counts_by_month
            if freqs.max()>0.05:
                freqs[total_counts_by_month<5] = np.nan
                ax.plot(time_points, freqs, label=clade, c=f"C{ci%10}")
            ax.set_ylim(0,1)

    plt.legend(ncol=3)
    plt.ylabel("frequency")
    fig.autofmt_xdate()
    return fig


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Separate strains by region and align specific genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str,
                        help="tsv files containing metadata")
    parser.add_argument('--nextclade', type=str,
                        help="tsv files containing nextclade output")
    parser.add_argument('--output-total-counts', help="file name to save figure to")
    parser.add_argument('--output-by-clade', help="file name to save figure to")
    parser.add_argument('--output-by-region', help="file name to save figure to")

    args=parser.parse_args()

    meta = pd.read_csv(args.metadata, index_col=0, sep='\t')
    meta["cdate"] = meta.date.apply(parse_date)
    nextclade = pd.read_csv(args.nextclade, index_col=0, sep='\t')["clade"]

    d = pd.concat([meta, nextclade],axis=1)

    d['year-month'] = d.cdate.apply(lambda x:(x.year, x.month) if x else None)

    start_date = (2020,1)
    today = (datetime.date.today().year, datetime.date.today().month)
    date_bins = []
    for year in range(start_date[0], today[0]+1):
        for month in range(start_date[1] if year==start_date[0] else 1, today[1]+1 if year==today[0] else 13):
            date_bins.append((year, month))

    if args.output_total_counts:
        fig = plot_totals(d, date_bins)
        fig.savefig(args.output_total_counts)

    if args.output_by_region:
        fig = clade_frequencies_by_region(d, date_bins)
        fig.savefig(args.output_by_region)
