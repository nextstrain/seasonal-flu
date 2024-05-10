#!/usr/bin/env python3
import argparse
from matplotlib import pyplot as plt
import pandas as pd


def load_lineage_dates(metadata_file, lineage):
    df = pd.read_csv(metadata_file, sep="\t").dropna(subset=["date"])
    df["lineage"] = lineage
    dates = df.loc[~(df["date"].str.contains("X")), ["lineage", "date", "region"]]

    return dates


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", nargs="+", help="one metadata file per lineage")
    parser.add_argument("--lineages", nargs="+", help="one lineage name per metadata file")
    parser.add_argument("--min-date", help="minimum date to plot counts by")
    parser.add_argument("--output-count-by-lineage", help="plot of metadata record counts per lineage and month")
    parser.add_argument("--output-h1n1pdm-count", help="plot of H1N1pdm counts per region and month")
    parser.add_argument("--output-h3n2-count", help="plot of H3N2 counts per region and month")
    parser.add_argument("--output-vic-count", help="plot of Vic counts per region and month")

    args = parser.parse_args()

    dates = pd.concat([
        load_lineage_dates(metadata, lineage)
        for metadata, lineage in zip(args.metadata, args.lineages)
    ]).query("(date != '36-09-05') & (date != '?')").copy()

    recent_dates = dates[dates["date"] > args.min_date].copy()
    recent_dates["date"] = pd.to_datetime(recent_dates["date"])

    regions = [
        region
        for region in sorted(recent_dates["region"].drop_duplicates().values)
        if region != "?"
    ]

    binned_counts = recent_dates.set_index(
        "date"
    ).groupby(
        "lineage"
    ).resample(
        "1MS"
    ).count().rename(
        columns={"lineage": "samples"}
    ).reset_index()

    # Plot counts per lineage per month.
    fig, ax = plt.subplots(1, 1, figsize=(10, 3), dpi=300)
    for lineage, lineage_counts in binned_counts.groupby("lineage"):
        ax.plot(
            lineage_counts["date"],
            lineage_counts["samples"],
            "o-",
            label=lineage,
            alpha=0.75,
        )

    ax.set_xlabel("Date")
    ax.set_ylabel("Number of samples in GISAID")
    ax.legend(
        frameon=False,
        loc='upper right',
        borderaxespad=0.0,
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(args.output_count_by_lineage)
    plt.close()

    # Bin counts by region and time per lineage.
    output_names = [
        args.output_h1n1pdm_count,
        args.output_h3n2_count,
        args.output_vic_count,
    ]
    for lineage, lineage_counts in binned_counts.groupby("lineage"):
        binned_region_counts = recent_dates.query(
            f"(lineage == '{lineage}') & (region != '?')"
        ).set_index(
            "date"
        ).groupby(
            "region"
        ).resample(
            "1MS"
        ).count().rename(
            columns={"region": "samples"}
        ).reset_index()

        fig, ax = plt.subplots(1, 1, figsize=(10, 3), dpi=300)
        for region, region_counts in binned_region_counts.groupby("region"):
            ax.plot(
                region_counts["date"],
                region_counts["samples"],
                "o-",
                label=region,
                alpha=0.75,
            )

        ax.set_xlabel("Date")
        ax.set_ylabel("Number of samples in GISAID")
        ax.legend(
            frameon=False,
            bbox_to_anchor=(1.05, 1),
            loc='upper left',
            borderaxespad=0.0,
        )

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        output_name = output_names.pop(0)
        plt.savefig(output_name)
        plt.close()
