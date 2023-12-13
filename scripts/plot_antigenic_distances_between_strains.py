#!/usr/bin/env python3
import matplotlib
# important to use a non-interactive backend, otherwise will crash on cluster
matplotlib.use('agg')

import argparse

from augur.utils import read_colors, read_strains
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--antigenic-distances", required=True, help="antigenic distances between strains")
    parser.add_argument("--auspice-config", required=True, help="Auspice config JSON with a custom color scale for clades.")
    parser.add_argument("--auspice-config-color-field", default="subclade", help="key name for the coloring in the Auspice config JSON to find colors for clades")
    parser.add_argument("--clades", help="a list of clades for which test strains should be plotted")
    parser.add_argument("--clade-color-field", default="clade_test", help="field in antigenic distance table to color measurements by")
    parser.add_argument("--references", help="a list of reference strains to plot in the order they should be displayed from top to bottom")
    parser.add_argument("--references-to-include", help="a list of reference strains to force include in plots")
    parser.add_argument("--references-to-exclude", help="a list of reference strains to exclude from plots")
    parser.add_argument("--top-references-to-keep", type=int, default=10, help="top N number of references to keep by number of titer measurements")
    parser.add_argument("--min-test-date", type=float, help="minimum numeric date for test strains to include in plots")
    parser.add_argument("--title", help="title to use for the figure")
    parser.add_argument("--plot-raw-data", action="store_true", help="plot raw data as points in the scatterplot in addition to the summary statistics of mean and confidence intervals")
    parser.add_argument("--output", required=True, help="plot of distances between strains")

    args = parser.parse_args()

    # Configure plotting style.
    sns.set_style("ticks")
    mpl.rcParams['figure.dpi'] = 300
    mpl.rcParams['font.size'] = 16
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['legend.fontsize'] = 14
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16

    # Load distances.
    df = pd.read_csv(
        args.antigenic_distances,
        sep="\t"
    )

    # Filter by minimum test strain date.
    if args.min_test_date:
        df = df[df["test_date"] >= args.min_test_date].copy()

    # Filter by clades.
    if args.clades:
        with open(args.clades, "r", encoding="utf-8") as fh:
            clades = [clade.strip() for clade in fh]

        df = df[df[args.clade_color_field].isin(clades)].copy()

    if df.shape[0] == 0:
        print("Could not find any measurements after filtering by clades and min date", file=sys.stderr)
        sys.exit(1)

    # Annotate reference names with clade.
    df["reference_name"] = df.apply(
        lambda row: f"{row['reference_strain']}\n({row['subclade_reference']})",
        axis=1
    )

    # Load colors.
    with open(args.auspice_config, "r", encoding="utf-8") as fh:
        auspice_config = json.load(fh)

    clade_colors = None
    for coloring in auspice_config["colorings"]:
        if coloring["key"] == args.auspice_config_color_field:
            clade_colors = dict(coloring["scale"])
            break

    if clade_colors is None:
        print("Could not load clade colors from the Auspice config JSON.", file=sys.stderr)
        sys.exit(1)

    reference_order = None
    if args.references:
        # If references are given in order, we take these as they are.
        with open(args.references, "r") as fh:
            references_to_use = [line.strip() for line in fh]

        record_by_strain = df.loc[
            :,
            ["reference_strain", "reference_name"]
        ].drop_duplicates().set_index("reference_strain").to_dict(orient="index")

        reference_order = [
            record_by_strain[strain]["reference_name"]
            for strain in references_to_use
            if strain in record_by_strain
        ]
    else:
        # Otherwise, look for references to populate the plot.
        # Load references to include or exclude.
        references_to_include = set()
        if args.references_to_include:
            references_to_include = read_strains(args.references_to_include)

            # Limit references to include to those present in the data.
            references_to_include &= set(df["reference_strain"].values)

        references_to_exclude = set()
        if args.references_to_exclude:
            references_to_exclude = read_strains(args.references_to_exclude)

        # Take the references with the most titer measurements, in addition to the
        # references requested by the user.
        top_references_to_keep = max(args.top_references_to_keep - len(references_to_include), 0)
        references_to_use = set(df["reference_strain"].value_counts().head(top_references_to_keep).index.values)

        # Force include specific references.
        references_to_use.update(references_to_include)

        # Exclude specific references.
        references_to_use -= references_to_exclude

    # Filter distances to those with requested references.
    filtered_df = df[df["reference_strain"].isin(references_to_use)].copy()

    # Use the user-defined clade order, when possible. If clades haven't been
    # provided, order test clades by global frequency in descending order.
    if args.clades:
        clade_order = clades
    else:
        clade_order = filtered_df.sort_values(
            "clade_frequency_test",
            ascending=False
        )[args.clade_color_field].drop_duplicates().values

    if reference_order is None:
        # Order reference strains by clade (in descending order of global
        # frequency) and mean log2 distance to test strains in ascending order.
        filtered_df["negative_clade_frequency_reference"] = -1 * filtered_df["clade_frequency_reference"]

        reference_order = filtered_df.groupby([
            "reference_name",
            "negative_clade_frequency_reference"
        ])["log2_titer"].mean().reset_index().sort_values([
            "negative_clade_frequency_reference",
            "log2_titer",
        ])["reference_name"].values

    # Assign colors to clades.
    number_of_clades = len(clade_order)
    color_by_clade = {
        clade: clade_colors.get(clade, "#CCCCCC")
        for clade in clade_order
    }

    # Initialize the figure
    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    sns.despine()

    if args.plot_raw_data:
        # Show each observation with a scatterplot
        sns.stripplot(
            x="log2_titer",
            y="reference_name",
            hue=args.clade_color_field,
            order=reference_order,
            hue_order=clade_order,
            data=filtered_df,
            palette=color_by_clade,
            dodge=True,
            size=8,
            alpha=0.4,
            jitter=0.2,
            zorder=1
        )

    # Plot conditional means.
    sns.pointplot(
        x="log2_titer",
        y="reference_name",
        hue=args.clade_color_field,
        order=reference_order,
        hue_order=clade_order,
        data=filtered_df,
        dodge=0.55,
        join=False,
        palette=color_by_clade,
        markers="d",
        scale=0.75,
        errorbar=("ci", 89),
    )

    # Draw a line at the origin to show where we expect effectively inhibited
    # viruses to appear.
    ax.axvline(
        x=0.0,
        color="#000000",
        alpha=0.75,
        zorder=-10,
    )

    # Draw a line at the traditional threshold used to denote antigenic drift.
    ax.axvline(
        x=2.0,
        color="#000000",
        alpha=0.25,
        zorder=-10,
        linestyle="dashed",
    )

    ax.set_xlabel("$\log_{2}$ normalized titer")
    ax.set_ylabel("Reference strain")

    # Improve the legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[:number_of_clades],
        labels[:number_of_clades],
        title="clade (test strains)",
        handletextpad=0,
        columnspacing=1,
        ncol=1,
        frameon=False
    )

    # Add zebra striping to rows to more clearly delineate which data belong to
    # which reference. Implementation based on:
    # https://stackoverflow.com/questions/2815455/matplotlib-zebra-stripe-a-figures-background-color
    yticks,_ = plt.yticks()
    xticks = plt.xticks()[0]
    left = xticks[0]
    right = xticks[-1]
    width = right - left
    height = yticks[1] - yticks[0]

    # create bars at yticks that are the length of our greatest xtick and have a height equal to our tick spacing
    ax.barh(
        yticks,
        [width] * len(yticks),
        height=height,
        left=left,
        color=["#FFFFFF", "#CCCCCC"],
        zorder=-20,
        alpha=0.25,
    )

    if args.title:
        plt.title(args.title)

    plt.tight_layout()
    plt.savefig(args.output)
