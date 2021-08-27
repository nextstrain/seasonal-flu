#!/usr/bin/env python3
import matplotlib
# important to use a non-interactive backend, otherwise will crash on cluster
matplotlib.use('agg')

import argparse

from augur.utils import read_strains
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

NAME_BY_LINEAGE = {
    "h3n2": "A/H3N2",
    "h1n1pdm": "A/H1N1pdm",
    "vic": "B/Victoria",
}

NAME_BY_SOURCE = {
    "cdc": "CDC",
    "crick": "Crick",
    "niid": "NIID",
    "vidrl": "VIDRL",
}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--antigenic-distances", required=True, help="antigenic distances between strains")
    parser.add_argument("--color-schemes", required=True, help="table of color schemes for different number of items")
    parser.add_argument("--references-to-include", help="a list of reference strains to force include in plots")
    parser.add_argument("--references-to-exclude", help="a list of reference strains to exclude from plots")
    parser.add_argument("--top-references-to-keep", type=int, default=10, help="top N number of references to keep by number of titer measurements")
    parser.add_argument("--min-test-date", type=float, help="minimum numeric date for test strains to include in plots")
    parser.add_argument("--min-clade-frequency", type=float, default=0.05, help="minimum frequency for clades whose test viruses should be plotted")
    parser.add_argument("--plot-raw-data", action="store_true", help="plot raw data as points in the scatterplot in addition to the summary statistics of mean and confidence intervals")
    parser.add_argument("--output", required=True, help="plot of distances between strains")

    args = parser.parse_args()

    # Configure plotting style.
    sns.set_style("ticks")
    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14

    # Load distances.
    df = pd.read_csv(
        args.antigenic_distances,
        sep="\t"
    )

    # Filter by minimum test strain date.
    if args.min_test_date:
        df = df[df["test_date"] >= args.min_test_date].copy()

    # Load color schemes.
    color_schemes = []
    with open(args.color_schemes, "r") as fh:
        for line in fh:
            color_schemes.append(line.strip().split("\t"))

    # Load references to include or exclude.
    references_to_include = set()
    if args.references_to_include:
        references_to_include = read_strains(args.references_to_include)

        # Limit references to include to those present in the data.
        references_to_include &= set(df["reference_strain"].values)

    references_to_exclude = set()
    if args.references_to_exclude:
        references_to_exclude = read_strains(args.references_to_exclude)

    # Filter distances to test viruses from clades at the requested minimum
    # frequency and reference viruses from clades that still exist.
    df = df[
        (df["clade_frequency_test"] > args.min_clade_frequency) & (df["clade_frequency_reference"] > 0.0)
    ].copy()

    # Take the references with the most titer measurements, in addition to the
    # references requested by the user.
    top_references_to_keep = max(args.top_references_to_keep - len(references_to_include), 0)
    recent_references = set(df["reference_strain"].value_counts().head(top_references_to_keep).index.values)

    # Force include specific references.
    recent_references.update(references_to_include)

    # Exclude specific references.
    recent_references -= references_to_exclude

    # Filter distances to those with requested references.
    filtered_df = df[df["reference_strain"].isin(recent_references)].copy()

    # Annotate reference names with clade.
    filtered_df["reference_name"] = filtered_df.apply(
        lambda row: f"{row['reference_strain']}\n({row['clade_reference']})",
        axis=1
    )

    # Order test clades by global frequency in descending order.
    clade_order = filtered_df.sort_values(
        "clade_frequency_test",
        ascending=False
    )["clade_test"].drop_duplicates().values

    # Order reference strains by clade (in descending order of global frequency)
    # and mean log2 distance to test strains in ascending order.
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
    clade_colors = color_schemes[number_of_clades - 1]
    color_by_clade = dict(zip(reversed(clade_order), clade_colors))

    # Get features of current dataset for display in the title.
    lineage = NAME_BY_LINEAGE[filtered_df["lineage"].values[0]]
    passage = filtered_df["passage"].values[0]
    assay = filtered_df["assay"].values[0].upper()
    sources = ", ".join([NAME_BY_SOURCE[source] for source in filtered_df["source"].drop_duplicates().sort_values().values])

    # Initialize the figure
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    sns.despine()

    if args.plot_raw_data:
        # Show each observation with a scatterplot
        sns.stripplot(
            x="log2_titer",
            y="reference_name",
            hue="clade_test",
            order=reference_order,
            hue_order=clade_order,
            data=filtered_df,
            palette=color_by_clade,
            dodge=True,
            size=8,
            alpha=0.5,
            jitter=0.2,
            zorder=1
        )

    # Plot conditional means.
    sns.pointplot(
        x="log2_titer",
        y="reference_name",
        hue="clade_test",
        order=reference_order,
        hue_order=clade_order,
        data=filtered_df,
        dodge=0.55,
        join=False,
        alpha=0.5,
        palette=color_by_clade,
        markers="d",
        scale=0.75,
        ci=89,
        legend=False
    )

    # Draw a line at the traditional threshold used to denote antigenic drift.
    ax.axvline(
        x=2.0,
        color="#000000",
        alpha=0.25,
        zorder=-10
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

    plt.title(f"{lineage} {passage} {assay} titers ({sources})")

    plt.tight_layout()
    plt.savefig(args.output)
