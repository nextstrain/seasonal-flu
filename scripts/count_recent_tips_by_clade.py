"""
Count recent tips by clade for reporting.
"""
import argparse
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Count the number of recent sequences per clade",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--recency", help="recency annotations in TSV format based on submission date of all sequences", required=True)
    parser.add_argument("--clades", help="clade annotations from a Nextclade TSV", required=True)
    parser.add_argument("--clade-column", default="subclade", help="name of the column to use for clade annotations in the given Nextclade TSV")
    parser.add_argument("--output", help="Markdown file with list of recent sequence counts per clade", required=True)
    args = parser.parse_args()

    # Load recency annotations.
    recency = pd.read_csv(
        args.recency,
        sep="\t",
    )
    recent_tips = set(recency["strain"].values)

    # Load clade labels per sequence.
    clades = pd.read_csv(
        args.clades,
        sep="\t",
        usecols=["strain", args.clade_column],
    )


    # Filter clade labels to recent non-low-quality sequences and count the
    # clade membership for each recent tip.
    count_by_clade = clades[
        (clades["strain"].isin(recent_tips))
    ].groupby(
        args.clade_column
    )["strain"].count().reset_index(
        name="count"
    ).rename(
        columns={args.clade_column: "clade"},
    )
    count_by_clade["count"] = count_by_clade["count"].astype(int)

    # Sort by count in descending order.
    count_by_clade = count_by_clade.sort_values(
        "count",
        ascending=False,
    )

    # Annotate the total number of tips.
    total_count = pd.DataFrame.from_dict({
        "clade": ["total"],
        "count": [count_by_clade["count"].sum()],
    })
    count_by_clade = pd.concat([
        count_by_clade,
        total_count,
    ])

    # Save Markdown table.
    with open(args.output, "w", encoding="utf-8") as oh:
        print(count_by_clade.to_markdown(index=False), file=oh)
