"""
Count recent tips by clade for reporting.
"""
import argparse
from augur.utils import read_node_data
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Count the number of recent sequences per clade",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--recency", help="recency annotations in node data JSON format based on submission date of all sequences", required=True)
    parser.add_argument("--clades", help="subclade annotations from a Nextclade TSV", required=True)
    parser.add_argument("--recency-values", nargs="+", default=["last week", "last month"], help="values in the 'recency' annotation to use for counting recent sequences")
    parser.add_argument("--output", help="Markdown file with list of recent sequence counts per clade", required=True)
    args = parser.parse_args()

    # Load recency annotations.
    recency = read_node_data(args.recency)

    # Find all nodes with recency values that match the requested values.
    recent_tips = {
        node_name
        for node_name, node_data in recency["nodes"].items()
        if node_data.get("recency") in args.recency_values
    }

    # Load clade labels per sequence.
    clades = pd.read_csv(
        args.clades,
        sep="\t",
        usecols=["seqName", "subclade", "qc.overallStatus"],
    )


    # Filter clade labels to recent non-low-quality sequences and count the
    # clade membership for each recent tip.
    count_by_clade = clades[
        (clades["qc.overallStatus"] != "bad") &
        (clades["seqName"].isin(recent_tips))
    ].groupby(
        "subclade"
    )["seqName"].count().reset_index(
        name="count"
    ).rename(
        columns={"subclade": "clade"},
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
