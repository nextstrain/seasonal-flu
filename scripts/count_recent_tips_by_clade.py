"""
Count recent tips by clade for reporting.
"""
import argparse
from augur.utils import read_node_data
from collections import Counter
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Annotate derived haplotypes per node from annotated clades and store as node data JSON",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--recency", help="recency annotations in node data JSON format based on submission date of sequences", required=True)
    parser.add_argument("--clades", help="clade annotations in node data JSON format", required=True)
    parser.add_argument("--recency-values", nargs="+", default=["last week", "last month"], help="values in the 'recency' annotation to use for counting recent sequences")
    parser.add_argument("--output", help="Markdown file with list of recent sequence counts per clade", required=True)
    args = parser.parse_args()

    # Load recency annotations.
    recency = read_node_data(args.recency)

    # Load clades.
    clades = read_node_data(args.clades)

    # Find all nodes with recency values that match the requested values.
    recent_tips = {
        node_name
        for node_name, node_data in recency["nodes"].items()
        if node_data.get("recency") in args.recency_values
    }

    # Count the clade membership for each recent tip.
    count_by_clade = pd.DataFrame.from_dict(
        Counter([
            node_data["subclade"]
            for node_name, node_data in clades["nodes"].items()
            if node_name in recent_tips
        ]),
        orient="index",
        columns=["count"],
    )

    count_by_clade.index.name = "clade"

    # Annotate the total number of tips.
    count_by_clade.loc["total", "count"] = count_by_clade["count"].sum()
    count_by_clade["count"] = count_by_clade["count"].astype(int)

    # Sort by count in descending order.
    count_by_clade = count_by_clade.sort_values(
        "count",
        ascending=False,
    )

    # Save Markdown table.
    with open(args.output, "w", encoding="utf-8") as oh:
        print(count_by_clade.to_markdown(), file=oh)
