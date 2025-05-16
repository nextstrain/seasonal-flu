"""
Convert a frequencies JSON into a single table of values from the last available timepoint labelled by tip or internal node status in a tree.
"""
import argparse
from augur.utils import json_to_tree
import Bio.Phylo
import json
import pandas as pd


def annotate_clades_for_tree(tree):
    # Track all clade memberships in a new attribute to properly handle nested
    # clades.
    for node in tree.find_clades():
        node.node_attrs["clades"] = set([node.node_attrs["emerging_haplotype"]["value"]])
        if node.parent is not None:
            node.node_attrs["clades"].update(node.parent.node_attrs["clades"])

    return tree


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert frequencies JSON to a data frame",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="auspice JSON file for the tree used to estimate the given frequencies")
    parser.add_argument("--frequencies", required=True, help="tip frequencies JSON")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the output table in the format of 'key=value' pairs")
    parser.add_argument("--root-clade", help="name of the clade to select as the root of the given tree such that frequencies are all relative to that root and normalized to sum to 1")
    parser.add_argument("--output", required=True, help="tab-delimited file with frequency per node per timepoint")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data associated with internal nodes in the output table")
    args = parser.parse_args()

    # Load tree.
    with open(args.tree, "r") as fh:
        tree_json = json.load(fh)

    tree = json_to_tree(tree_json)
    tree = annotate_clades_for_tree(tree)

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)

    pivots = frequencies_json.pop("pivots")
    projection_pivot = frequencies_json.pop("projection_pivot")
    frequencies = {
        node_name: region_frequencies["frequencies"]
        for node_name, region_frequencies in frequencies_json.items()
        if node_name not in ["counts", "generated_by"]
    }

    # Zoom into the root clade of the tree.
    if args.root_clade is not None:
        for node in tree.find_clades(terminal=False):
            if node.node_attrs["emerging_haplotype"]["value"] == args.root_clade:
                tree = node
                print(f"Zoom into {args.root_clade} at node {node.name}")
                break

    # Extract frequencies per node for the selected root clade.
    records = []
    for node in tree.find_clades():
        if args.include_internal_nodes or node.is_terminal():
            for i, pivot in enumerate(pivots):
                for clade in node.node_attrs["clades"]:
                    records.append({
                        "strain": node.name,
                        "pivot": pivot,
                        "frequency": float(frequencies[node.name][i]),
                        "is_terminal": node.is_terminal(),
                        "clade": clade,
                        "observed": pivot <= projection_pivot
                    })

    # Convert frequencies data into a data frame.
    df = pd.DataFrame(records)

    # Normalize frequencies to sum to 1 per timepoint.
    if args.root_clade is not None:
        strains_and_pivots = df.loc[:, ["pivot", "strain", "frequency"]].drop_duplicates()
        pivot_sum = strains_and_pivots.groupby("pivot")["frequency"].sum().reset_index()
        df = df.merge(pivot_sum, on="pivot", how="left", suffixes=["", "_sum"])
        df["frequency"] = df["frequency"] / df["frequency_sum"]

    # Add any additional annotations requested by the user in the format of
    # "key=value" pairs where each key becomes a new column with the given
    # value.
    if args.annotations:
        for annotation in args.annotations:
            key, value = annotation.split("=")
            df[key] = value

    # Save the table.
    df.to_csv(args.output, sep="\t", float_format="%.6f", index=False, header=True)
