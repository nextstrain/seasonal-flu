"""
Convert a frequencies JSON into a single table of values from the last available timepoint labelled by tip or internal node status in a tree.
"""
import argparse
import Bio.Phylo
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert frequencies JSON to a data frame",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick file for the tree used to estimate the given frequencies")
    parser.add_argument("--frequencies", required=True, help="frequencies JSON")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the output table in the format of 'key=value' pairs")
    parser.add_argument("--output", required=True, help="tab-delimited file with frequency per node at the last available timepoint")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data associated with internal nodes in the output table")
    parser.add_argument("--minimum-frequency", type=float, default=1e-5, help="minimum frequency to keep below which values will be zeroed and all others renormalized to sum to one")
    args = parser.parse_args()

    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)

    frequencies = {
        node_name: region_frequencies["global"] if "global" in region_frequencies else region_frequencies["frequencies"]
        for node_name, region_frequencies in frequencies_json.items()
        if node_name not in ["pivots", "counts", "generated_by"]
    }

    # Collect the last frequency for each node keeping only terminal nodes
    # (tips) unless internal nodes are also requested.
    records = [
        {
            "strain": node.name,
            "frequency": float(frequencies[node.name][-1]),
            "is_terminal": node.is_terminal()
        }
        for node in tree.find_clades()
        if args.include_internal_nodes or node.is_terminal()
    ]

    # Convert frequencies data into a data frame.
    df = pd.DataFrame(records)

    # Replace records whose frequency values are below the requested minimum
    # with zeros and renormalize the remaining records to sum to one.
    to_zero = df["frequency"] < args.minimum_frequency
    not_to_zero = ~to_zero
    df.loc[to_zero, "frequency"] = 0.0
    df.loc[not_to_zero, "frequency"] = df.loc[not_to_zero, "frequency"] / df.loc[not_to_zero, "frequency"].sum()

    # Add any additional annotations requested by the user in the format of
    # "key=value" pairs where each key becomes a new column with the given
    # value.
    if args.annotations:
        for annotation in args.annotations:
            key, value = annotation.split("=")
            df[key] = value

    # Save the table.
    df.to_csv(args.output, sep="\t", float_format="%.6f", index=False, header=True)
