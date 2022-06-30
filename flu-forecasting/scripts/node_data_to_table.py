"""
Convert one or more augur node data JSONs into a single table of values labelled by tip or internal node status in a tree.
"""
import argparse
from augur.utils import read_node_data
import Bio.Phylo
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert node data JSONs to a data frame",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick file for the tree used to construct the given node data JSONs")
    parser.add_argument("--metadata", help="file with metadata associated with viral sequences, one for each segment")
    parser.add_argument("--jsons", nargs="+", required=True, help="node data JSON(s) from augur")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the output table in the format of 'key=value' pairs")
    parser.add_argument("--excluded-fields", nargs="+", help="names of columns to omit from output table")
    parser.add_argument("--output", required=True, help="tab-delimited file collecting all given node data")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data associated with internal nodes in the output table")
    args = parser.parse_args()

    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Load metadata for samples.
    metadata = pd.read_csv(args.metadata, sep="\t")

    # Load one or more node data JSONs into a single dictionary indexed by node name.
    node_data = read_node_data(args.jsons)

    # Convert node data into a data frame.
    # Data are initially loaded with one column per node.
    # Transposition converts the table to the expected one row per node format.
    df = pd.DataFrame(node_data["nodes"]).T.rename_axis("strain").reset_index()

    # Annotate node data with per sample metadata.
    df = df.merge(metadata, on="strain", suffixes=["", "_metadata"])

    # Remove excluded fields if they are in the data frame.
    if args.excluded_fields:
        df = df.drop(columns=[field for field in args.excluded_fields if field in df.columns])

    # Annotate the tip/internal status of each node using the tree.
    node_terminal_status_by_name = {node.name: node.is_terminal() for node in tree.find_clades()}
    df["is_terminal"] = df["strain"].map(node_terminal_status_by_name)

    # Eliminate internal nodes if they have not been requested.
    if not args.include_internal_nodes:
        df = df[df["is_terminal"]].copy()

    # Add any additional annotations requested by the user in the format of
    # "key=value" pairs where each key becomes a new column with the given
    # value.
    if args.annotations:
        for annotation in args.annotations:
            key, value = annotation.split("=")
            df[key] = value

    # Save the table.
    df.to_csv(args.output, sep="\t", index=False, header=True)
