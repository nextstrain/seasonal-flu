"""Annotate number of titers per node in node data JSON format for auspice.
"""
import argparse
from augur.titer_model import TiterCollection
from augur.utils import read_tree, write_json
import pandas as pd


def get_categorical_range_for_count(titer_count):
    if titer_count == 0:
        return "0"
    elif titer_count <= 5:
        return "1-5"
    elif titer_count <= 10:
        return "6-10"
    else:
        return ">10"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Newick tree with strains to annotate number of titers")
    parser.add_argument("--titers", help="tab-delimited file of raw titer measurements")
    parser.add_argument("--attribute-name", help="attribute name to use for number of titers per strain")
    parser.add_argument("--include-internal-nodes", action="store_true", help="calculate total measurements per internal node in addition to tips")
    parser.add_argument("--use-categorical-ranges", action="store_true", help="annotate nodes with categorical ranges of titer counts to enable better control over coloring in auspice")
    parser.add_argument("--output", help="JSON in node data format for use by augur export")

    args = parser.parse_args()

    tree = read_tree(args.tree)

    # Count titer measurements per test strain.
    titers, strains, sources = TiterCollection.load_from_file(args.titers)
    titer_count_by_strain = pd.Series([record[0] for record in titers.keys()]).value_counts().to_dict()

    # Annotate number of measurements per node such that internal nodes are
    # annotated with the sum of all measurements for their descendant tips.
    node_data = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            if node.name in titer_count_by_strain:
                node_data[node.name] = {
                    args.attribute_name: titer_count_by_strain[node.name]
                }
        elif args.include_internal_nodes:
            node_data[node.name] = {
                args.attribute_name: sum([
                    node_data[child.name][args.attribute_name]
                    for child in node.clades
                    if child.name in node_data
                ])
            }

    # Assign categorical counts. These ranges are hardcoded for now.
    if args.use_categorical_ranges:
        for node in tree.find_clades():
            if node.name in node_data:
                node_data[node.name][args.attribute_name] = get_categorical_range_for_count(
                    node_data[node.name][args.attribute_name]
                )
            else:
                # Get the categorical value for zero, if no counts are assigned to this node.
                node_data[node.name] = {
                    args.attribute_name: get_categorical_range_for_count(0)
                }

    # Save titers per strain in node data format.
    write_json({"nodes": node_data}, args.output)
