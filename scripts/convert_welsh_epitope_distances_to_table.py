import argparse
from augur.utils import read_node_data, read_tree
import csv
import json


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree", required=True, help="Newick tree used to identify leaf nodes in distances JSON")
    parser.add_argument("--distances", required=True, help="distances in node data JSON format")
    parser.add_argument("--distance-map", help="distance map with metadata stored in specific top-level attributes")
    parser.add_argument("--distance-map-attributes", nargs="*", help="optional list of top-level attributes to extract as metadata from the distance map")
    parser.add_argument("--output", required=True, help="distances in TSV format")

    args = parser.parse_args()

    # Load tree.
    tree = read_tree(args.tree)

    # Load distances.
    distances_per_node = read_node_data(args.distances)["nodes"]

    # Load distance map, if attributes requested.
    attributes = {}
    if args.distance_map and args.distance_map_attributes:
        with open(args.distance_map, encoding="utf-8") as fh:
            distance_map = json.load(fh)

        attributes = [
            distance_map[attribute]
            for attribute in args.distance_map_attributes
        ]

    # Write out distances in TSV format.
    with open(args.output, "w", encoding="utf-8") as oh:
        writer = csv.writer(oh, delimiter="\t", lineterminator="\n")
        writer.writerow(["strain", "welsh_escape_for_serum"] + args.distance_map_attributes)

        for node in tree.find_clades(terminal=True):
            record = [
                node.name,
                distances_per_node[node.name]["welsh_escape_for_serum"],
            ] + attributes
            writer.writerow(record)
