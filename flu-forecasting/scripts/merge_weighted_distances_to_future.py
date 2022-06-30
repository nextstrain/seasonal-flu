#!/usr/bin/env python3

import argparse
import json


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--distances", nargs="+", required=True, help="node data JSON files of distances to be merged into a single output node data JSON")
    parser.add_argument("--output", required=True, help="merged node data JSON file")

    args = parser.parse_args()

    # Start with a single base JSON in node data format.
    with open(args.distances[0], "r", encoding="utf-8") as fh:
        base_json = json.load(fh)

    # Update the base JSON with each subsequent model's distances to the future.
    for json_file in args.distances[1:]:
        with open(json_file, "r", encoding="utf-8") as fh:
            model_json = json.load(fh)
            for strain, distances in model_json["nodes"].items():
                base_json["nodes"][strain].update(distances)

    # Save merged data.
    with open(args.output, "w", encoding="utf-8") as oh:
        json.dump(base_json, oh)
