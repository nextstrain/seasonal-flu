#!/usr/bin/env python3
import argparse
import json
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--auspice-config", required=True, help="Auspice config JSON with coloring entry to have scale added to")
    parser.add_argument("--coloring-field", required=True, help="name of the coloring field in the Auspice config JSON")
    parser.add_argument("--values", nargs="+", required=True, help="list of values to assign colors to")
    parser.add_argument("--color-schemes", required=True, help="file with color schemes with N tab-delimited colors on row N")
    parser.add_argument("--output", required=True, help="Auspice config JSON with scale added to the requested coloring")
    args = parser.parse_args()

    n_colors = len(args.values)
    sorted_values = sorted(args.values)

    # Load and validate the input Auspice config JSON.
    with open(args.auspice_config, "r", encoding="utf-8") as fh:
        auspice_config = json.load(fh)

    if not "colorings" in auspice_config:
        print(f"The Auspice config JSON, '{args.auspice_config}', does not have a 'colorings' entry.", file=sys.stderr)
        sys.exit(1)

    if not any(args.coloring_field == coloring["key"] for coloring in auspice_config["colorings"]):
        print(f"The Auspice config JSON, '{args.auspice_config}', does not have the requested coloring field, '{args.coloring_field}'.", file=sys.stderr)
        sys.exit(1)

    # Load the required number of colors.
    with open(args.color_schemes, "r", encoding="utf-8") as fh:
        for line in fh:
            if len(line.split("\t")) == n_colors:
                colors = [color.strip() for color in line.split("\t")]
                break

    # Add the color scale to the config.
    color_scale = [
        [value, color]
        for value, color in zip(sorted_values, colors)
    ]
    for coloring in auspice_config["colorings"]:
        if coloring["key"] == args.coloring_field:
            coloring["scale"] = color_scale
            break

    # Save the new Auspice config JSON.
    with open(args.output, "w", encoding="utf-8") as oh:
        json.dump(auspice_config, oh, indent=2)
