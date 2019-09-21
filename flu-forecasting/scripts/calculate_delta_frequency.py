"""Calculate the change in frequency for clades over time (aka the delta frequency or dfreq).
"""
import argparse
from augur.frequency_estimators import TreeKdeFrequencies
from augur.utils import write_json
import Bio.Phylo
from collections import defaultdict
import json
import numpy as np


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the change in frequency for clades over time",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick tree")
    parser.add_argument("--frequencies", required=True, help="frequencies JSON")
    parser.add_argument("--frequency-method", required=True, choices=["kde", "diffusion"], help="method used to estimate frequencies")
    parser.add_argument("--clades", help="JSON of clade annotations for nodes in the given tree")
    parser.add_argument("--delta-pivots", type=int, default=1, help="number of frequency pivots to look back in time for change in frequency calculation")
    parser.add_argument("--output", required=True, help="JSON of delta frequency annotations for nodes in the given tree")

    args = parser.parse_args()

    # Load the tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)

    if args.frequency_method == "kde":
        kde_frequencies = TreeKdeFrequencies.from_json(frequencies_json)
        frequencies = kde_frequencies.frequencies

        # Load clades.
        with open(args.clades, "r") as fh:
            clades_json = json.load(fh)

        clades_by_node = {
            key: value["clade_membership"]
            for key, value in clades_json["nodes"].items()
        }

        # Calculate the total frequency per clade at the most recent timepoint and
        # requested timepoint in the past using non-zero tip frequencies.
        current_clade_frequencies = defaultdict(float)
        previous_clade_frequencies = defaultdict(float)

        for tip in tree.find_clades(terminal=True):
            # Add tip to current clade frequencies.
            current_clade_frequencies[clades_by_node[tip.name]] += frequencies[tip.name][-1]

            # Add tip to previous clade frequencies.
            previous_clade_frequencies[clades_by_node[tip.name]] += frequencies[tip.name][-(args.delta_pivots + 1)]

        # Determine the total time that elapsed between the current and past timepoint.
        delta_time = kde_frequencies.pivots[-1] - kde_frequencies.pivots[-(args.delta_pivots + 1)]

        # Calculate the change in frequency over time elapsed for each clade.
        delta_frequency_by_clade = {}
        for clade, current_frequency in current_clade_frequencies.items():
            # If the current clade was not observed in the previous timepoint, it
            # will have a zero frequency.
            delta_frequency_by_clade[clade] = (current_frequency - previous_clade_frequencies.get(clade, 0.0)) / delta_time

        # Assign clade delta frequencies to all corresponding tips and internal nodes.
        delta_frequency = {}
        for node in tree.find_clades(terminal=True):
            delta_frequency[node.name] = {
                "delta_frequency": delta_frequency_by_clade.get(clades_by_node[node.name], 0.0)
            }
    else:
        frequencies = frequencies_json

        # Determine the total time that elapsed between the current and past timepoint.
        delta_time = frequencies["pivots"][-1] - frequencies["pivots"][-(args.delta_pivots + 1)]

        delta_frequency = {}
        for node in tree.find_clades(terminal=True):
            delta_frequency[node.name] = {
                "delta_frequency": (frequencies[node.name]["global"][-1] - frequencies[node.name]["global"][-(args.delta_pivots + 1)]) / delta_time
            }

    # Write out the node annotations.
    write_json({"nodes": delta_frequency}, args.output)
