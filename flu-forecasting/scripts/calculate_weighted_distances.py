"""Calculate weighted distances to the future between a given set of tip attributes and a given set of forecasted tip attributes.
"""
import argparse
from augur.utils import write_json
import numpy as np
import pandas as pd
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate pairwise distances between samples at adjacent timepoints (t and t - delta months)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tip-attributes", required=True, help="a tab-delimited file describing tip attributes at one or more timepoints and including sequences")
    parser.add_argument("--forecasts", required=True, help="a tab-delimited file describing forecasts for tips and their sequences")
    parser.add_argument("--output", help="node data JSON for the tips in the given tip attributes table", required=True)
    parser.add_argument("--distance-attribute-name", default="weighted_distance_to_future", help="attribute name for weighted distances")
    parser.add_argument("--sequence-attribute-name", default="aa_sequence", help="attribute name of sequences to compare")
    args = parser.parse_args()

    # Load tip attributes.
    tips = pd.read_csv(args.tip_attributes, sep="\t", parse_dates=["timepoint"])
    if args.sequence_attribute_name not in tips.columns:
        print("ERROR: missing sequence column '%s' in tip attributes file '%s'" % (args.sequence_attribute_name, args.tip_attributes), file=sys.stderr)
        sys.exit(1)

    tip_sequence_by_name = dict(tips.loc[:, ["strain", args.sequence_attribute_name]].values)

    # Load forecasts.
    forecasts = pd.read_csv(args.forecasts, sep="\t")
    if args.sequence_attribute_name not in forecasts.columns:
        print("ERROR: missing sequence column '%s' in forecasts file '%s'" % (args.sequence_attribute_name, args.forecasts), file=sys.stderr)
        sys.exit(1)

    future_tip_sequence_by_name = dict(forecasts.loc[:, ["strain", args.sequence_attribute_name]].values)
    future_tip_frequency_by_name = dict(forecasts.loc[:, ["strain", "projected_frequency"]].values)

    # Convert future tip sequences as arrays once for pairwise comparisons.
    for tip_name in future_tip_sequence_by_name.keys():
        future_tip_sequence_by_name[tip_name] = np.frombuffer(
            future_tip_sequence_by_name[tip_name].encode(),
            dtype="S1"
        )

    # Calculate weighted distances between given tips and forecasts and store in
    # node data JSON format.
    distances = {}
    for tip_name, tip_sequence in tip_sequence_by_name.items():
        current_tip_sequence_array = np.frombuffer(tip_sequence.encode(), dtype="S1")
        weighted_distance_to_future = 0.0

        for future_tip_name in future_tip_sequence_by_name.keys():
            distance = (current_tip_sequence_array != future_tip_sequence_by_name[future_tip_name]).sum()
            weighted_distance_to_future += future_tip_frequency_by_name[future_tip_name] * distance

        distances[tip_name] = {args.distance_attribute_name: weighted_distance_to_future}

    # Export distances to JSON.
    write_json({"nodes": distances}, args.output)
