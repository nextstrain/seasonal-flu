"""
Calculate pairwise distances between samples at adjacent timepoints (t and t - delta months).
"""
import argparse
import numpy as np
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate pairwise distances between samples at adjacent timepoints (t and t - delta months)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--tip-attributes", required=True, help="a tab-delimited file describing tip attributes at one or more timepoints")
    parser.add_argument("--delta-months", required=True, nargs="+", type=int, help="number of months between timepoints to be compared")
    parser.add_argument("--output", help="tab-delimited file of pairwise distances between tips in timepoints separate by the given delta time", required=True)
    parser.add_argument("--sequence-attribute-name", default="aa_sequence", help="attribute name of sequences to compare")
    args = parser.parse_args()

    # Load tip attributes.
    tips = pd.read_csv(args.tip_attributes, sep="\t", parse_dates=["timepoint"])

    # Calculate pairwise distances between all tips within a timepoint and at
    # the next timepoint as defined by the given delta months.
    distances = []

    for timepoint, timepoint_df in tips.groupby("timepoint"):
        current_tips = [
            tuple(values)
            for values in timepoint_df.loc[:, ["strain", args.sequence_attribute_name]].values.tolist()
        ]
        comparison_tips = current_tips

        for delta_month in args.delta_months:
            future_timepoint_df = tips[tips["timepoint"] == (timepoint + pd.DateOffset(months=delta_month))]
            future_tips = [
                tuple(values)
                for values in future_timepoint_df.loc[:, ["strain", args.sequence_attribute_name]].values.tolist()
            ]
            comparison_tips = comparison_tips + future_tips

        comparison_tips = list(set(comparison_tips))
        for current_tip, current_tip_sequence in current_tips:
            current_tip_sequence_array = np.frombuffer(current_tip_sequence.encode(), dtype="S1")

            for future_tip, future_tip_sequence in comparison_tips:
                future_tip_sequence_array = np.frombuffer(future_tip_sequence.encode(), dtype="S1")
                distance = (current_tip_sequence_array != future_tip_sequence_array).sum()
                distances.append({
                    "sample": current_tip,
                    "other_sample": future_tip,
                    "distance": distance
                })

    # Save distances table.
    distances = pd.DataFrame(distances)
    distances.to_csv(args.output, sep="\t", index=False)
