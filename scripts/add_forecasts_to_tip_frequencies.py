"""Add forecasts to tip frequencies.
"""
import argparse
import json
import numpy as np
import pandas as pd
from treetime.utils import numeric_date


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add forecasts to tip frequencies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--frequencies", required=True, help="an auspice tip frequencies JSON")
    parser.add_argument("--forecasts", required=True, help="a tab-delimited file containing columns for 'strain' and 'projected_frequency'")
    parser.add_argument("--output", required=True, help="an auspice tip frequencies JSON with forecast frequencies and their pivots added")
    args = parser.parse_args()

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies = json.load(fh)

    # Load forecasts.
    forecasts = pd.read_csv(args.forecasts, sep="\t", parse_dates=["timepoint", "future_timepoint"])

    # Annotate numeric dates as pivots.
    forecasts["pivot"] = forecasts["future_timepoint"].apply(lambda date: round(numeric_date(date), 2))
    forecasts = forecasts.sort_values(["strain", "pivot"]).copy()

    # Get distinct new pivots and add these to the existing frequencies.
    new_pivots = forecasts["pivot"].unique().tolist()
    frequencies["pivots"].extend(new_pivots)

    # Annotate new frequencies for all strains.
    for strain, pivot, projected_frequency in forecasts.loc[:, ["strain", "pivot", "projected_frequency"]].values:
        frequencies[strain]["frequencies"].append(projected_frequency)

    # Add zero frequencies for all strains that were not included in forecasts.
    for strain in frequencies.keys():
        if strain != "pivots" and len(frequencies[strain]["frequencies"]) < len(frequencies["pivots"]):
            frequencies[strain]["frequencies"].extend(list(np.zeros_like(new_pivots)))

    # Annotate the projection pivot.
    projection_pivot = forecasts["timepoint"].drop_duplicates().apply(lambda date: round(numeric_date(date), 2)).values[0]
    frequencies["projection_pivot"] = projection_pivot

    # Save projected frequencies and pivots.
    with open(args.output, "w") as oh:
        json.dump(frequencies, oh, indent=1)
