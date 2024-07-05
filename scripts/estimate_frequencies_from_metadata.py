#!/usr/bin/env python3
import argparse
import numpy as np

from augur.dates import get_numerical_dates, numeric_date_type
from augur.frequencies import format_frequencies
from augur.frequency_estimators import get_pivots, KdeFrequencies
from augur.io import read_metadata
from augur.utils import write_json


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Estimate sequence frequencies from metadata with collection dates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--metadata", required=True, help="TSV file of metadata with at least 'strain' and 'date' columns")
    parser.add_argument("--narrow-bandwidth", required=True, type=float, help="narrow bandwidth for KDE frequencies")
    parser.add_argument("--proportion-wide", type=float, default=0.0, help="proportion of wide bandwidth to use for KDE frequencies")
    parser.add_argument("--pivot-interval", type=int, default=4, help="interval between pivots in weeks")
    parser.add_argument("--min-date", type=numeric_date_type, help="minimum date to estimate frequencies for")
    parser.add_argument("--max-date", type=numeric_date_type, help="maximum date to estimate frequencies for")
    parser.add_argument("--output", required=True, help="JSON file in tip-frequencies format")
    args = parser.parse_args()

    columns_to_load = ["strain", "date"]
    metadata = read_metadata(
        args.metadata,
        columns=columns_to_load,
        dtype="string",
    )
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')

    observations = [
        np.mean(dates[strain])
        for strain in metadata.index.values
    ]
    pivots = get_pivots(
        observations,
        args.pivot_interval,
        args.min_date,
        args.max_date,
        "weeks",
    )

    frequencies = KdeFrequencies(
        sigma_narrow=args.narrow_bandwidth,
        proportion_wide=args.proportion_wide,
        pivot_frequency=args.pivot_interval,
        start_date=args.min_date,
        end_date=args.max_date,
    )
    frequency_matrix = frequencies.estimate_frequencies(
        observations,
        pivots,
    )
    frequencies.frequencies = {
        strain: frequency_matrix[index]
        for index, strain in enumerate(metadata.index.values)
        if frequency_matrix[index].sum() > 0
    }
    frequencies.pivots = pivots

    frequency_dict = frequencies.to_json()
    write_json(frequency_dict, args.output)
