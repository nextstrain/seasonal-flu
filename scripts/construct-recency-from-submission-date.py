import argparse
from augur.utils import read_metadata, write_json
import datetime
import json
import numpy as np
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign each sequence a field that specifies when it was added",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--metadata',
        type=str,
        required=True,
        help="tab-delimited metadata file"
    )
    parser.add_argument(
        '--submission-date-field',
        default="date_submitted",
        help="field in the metadata with date that each sequences was submitted to its database"
    )
    parser.add_argument(
        '--date-bins',
        type=int,
        nargs="+",
        help="edges of bins to group sequences by days since they were submitted. Sequences submitted farther back than the maximum bin's number of days in the past will be assigned to an 'older' bin.",
        default=[7, 30, 120]
    )
    parser.add_argument(
        '--date-bin-labels',
        nargs="+",
        help="names to use for bins in node JSON output",
        default=["last week", "last month", "last quarter"]
    )
    parser.add_argument(
        '--upper-bin-label',
        help="name to use for the upper bin for sequences that were submitted farther back than the maximum given bin number of days in the past.",
        default="older"
    )
    parser.add_argument(
        '--output-field-name',
        default="recency",
        help="name of the field to store recency information in the node data JSON"
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help="node data JSON with recency annotations"
    )
    args = parser.parse_args()

    meta, columns = read_metadata(args.metadata)
    meta = pd.DataFrame.from_dict(meta, "index")
    meta[args.submission_date_field] = pd.to_datetime(meta[args.submission_date_field])

    node_data = {'nodes': {}}

    # Calculate days since sequences were submitted.
    today = pd.to_datetime(datetime.date.today())
    meta["_days_since_submission"] = (today - meta[args.submission_date_field]).dt.days

    # Create bins to use for day intervals.
    bins = args.date_bins

    # Bins need to start with zero.
    if 0 not in bins:
        bins.insert(0, 0)

    # The last bin needs to include the maximum possible value.
    bins.append(np.inf)

    # Build a list of bin labels.
    bin_labels = args.date_bin_labels
    bin_labels.append(args.upper_bin_label)

    # Bin sequences by relevant submission delay intervals.
    meta["_day_bins"] = pd.cut(
        meta["_days_since_submission"],
        bins=bins,
        labels=bin_labels,
        include_lowest=True
    )

    # Create node data annotations of recency per strain.
    recency_by_strain = meta["_day_bins"].to_dict()
    for strain, recency in recency_by_strain.items():
        node_data['nodes'][strain] = {args.output_field_name: recency}

    write_json(node_data, args.output)
