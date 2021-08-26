#!/usr/bin/env python3
import argparse
from augur.io import read_metadata
from augur.utils import read_strains, write_json
import epiweeks
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        usage="Calculate epiweeks for dates in the given metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--metadata", required=True, help="metadata with a 'date' column")
    parser.add_argument("--strains", required=True, help="a list of strains to annotate epiweeks for using the given metadata")
    parser.add_argument("--attribute-name", default="epiweek", help="name to store annotations of epiweeks in JSON output")
    parser.add_argument("--output-node-data", required=True, help="node data JSON with epiweek annotations")

    args = parser.parse_args()

    # Read metadata.
    metadata = read_metadata(args.metadata)

    # Read strains.
    strains = read_strains(args.strains)

    # Find metadata for requested strains.
    # Use `metadata.index.isin` in case listed strain does not exist in metadata
    metadata_for_strains = metadata.loc[metadata.index.isin(strains)].copy()

    # Find records with unambiguous dates.
    metadata_with_dates = metadata_for_strains.loc[~metadata["date"].str.contains("X"), ["date"]].copy()

    # Convert date strings to timestamps.
    metadata_with_dates["date"] = pd.to_datetime(metadata_with_dates["date"])

    # Calculate epiweeks from date objects as a new annotation.
    metadata_with_dates["epiweek"] = metadata_with_dates["date"].apply(lambda date: epiweeks.Week.fromdate(date).cdcformat())

    # Create a node data object with epiweeks.
    node_data = {}
    for strain, record in metadata_with_dates.to_dict(orient="index").items():
        node_data[strain] = {
            args.attribute_name: record["epiweek"],
        }

    # Save node data.
    write_json({"nodes": node_data}, args.output_node_data)
