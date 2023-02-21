#!/usr/bin/env python3
import argparse
from augur.validate import measurements_collection_config
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--collection", required=True, help="collection TSV that will be passed to augur measurements export")
    parser.add_argument("--groupings", nargs="+", required=True, help="grouping columns in the given TSV for which an order will be generated")
    parser.add_argument("--fields", nargs="+", help="fields from the TSV file to display in the tooltip in the order given to this argument")
    parser.add_argument("--output", required=True, help="configuration JSON for measurements export with ordering specified for the requested grouping")

    args = parser.parse_args()

    # Read collection.
    collection_df = pd.read_csv(args.collection, sep="\t", usecols=args.groupings + ["reference_date", "clade_reference"])

    # Find earliest reference strain within each clade.
    min_reference_date_by_reference_clade = collection_df.groupby("clade_reference")["reference_date"].min().reset_index().rename(
        columns={"reference_date": "min_reference_date"}
    )

    # Annotate earliest reference date to collection.
    collection_df = collection_df.merge(
        min_reference_date_by_reference_clade,
        on="clade_reference",
        how="left",
    )

    # Sort collection by earliest reference date and clade.
    sorted_df = collection_df.sort_values(
        ["min_reference_date", "clade_reference", "reference_date"],
        ascending=[False, True, False],
    )

    # Extract sorted values for each grouping column.
    groupings_config = []
    for grouping in args.groupings:
        # Skip sorting of "source" values by clade, since it doesn't make sense.
        # Omitting the "order" key from its config also demonstrates how the
        # default ordering by record count works.
        if grouping == "source":
            groupings_config.append({
                "key": grouping,
            })
        else:
            sorted_grouping_values = sorted_df[grouping].drop_duplicates().tolist()
            groupings_config.append({
                "key": grouping,
                "order": sorted_grouping_values,
            })

    # Build the configuration JSON entry for the requested grouping.
    config = {
        "groupings": groupings_config,
    }

    # Specify fields and their order, if given.
    if args.fields:
       config["fields"] = [
           {"key": field}
           for field in args.fields
       ]

    # Save configuration JSON.
    with open(args.output, "w", encoding="utf-8") as oh:
        json.dump(config, oh)

    # Validate configuration JSON.
    measurements_collection_config(args.output)
