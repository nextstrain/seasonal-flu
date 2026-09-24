#!/usr/bin/env python3
import argparse
from augur.utils import read_node_data
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--titer-model", required=True, help="node data JSON from titer model with inferred titers annotated in 'nodes' key by field ending with 'cTiterSub'")
    parser.add_argument("--titers", required=True, help="TSV of titers used to fit the given model")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the output table in the format of 'key=value' pairs")
    parser.add_argument("--output", required=True, help="table of antigenic distances in log2 titers between reference and test strains")

    args = parser.parse_args()

    # Load raw titers to get the reference name.
    raw_titers = pd.read_csv(
        args.titers,
        sep="\t",
        nrows=2,
    )
    reference = raw_titers["serum_strain"].values[0]

    # Load titer model data.
    titer_data = read_node_data(args.titer_model)["nodes"]

    # Convert titer data to a data frame.
    titer_records = []
    for test_strain, test_strain_values in titer_data.items():
        for key, value in test_strain_values.items():
            if key.endswith("cTiterSub"):
                titer_records.append({
                    "reference_strain": reference,
                    "test_strain": test_strain,
                    "log2_titer": value,
                })

    titer_table = pd.DataFrame(titer_records)

    # Add any additional annotations requested by the user in the format of
    # "key=value" pairs where each key becomes a new column with the given
    # value.
    if args.annotations:
        for annotation in args.annotations:
            key, value = annotation.split("=")
            titer_table[key] = value

    # Save the annotated table.
    titer_table.to_csv(
        args.output,
        sep="\t",
        index=False,
        float_format="%.4f"
    )
