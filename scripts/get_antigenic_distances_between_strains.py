#!/usr/bin/env python3
import argparse
from augur.utils import read_node_data
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--titers", required=True, help="titer model JSON with raw and normalized titers annotated in 'titers' key")
    parser.add_argument("--clades", required=True, help="clade annotations in a node data JSON")
    parser.add_argument("--branch-lengths", required=True, help="branch length annotations including `numdate` calculated by TreeTime")
    parser.add_argument("--frequencies", required=True, help="tip frequencies JSON from augur frequencies")
    parser.add_argument("--output", required=True, help="table of antigenic distances in log2 titers between reference and test strains with clade annotations")
    args = parser.parse_args()

    # Load branch lengths.
    branch_lengths = read_node_data(args.branch_lengths)["nodes"]
    date_by_strain = {
        strain: node_data["numdate"]
        for strain, node_data in branch_lengths.items()
    }

    # Load frequencies and get the most recent frequency per strain.
    with open(args.frequencies) as fh:
        frequencies_data = json.load(fh)

    current_frequency_by_strain = {
        strain: strain_data["frequencies"][-1]
        for strain, strain_data in frequencies_data.items()
        if "frequencies" in strain_data
    }

    # Load titer data.
    with open(args.titers, "r") as fh:
        titer_data = json.load(fh)

    titers = titer_data.pop("titers")
    potencies = titer_data.pop("potency")

    # Convert titer data to a data frame.
    titer_records = []
    for reference_strain, reference_titers in titers.items():
        reference_date = date_by_strain[reference_strain]

        for test_strain, test_titers in reference_titers.items():
            test_date = date_by_strain[test_strain]
            test_frequency = current_frequency_by_strain[test_strain]

            for serum, serum_titers in test_titers.items():
                log2_titer, raw_titer = serum_titers

                titer_records.append({
                    "reference_strain": reference_strain,
                    "reference_date": reference_date,
                    "test_strain": test_strain,
                    "test_date": test_date,
                    "serum": serum,
                    "log2_titer": log2_titer,
                    "raw_titer": raw_titer,
                    "test_frequency": test_frequency
                })

    titer_table = pd.DataFrame(titer_records)

    # Convert potencies to a data frame.
    potency_table = pd.DataFrame([
        {
            "reference_strain": strain,
            "potency": strain_potencies["mean_potency"]
        }
        for strain, strain_potencies in potencies.items()
    ])

    # Annotate titers with potencies.
    titer_table = titer_table.merge(
        potency_table,
        on="reference_strain"
    )

    # Load clade data.
    clades = read_node_data(args.clades)

    # Convert clade data to a data frame.
    clade_table = pd.DataFrame([
        {
            "strain": strain,
            "clade": strain_data["clade_membership"]
        }
        for strain, strain_data in clades["nodes"].items()
        if not strain.startswith("NODE")
    ])

    # Annotate titers with clades for reference and test strains.
    # First, annotate reference clades.
    titer_table = titer_table.merge(
        clade_table,
        left_on="reference_strain",
        right_on="strain"
    ).drop(
        columns=["strain"]
    )

    # Then, annotate test clades.
    titer_table = titer_table.merge(
        clade_table,
        left_on="test_strain",
        right_on="strain",
        suffixes=["_reference", "_test"]
    ).drop(
        columns=["strain"]
    )

    # Save the annotated table.
    titer_table.to_csv(
        args.output,
        sep="\t",
        index=False,
        float_format="%.4f"
    )
