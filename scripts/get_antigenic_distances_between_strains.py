#!/usr/bin/env python3
import argparse
from augur.utils import annotate_parents_for_tree, read_node_data, read_tree
from collections import defaultdict
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--titer-model", required=True, help="titer model JSON with raw and normalized titers annotated in 'titers' key")
    parser.add_argument("--titers", required=True, help="raw titers table with information about the source of each titer")
    parser.add_argument("--tree", required=True, help="tree used to identify the given clades")
    parser.add_argument("--clades", required=True, help="clade annotations in a node data JSON")
    parser.add_argument("--branch-lengths", required=True, help="branch length annotations including `numdate` calculated by TreeTime")
    parser.add_argument("--frequencies", required=True, help="tip frequencies JSON from augur frequencies")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the output table in the format of 'key=value' pairs")
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

    # Load raw titer data, to get the original source for each tuple of test
    # virus/reference virus/ferret.
    raw_titers = pd.read_csv(
        args.titers,
        sep="\t",
        usecols=("virus_strain", "serum_strain", "serum_id", "source"),
    )

    # The source column starts with the original collaborating center (e.g.,
    # "cdc" or "crick_something"), so we split on the underscore and take the
    # first piece as the center.
    raw_titers["source"] = raw_titers["source"].apply(lambda source: source.split("_")[0])

    # Load titer model data.
    with open(args.titer_model, "r") as fh:
        titer_data = json.load(fh)

    titers = titer_data.pop("titers")
    potencies = titer_data.pop("potency")

    # Convert titer data to a data frame.
    titer_records = []
    for reference_strain, reference_titers in titers.items():
        reference_date = date_by_strain[reference_strain]
        reference_titer_index = (raw_titers["serum_strain"] == reference_strain)

        for test_strain, test_titers in reference_titers.items():
            test_date = date_by_strain[test_strain]
            test_frequency = current_frequency_by_strain[test_strain]
            test_titer_index = (raw_titers["virus_strain"] == test_strain)

            for serum, serum_titers in test_titers.items():
                log2_titer, raw_titer = serum_titers
                serum_index = (raw_titers["serum_id"] == serum)

                sources = raw_titers.loc[
                    reference_titer_index & test_titer_index & serum_index,
                    "source"
                ].values
                source = sources[0] if len(sources) > 0 else "unknown"

                titer_records.append({
                    "reference_strain": reference_strain,
                    "reference_date": reference_date,
                    "test_strain": test_strain,
                    "test_date": test_date,
                    "serum": serum,
                    "log2_titer": log2_titer,
                    "raw_titer": raw_titer,
                    "test_frequency": test_frequency,
                    "source": source,
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

    # Load tree.
    tree = read_tree(args.tree)
    tree = annotate_parents_for_tree(tree)

    # Load clade data.
    clades = read_node_data(args.clades)

    # Track all clade memberships in a new attribute to properly handle nested
    # clades.
    clades_by_name = defaultdict(set)
    for node in tree.find_clades():
        clades_by_name[node.name].add(clades["nodes"][node.name]["clade_membership"])
        if node.parent is not None:
            clades_by_name[node.name].update(
                clades_by_name[node.parent.name]
            )

    # Calculate clade frequencies.
    frequency_by_clade = defaultdict(float)
    for node in tree.find_clades(terminal=True):
        for clade in clades_by_name[node.name]:
            frequency_by_clade[clade] += current_frequency_by_strain[node.name]

    # Convert clade data to a data frame.
    clade_table = pd.DataFrame([
        {
            "strain": strain,
            "clade": strain_data["clade_membership"],
            "clade_frequency": frequency_by_clade[strain_data["clade_membership"]],
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
