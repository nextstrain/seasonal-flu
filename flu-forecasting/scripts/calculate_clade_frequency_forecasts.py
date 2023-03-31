#!/usr/bin/env python3
import argparse
import pandas as pd

from augur.utils import read_node_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--forecasts", required=True, help="TSV of forecasts per strain")
    parser.add_argument("--clades", required=True, help="node data JSON of clades per strain")
    parser.add_argument("--output", required=True, help="TSV of forecasts per clade")

    args = parser.parse_args()

    # Load forecasts.
    forecasts = pd.read_csv(
        args.forecasts,
        sep="\t",
        usecols=("timepoint", "strain", "projected_frequency"),
    )

    # Load clades.
    clades = read_node_data(args.clades)
    clade_by_strain = {
        name: data["clade_membership"]
        for name, data in clades["nodes"].items()
    }

    # Assign clades to strains in the forecasts table.
    forecasts["clade"] = forecasts["strain"].map(clade_by_strain)

    # Calculate projected frequency per clade.
    forecasts_by_clade = forecasts.groupby(["timepoint", "clade"])["projected_frequency"].sum().reset_index()

    # Save forecasts by clade.
    forecasts_by_clade.to_csv(
        args.output,
        sep="\t",
        header=True,
        index=False,
    )
