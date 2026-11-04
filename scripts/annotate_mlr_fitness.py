#!/usr/bin/env python3
import argparse

from augur.io import read_metadata, write_json
from augur.utils import read_tree
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree", required=True, help="Newick file of a tree with named internal nodes for which antigenic advance should be calculated")
    parser.add_argument("--metadata", required=True, help="TSV file of strain metadata with haplotypes in `subclade_haplotype_ha` field")
    parser.add_argument("--fitnesses", required=True, nargs="+", help="TSV file(s) of haplotype fitnesses with haplotypes in the `variant` field, fitnesses in the `mean` field, and a `location` field with a 'hierarchical' entry")
    parser.add_argument("--variant-columns", required=True, nargs="+", help="names of the columns in the metadata with variant classifications corresponding to the variants in the fitness files provided. Column names should be in the same order as the corresponding fitness TSV files.")
    parser.add_argument("--output-node-data", required=True, help="node data JSON file with MLR fitness per node")

    args = parser.parse_args()

    tree = read_tree(args.tree)

    metadata = read_metadata(args.metadata)

    # Create a mapping of each variant column (e.g., "subclade_haplotype_ha")
    # and the fitnesses per variant in that column. Each variant column
    # corresponds to its own fitness table from the MLR.
    fitness_by_variant_column_and_name = {}
    for variant_column, fitness_path in zip(args.variant_columns, args.fitnesses):
        fitness_by_variant_column_and_name[variant_column] = dict(
            pd.read_csv(
                fitness_path,
                sep="\t",
            ).query("location == 'hierarchical'").loc[
                :,
                ["variant", "mean"]
            ].values
        )

    # For each tip in the tree and for each variant column, find the variant for
    # the tip and its corresponding fitness. Do not create annotations for tips
    # whose variant lacks a fitness annotation, since those variants were not
    # present in the corresponding MLR model.
    node_data = {}
    for node in tree.find_clades(terminal=True):
        for variant_column in args.variant_columns:
            variant = metadata.at[node.name, variant_column]

            fitness = fitness_by_variant_column_and_name.get(
                variant_column,
                {}
            ).get(
                variant
            )

            if fitness is not None:
                if node.name not in node_data:
                    node_data[node.name] = {}

                node_data[node.name][f"mlr_fitness_{variant_column}"] = fitness

    write_json({"nodes": node_data}, args.output_node_data)
