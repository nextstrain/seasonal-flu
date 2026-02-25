#!/usr/bin/env python3
import argparse

from augur.reconstruct_sequences import load_alignments
from augur.titer_model import SubstitutionModel
from augur.io import write_json
from augur.utils import read_tree
from collections import defaultdict
import numpy as np
from numpyro.diagnostics import hpdi
import pandas as pd


def annotate_tree(titer_model, tree, substitution_effect, n_samples, groups):
    """Annotates antigenic advance attributes to nodes of a given tree built from
    the same sequences used to train the model.
    """
    posterior_dTiterSub = {
        tree.root.name: np.zeros((n_samples, len(groups))),
    }
    posterior_cTiterSub = {
        tree.root.name: np.zeros((n_samples, len(groups))),
    }

    for node in tree.find_clades():
        for child in node.clades:
            # Get mutations between the current node and its parent.
            mutations = titer_model.get_mutations(node.name, child.name)

            # Calculate titer drop on the branch to the current node.
            posterior_dTiterSub[child.name] = np.zeros_like(posterior_dTiterSub[tree.root.name])
            for gene, mutation in mutations:
                gene_mutation = f"{gene}:{mutation}"
                if gene_mutation in substitution_effect:
                    for group_index, group in enumerate(groups):
                        posterior_dTiterSub[child.name][:, group_index] += substitution_effect[gene_mutation][group]

            # Calculate the cumulative titer drop from the root to the current node.
            posterior_cTiterSub[child.name] = posterior_cTiterSub[node.name] + posterior_dTiterSub[child.name]

    return posterior_dTiterSub, posterior_cTiterSub


def main(args):
    """Calculate antigenic advance.
    """
    # Load substitution weights.
    substitution_effects = pd.read_parquet(args.substitution_weights)

    group_by = ["substitution"]
    if args.group_by:
        group_by.append(args.group_by)

    if len(group_by) == 1:
        group_by = group_by[0]

    effects_by_substitution = defaultdict(dict)
    distinct_groups = set()
    for group_key, group_df in substitution_effects.groupby(group_by, sort=False):
        if args.group_by:
            effects_by_substitution[group_key[0]][group_key[1]] = group_df["value"].values
            distinct_groups.add(group_key[1])
        else:
            effects_by_substitution[group_key]["all"] = group_df["value"].values
            distinct_groups.add("all")

    distinct_groups = sorted(distinct_groups)

    n_samples = substitution_effects["posterior_sample"].drop_duplicates().shape[0]

    # Load tree.
    tree = read_tree(args.tree)

    # Load alignments.
    alignments = load_alignments(args.alignment, args.gene_names)

    # Setup a titer model with alignments but no titer data, so we can find
    # substitutions between nodes.
    titer_model = SubstitutionModel(alignments, {})

    # Annotate tree with effects per substitution.
    posterior_dTiterSub, posterior_cTiterSub  = annotate_tree(
        titer_model,
        tree,
        effects_by_substitution,
        n_samples,
        distinct_groups,
    )

    node_data = {}
    for node in tree.find_clades():
        node_data[node.name] = {}

        # Weights per branch (dTiterSub) and cumulative from the root
        # (cTiterSub) are stored as a matrix per node with posterior samples in
        # rows and serum groups in columns. The default serum group is "all"
        # when no group has been requested for the given weights input.
        # Calculate the mean and HPDIs across all posterior samples for each
        # serum group.
        dTiterSub = np.mean(posterior_dTiterSub[node.name], axis=0)
        dTiterSubLowerHPDI, dTiterSubUpperHPDI = hpdi(posterior_dTiterSub[node.name], 0.95, axis=0)

        cTiterSub = np.mean(posterior_cTiterSub[node.name], axis=0)
        cTiterSubLowerHPDI, cTiterSubUpperHPDI = hpdi(posterior_cTiterSub[node.name], 0.95, axis=0)

        # Save each serum group's mean and HPDIs in its own set of attributes
        # for visualization in Auspice.
        for group_index, group in enumerate(distinct_groups):
            # When we only have the default group, don't use the group name.
            # Otherwise, include the user-defined attribute prefix followed by
            # the serum group name.
            if group == "all" and len(distinct_groups) == 1:
                prefix = args.attribute_prefix
            else:
                prefix = f"{args.attribute_prefix}{group}_"

            node_data[node.name].update({
                f"{prefix}dTiterSub": dTiterSub[group_index],
                f"{prefix}dTiterSubLowerHPDI": dTiterSubLowerHPDI[group_index],
                f"{prefix}dTiterSubUpperHPDI": dTiterSubUpperHPDI[group_index],
                f"{prefix}cTiterSub": cTiterSub[group_index],
                f"{prefix}cTiterSubLowerHPDI": cTiterSubLowerHPDI[group_index],
                f"{prefix}cTiterSubUpperHPDI": cTiterSubUpperHPDI[group_index],
            })

    # Calculate mean substitution effects for the "all" group.
    mean_effect_by_substitution = {
        substitution: effects_by_substitution[substitution]["all"].mean().round(4)
        for substitution in effects_by_substitution.keys()
    }

    # Export the antigenic advance per node.
    write_json(
        {
            "nodes": node_data,
            "substitution": mean_effect_by_substitution,
        },
        args.output_node_data
    )

    # Optionally, export the posterior distributions of antigenic advance values per node.
    if args.output_posterior:
        posterior = []
        posterior_sample_indexes = np.arange(0, posterior_cTiterSub[tree.root.name].shape[0])
        for node in tree.find_clades():
            for group_index, group in enumerate(distinct_groups):
                posterior.append(
                    pd.DataFrame({
                        "node": node.name,
                        "group": group,
                        "posterior_sample": posterior_sample_indexes,
                        f"{args.attribute_prefix}cTiterSub": posterior_cTiterSub[node.name][:, group_index],
                    })
                )

        posterior = pd.concat(posterior, ignore_index=True)
        posterior.to_parquet(args.output_posterior)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--substitution-weights", required=True, help="parquet file with weights per substitution and posterior sample defined in columns named 'posterior_sample', 'substitution', and 'value'")
    parser.add_argument("--tree", required=True, help="Newick file of a tree with named internal nodes for which antigenic advance should be calculated")
    parser.add_argument("--alignment", nargs="+", required=True, help="FASTA file(s) of sequences used to fit the given substitution weights with one sequence per node in the given tree")
    parser.add_argument("--gene-names", nargs="+", required=True, help="names of the genes represented by each alignment, given in the same order")
    parser.add_argument("--attribute-prefix", default="", help="prefix for node attributes in the JSON output including cumulative titer drop ('cTiterSub') and per-substitution titer drop ('dTiterSub'). Set a prefix to disambiguate annotations from multiple substitution model JSONs in the final Auspice JSON.")
    parser.add_argument("--group-by", help="name of a column to group substitution weights by")
    parser.add_argument("--output-node-data", required=True, help="node data JSON file with antigenic advance per node")
    parser.add_argument("--output-posterior", help="parquet file of antigenic advance per node and posterior sample")

    args = parser.parse_args()

    main(args)
