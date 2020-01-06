"""Calculate the distance between sequences between seasons.
"""
import argparse
from augur.frequency_estimators import TreeKdeFrequencies
from augur.utils import annotate_parents_for_tree, read_node_data, write_json
import Bio.Phylo
import json
import numpy as np


def get_titer_distance_between_nodes(tree, past_node, current_node, titer_attr="dTiter"):
    ancestor = tree.common_ancestor([past_node, current_node])
    nodes_between_tips = ancestor.get_path(past_node) + ancestor.get_path(current_node)
    return sum([node.attr[titer_attr] for node in nodes_between_tips])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--frequencies", help="frequencies JSON", required=True)
    parser.add_argument("--model-attribute-name", help="name of attribute to use from titer model file", default="dTiter")
    parser.add_argument("--attribute-name", help="name to store distances", required=True)
    parser.add_argument("--model", help="JSON providing the titer tree model", required=True)
    parser.add_argument("--date-annotations", help="JSON of branch lengths and date annotations from augur refine for samples in the given tree; required for comparisons to earliest or latest date", required=True)
    parser.add_argument("--years-back-to-compare", type=int, help="number of years prior to the current season to search for samples to calculate pairwise comparisons with", required=True)
    parser.add_argument("--max-past-samples", type=int, default=200, help="maximum number of past samples to randomly select for comparison to current samples")
    parser.add_argument("--min-frequency", type=float, default=0.0, help="minimum frequency to consider a sample alive")
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)

    args = parser.parse_args()

    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)

    frequencies = TreeKdeFrequencies.from_json(frequencies_json)
    pivots = frequencies.pivots

    # Identify pivots that belong within our search window for past samples.
    past_pivot_indices = (pivots < pivots[-1]) & (pivots >= pivots[-1] - args.years_back_to_compare)

    # Load date and titer model annotations and annotate tree with them.
    annotations = read_node_data([args.date_annotations, args.model])
    for node in tree.find_clades():
        node.attr = annotations["nodes"][node.name]
        node.attr["num_date"] = node.attr["numdate"]

    # Identify samples to compare including those in the current timepoint
    # (pivot) and those in previous timepoints.
    current_samples = []
    past_samples = []
    date_by_sample = {}
    tips_by_sample = {}
    for tip in tree.find_clades(terminal=True):
        # Samples with nonzero frequencies in the last timepoint are current
        # samples. Those with one or more nonzero frequencies in the search
        # window of the past timepoints are past samples.
        if frequencies.frequencies[tip.name][-1] > args.min_frequency:
            current_samples.append(tip.name)
            tips_by_sample[tip.name] = tip
        elif (frequencies.frequencies[tip.name][past_pivot_indices] > args.min_frequency).sum() > 0:
            past_samples.append(tip.name)
            tips_by_sample[tip.name] = tip

        date_by_sample[tip.name] = tip.attr["numdate"]

    print("Expecting %i comparisons for %i current and %i past samples" % (len(current_samples) * args.max_past_samples, len(current_samples), len(past_samples)))
    distances_by_node = {}
    comparisons = 0

    for current_sample in current_samples:
        if not current_sample in distances_by_node:
            distances_by_node[current_sample] = {}

        if not args.attribute_name in distances_by_node[current_sample]:
            distances_by_node[current_sample][args.attribute_name] = {}

        for past_sample in np.random.choice(past_samples, size=min(args.max_past_samples, len(past_samples)), replace=False):
            # The past is in the past.
            if date_by_sample[past_sample] < date_by_sample[current_sample]:
                distances_by_node[current_sample][args.attribute_name][past_sample] = get_titer_distance_between_nodes(
                    tree,
                    tips_by_sample[past_sample],
                    tips_by_sample[current_sample],
                    args.model_attribute_name
                )

            comparisons += 1
            if comparisons % 1000 == 0:
                print("Completed", comparisons, "comparisons, with last distance of", distances_by_node[current_sample][args.attribute_name][past_sample], flush=True)

    print("Calculated %i comparisons" % comparisons)
    # Prepare params for export.
    params = {
        "attribute": args.attribute_name,
        "years_back_to_compare": args.years_back_to_compare
    }

    # Export distances to JSON.
    write_json({"params": params, "nodes": distances_by_node}, args.output)
