"""Calculate the distance between sequences between seasons.
"""
import argparse
from augur.frequency_estimators import float_to_datestring, timestamp_to_float
from augur.utils import annotate_parents_for_tree, read_node_data, write_json
import Bio.Phylo
import json
import numpy as np
import pandas as pd


def get_titer_distance_between_nodes(tree, past_node, current_node, titer_attr="dTiter"):
    # Find MRCA of tips from one tip up. Sum the titer attribute of interest
    # while walking up to the MRCA, to avoid an additional pass later. The loop
    # below stops when the past node is found in the list of the candidate
    # MRCA's terminals. This test should always evaluate to true when the MRCA
    # is the root node, so we should not have to worry about trying to find the
    # parent of the root.
    current_node_branch_sum = 0.0
    mrca = current_node
    while past_node.name not in mrca.terminals:
        current_node_branch_sum += mrca.attr[titer_attr]
        mrca = mrca.parent

    # Sum the node weights for the other tip from the bottom up until we reach
    # the MRCA. The value of the MRCA is intentionally excluded here, as it
    # would represent the branch leading to the MRCA and would be outside the
    # path between the two tips.
    past_node_branch_sum = 0.0
    current_node = past_node
    while current_node != mrca:
        past_node_branch_sum += current_node.attr[titer_attr]
        current_node = current_node.parent

    final_sum = past_node_branch_sum + current_node_branch_sum
    return final_sum


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--frequencies", help="frequencies JSON", required=True)
    parser.add_argument("--model-attribute-name", help="name of attribute to use from titer model file", default="dTiter")
    parser.add_argument("--attribute-name", help="name to store distances", required=True)
    parser.add_argument("--model", help="JSON providing the titer tree model", required=True)
    parser.add_argument("--date-annotations", help="JSON of branch lengths and date annotations from augur refine for samples in the given tree; required for comparisons to earliest or latest date", required=True)
    parser.add_argument("--months-back-for-current-samples", type=int, help="number of months prior to the last date with estimated frequencies to include samples as current", required=True)
    parser.add_argument("--years-back-to-compare", type=int, help="number of years prior to the current season to search for samples to calculate pairwise comparisons with", required=True)
    parser.add_argument("--max-past-samples", type=int, default=200, help="maximum number of past samples to randomly select for comparison to current samples")
    parser.add_argument("--min-frequency", type=float, default=0.0, help="minimum frequency to consider a sample alive")
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)

    args = parser.parse_args()

    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Make a single pass through the tree in postorder to store a set of all
    # terminals descending from each node. This uses more memory, but it allows
    # faster identification of MRCAs between any pair of tips in the tree and
    # speeds up pairwise distance calculations by orders of magnitude.
    for node in tree.find_clades(order="postorder"):
        node.terminals = set()
        for child in node.clades:
            if child.is_terminal():
                node.terminals.add(child.name)
            else:
                node.terminals.update(child.terminals)

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies = json.load(fh)

    pivots = np.array(frequencies.pop("pivots"))

    # Identify pivots that belong within our search window for past samples.
    # First, calculate dates associated with the interval for current samples
    # based on the number of months back requested. Then, calculate interval for
    # past samples with an upper bound based on the earliest current samples and
    # a lower bound based on the years back requested.
    last_pivot_datetime = pd.to_datetime(float_to_datestring(pivots[-1]))
    last_current_datetime = last_pivot_datetime - pd.DateOffset(months=args.months_back_for_current_samples)
    last_past_datetime = last_pivot_datetime - pd.DateOffset(years=args.years_back_to_compare)

    # Find the pivot indices that correspond to the current and past pivots.
    current_pivot_indices = np.array([
        pd.to_datetime(float_to_datestring(pivot)) > last_current_datetime
        for pivot in pivots
    ])
    past_pivot_indices = np.array([
        ((pd.to_datetime(float_to_datestring(pivot)) >= last_past_datetime) &
         (pd.to_datetime(float_to_datestring(pivot)) <= last_current_datetime))
        for pivot in pivots
    ])

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
        frequencies[tip.name]["frequencies"] = np.array(frequencies[tip.name]["frequencies"])
        if (frequencies[tip.name]["frequencies"][current_pivot_indices] > args.min_frequency).sum() > 0:
            current_samples.append(tip.name)
            tips_by_sample[tip.name] = tip
        elif (frequencies[tip.name]["frequencies"][past_pivot_indices] > args.min_frequency).sum() > 0:
            past_samples.append(tip.name)
            tips_by_sample[tip.name] = tip

        date_by_sample[tip.name] = tip.attr["numdate"]

    print("Expecting %i comparisons for %i current and %i past samples" % (len(current_samples) * len(past_samples), len(current_samples), len(past_samples)))
    distances_by_node = {}
    comparisons = 0

    for current_sample in current_samples:
        if not current_sample in distances_by_node:
            distances_by_node[current_sample] = {}

        if not args.attribute_name in distances_by_node[current_sample]:
            distances_by_node[current_sample][args.attribute_name] = {}

        for past_sample in past_samples:
            # The past is in the past.
            if date_by_sample[past_sample] < date_by_sample[current_sample]:
                distances_by_node[current_sample][args.attribute_name][past_sample] = np.around(get_titer_distance_between_nodes(
                    tree,
                    tips_by_sample[past_sample],
                    tips_by_sample[current_sample],
                    args.model_attribute_name
                ), 4)

                comparisons += 1
                if comparisons % 10000 == 0:
                    print("Completed", comparisons, "comparisons, with last distance of", distances_by_node[current_sample][args.attribute_name][past_sample], flush=True)

    print("Calculated %i comparisons" % comparisons)
    # Prepare params for export.
    params = {
        "attribute": args.attribute_name,
        "years_back_to_compare": args.years_back_to_compare
    }

    # Export distances to JSON.
    write_json({"params": params, "nodes": distances_by_node}, args.output, indent=None)
