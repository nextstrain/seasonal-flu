#!/usr/bin/env python3
import argparse

from augur.io import write_json
from augur.utils import annotate_parents_for_tree, read_node_data, read_tree
from bisect import bisect_right
import pandas as pd


def annotate_tree(mutations_by_node, tree, weights, starts, ends, dates):
    """Annotates log convergence ratio attributes to nodes of a given tree.
    """
    lcr_sum = {tree.root.name: 0.0}
    total = 0
    missing = 0
    outside_windows = 0
    after_last_window = 0

    for node in tree.find_clades():
        # Get mutations between current node and its parent.
        mutations = [
            (gene, mutation)
            for gene, gene_mutations in mutations_by_node[node.name].items()
            for mutation in gene_mutations
        ]

        # Find the window to use for scoring the current node (earliest end date after the branch date)
        date = dates[node.name]
        window = bisect_right(ends, date)

        # Calculate the log convergence ratio sum on the branch to the current node.
        total += len(mutations)
        branch = 0.0

        if window == len(ends):
            # Past every window's end, so score with the last window.
            after_last_window += len(mutations)
            window = len(ends) - 1

        if starts[window] > date:
            outside_windows += len(mutations)
        else:
            for gene, mutation in mutations:
                weight = weights.get((f"{gene}:{mutation}", window))
                if weight is None:
                    missing += 1
                    continue

                branch += weight

        # Calculate cumulative log convergence ratio from root to current node.
        lcr_sum[node.name] = branch
        if getattr(node, "parent"):
            lcr_sum[node.name] += lcr_sum.get(node.parent.name, 0.0)


    return lcr_sum, total, missing, outside_windows, after_last_window


def main(args):
    """Calculate log convergence ratio.
    """
    # Load log convergence ratios &  time windows
    weight_table = pd.read_csv(args.weights, sep="\t", comment="#")

    # Sort windows by end date to help find first window with end date after current node in annotate_tree
    windows = pd.read_csv(args.windows, sep="\t").sort_values("window_end_numeric", ignore_index=True)
    starts = windows["window_start_numeric"].tolist()
    ends = windows["window_end_numeric"].tolist()
    window_index = {end: index for index, end in enumerate(windows["window_end"])}

    weights = {
        (f"{row.gene}:{row.from_aa}{row.site}{row.derived_aa}", window_index[row.window_end]): float(row.lcr)
        for row in weight_table.itertuples()
    }

    # Load node dates.
    dates = {
        name: node["numdate"]
        for name, node in read_node_data(args.branch_lengths)["nodes"].items()
    }

    # Load tree.
    tree = read_tree(args.tree)
    tree = annotate_parents_for_tree(tree)

    # Load amino acid mutations per gene on the branch to each node.
    mutations_by_node = {
        name: node["aa_muts"]
        for name, node in read_node_data(args.mutations)["nodes"].items()
    }

    lcr_sum, total, missing, outside_windows, after_last_window = annotate_tree(
        mutations_by_node,
        tree,
        weights,
        starts,
        ends,
        dates,
    )

    node_data = {
        node.name: {"lcr_sum": lcr_sum[node.name]}
        for node in tree.find_clades()
    }

    # Export log convergence ratio sum per node.
    write_json(
        {
            "nodes": node_data,
            "params": {
                "mutations": total,
                "mutations_missing_from_weights": missing,
                "mutations_outside_windows": outside_windows,
                "mutations_after_last_window": after_last_window,
            },
        },
        args.output_node_data
    )

    # print to log
    print(f"{total} mutations, {missing} missing from the weights table, {outside_windows} in no window, {after_last_window} after the last window")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--weights", required=True, help="TSV file with one log convergence ratio per substitution per time window, defined in columns 'gene', 'site', 'from_aa', 'derived_aa', 'window_end', and 'lcr'. 'window_end' is the last day included in the time window, and is used to identify the window in the windows file")
    parser.add_argument("--windows", required=True, help="TSV file with the windows used in the weights file, containing 'window_start', 'window_end', 'window_start_numeric', and 'window_end_numeric'")
    parser.add_argument("--tree", required=True, help="Newick file of a tree with named internal nodes for which log convergence ratio should be calculated")
    parser.add_argument("--branch-lengths", required=True, help="node data JSON from augur refine with a numeric date ('numdate') per node in the given tree")
    parser.add_argument("--mutations", required=True, help="node data JSON from augur ancestral with amino acid mutations ('aa_muts') per node in the given tree")
    parser.add_argument("--output-node-data", required=True, help="node data JSON file with log convergence ratio per node")

    args = parser.parse_args()

    main(args)
