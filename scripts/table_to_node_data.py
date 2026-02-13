"""Create Augur-compatible node data JSON from a pandas data frame.
"""
import argparse
import sys

import pandas as pd

from augur.utils import annotate_parents_for_tree, read_tree, write_json


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, help="table to convert to a node data JSON")
    parser.add_argument("--tree", help="tree with named internal nodes that match the index column values in the given table. Only required for assigning branch labels.")
    parser.add_argument("--index-column", default="strain", help="name of the column to use as an index")
    parser.add_argument("--columns", nargs="+", help="optional list of columns in the given table to include in the output. If not provided, all columns will be included.")
    parser.add_argument("--node-name", default="nodes", help="name of the node data attribute in the JSON output")
    parser.add_argument("--branch-labels", nargs="+", help="optional map of column names to branch labels. You must supply this for each column you would like to label. By default the branch label key the same as the column name, but you may customise this via the COLUMN=LABEL syntax.")
    parser.add_argument("--column-to-node-attribute", nargs="+", help="optional map of column names to node attributes in the node data JSON in the format of COLUMN=ATTRIBUTE")
    parser.add_argument("--error-on-polyphyletic-labels", action="store_true", help="exit with an error if the script finds any polyphyletic branch labels")
    parser.add_argument("--output", required=True, help="node data JSON file")

    args = parser.parse_args()

    delimiter = "," if args.table.endswith(".csv") else "\t"

    table = pd.read_csv(
        args.table,
        sep=delimiter,
        index_col=args.index_column,
        dtype=str,
    )

    if args.columns:
        table = table.loc[:, args.columns].copy()

    # Optionally rename columns to new node attribute names.
    column_to_node_attribute = {}
    if args.column_to_node_attribute:
        column_to_node_attribute = dict(
            tuple(value.split("="))
            for value in args.column_to_node_attribute
        )

    table_dict = table.rename(columns=column_to_node_attribute).transpose().to_dict()
    node_data = {
        args.node_name: table_dict,
    }

    # Optionally annotate branch labels for internal nodes, if a tree is given.
    found_polyphyletic_labels = False
    if args.branch_labels:
        if not args.tree:
            print(
                "ERROR: You must provide a Newick tree with named internal nodes (e.g., from augur refine) to assign branch labels.",
                file=sys.stderr,
            )
            sys.exit(1)

        # Load the tree.
        tree = annotate_parents_for_tree(read_tree(args.tree))

        # Parse branch columns and labels.
        branch_label_by_column = {}
        for branch_label in args.branch_labels:
            if "=" in branch_label:
                column, label = branch_label.split("=")
            else:
                column = label = branch_label

            branch_label_by_column[column] = label

        # For each branch label column, find all distinct values of the column
        # and then find the first node in the tree with each value.
        branches = {}
        for column, label in branch_label_by_column.items():
            # Get distinct values for this column.
            branch_values = set(table[column].drop_duplicates().values)
            labeled_branch_values = set()

            # Map each node to its value.
            value_by_node = dict(table[column].reset_index().values)

            # Using a preorder traversal, find the first node in the tree with
            # each distinct value.
            for node in tree.find_clades():
                node_value = value_by_node.get(node.name)
                parent_node_value = value_by_node.get(getattr(node.parent, "name", None))

                # Check whether the current node value is among the values we
                # want to label.
                if node_value in branch_values and node_value != parent_node_value:
                    # Check whether the current value has already been assigned
                    # to a branch. If it hasn't, assign the label to the branch.
                    # If it has been assigned and the current node represents
                    # another clade with the value, let the user know that the
                    # value could not be assigned monophyletically.
                    if node.name not in branches:
                        branches[node.name] = {"labels": {}}

                    if node_value not in labeled_branch_values:
                        branches[node.name]["labels"][label] = node_value
                        labeled_branch_values.add(node_value)
                    else:
                        print(
                            f"WARNING: The value {node_value} for the column {column} is not monophyletic ({node.name})",
                            file=sys.stderr,
                        )

                    branches[node.name]["labels"][label + "_by_node"] = node_value + ":" + node.name
                    # if node_value not in labeled_branch_values:
                    #     if node.name not in branches:
                    #         branches[node.name] = {"labels": {}}

                    #     branches[node.name]["labels"][label] = node_value
                    #     labeled_branch_values.add(node_value)
                    # elif node_value != parent_node_value:
                    #     found_polyphyletic_labels = True
                    #     print(
                    #         f"WARNING: The value {node_value} for the column {column} is not monophyletic ({node.name})",
                    #         file=sys.stderr,
                    #     )

        node_data["branches"] = branches

    if found_polyphyletic_labels and args.error_on_polyphyletic_labels:
        sys.exit(1)

    write_json(node_data, args.output)
