"""
Author: John Huddleston <huddlej at gmail.com>
"""
import argparse
import Bio.Phylo
import json
import pandas as pd
import sys


def annotate_parents_for_tree(tree):
    """Annotate each node in the given tree with its parent.

    Examples
    --------
    >>> import io
    >>> tree = Bio.Phylo.read(io.StringIO("(A, (B, C))"), "newick")
    >>> not any([hasattr(node, "parent") for node in tree.find_clades()])
    True
    >>> tree = annotate_parents_for_tree(tree)
    >>> tree.root.parent is None
    True
    >>> all([hasattr(node, "parent") for node in tree.find_clades()])
    True
    """
    tree.root.parent = None
    for node in tree.find_clades(order="level"):
        for child in node.clades:
            child.parent = node

    # Return the tree.
    return tree


def json_to_tree(json_dict, root=True, parent_cumulative_branch_length=None):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.

    Assigns links back to parent nodes for the root of the tree.

    Test opening a JSON from augur export v1.

    >>> import json
    >>> json_fh = open("tests/data/json_tree_to_nexus/flu_h3n2_ha_3y_tree.json", "r")
    >>> json_dict = json.load(json_fh)
    >>> tree = json_to_tree(json_dict)
    >>> tree.name
    'NODE_0002020'
    >>> len(tree.clades)
    2
    >>> tree.clades[0].name
    'NODE_0001489'
    >>> hasattr(tree, "attr")
    True
    >>> "dTiter" in tree.attr
    True
    >>> tree.clades[0].parent.name
    'NODE_0002020'
    >>> tree.clades[0].branch_length > 0
    True

    Test opening a JSON from augur export v2.

    >>> json_fh = open("tests/data/zika.json", "r")
    >>> json_dict = json.load(json_fh)
    >>> tree = json_to_tree(json_dict)
    >>> hasattr(tree, "name")
    True
    >>> len(tree.clades) > 0
    True
    >>> tree.clades[0].branch_length > 0
    True

    Branch lengths should be the length of the branch to each node and not the
    length from the root. The cumulative branch length from the root gets its
    own attribute.

    >>> tip = [tip for tip in tree.find_clades(terminal=True) if tip.name == "USA/2016/FLWB042"][0]
    >>> round(tip.cumulative_branch_length, 6)
    0.004747
    >>> round(tip.branch_length, 6)
    0.000186

    """
    # Check for v2 JSON which has combined metadata and tree data.
    if root and "meta" in json_dict and "tree" in json_dict:
        json_dict = json_dict["tree"]

    node = Bio.Phylo.Newick.Clade()

    # v1 and v2 JSONs use different keys for strain names.
    if "name" in json_dict:
        node.name = json_dict["name"]
    else:
        node.name = json_dict["strain"]

    # Assign all non-children attributes.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    # Only v1 JSONs support a single `attr` attribute.
    if hasattr(node, "attr"):
        node.numdate = node.attr.get("num_date")
        node.cumulative_branch_length = node.attr.get("div")

        if "translations" in node.attr:
            node.translations = node.attr["translations"]
    elif hasattr(node, "node_attrs"):
        node.cumulative_branch_length = node.node_attrs.get("div")

    node.branch_length = 0.0
    if parent_cumulative_branch_length is not None and hasattr(node, "cumulative_branch_length"):
        node.branch_length = node.cumulative_branch_length - parent_cumulative_branch_length

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [
            json_to_tree(
                child,
                root=False,
                parent_cumulative_branch_length=node.cumulative_branch_length
            )
            for child in json_dict["children"]
        ]

    if root:
        node = annotate_parents_for_tree(node)

    return node


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, help="auspice tree JSON")
    parser.add_argument("--output-metadata", help="tab-delimited file of attributes per node of the given tree")
    parser.add_argument("--output-tree", help="Newick version of the given tree")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data from internal nodes in metadata output")
    parser.add_argument("--attributes", nargs="+", help="names of attributes to export from the given tree in the metadata output")

    args = parser.parse_args()

    # Load tree from JSON.
    with open(args.tree, "r", encoding="utf-8") as fh:
        tree_json = json.load(fh)

    tree = json_to_tree(tree_json)

    # Output the tree, if requested.
    if args.output_tree:
        Bio.Phylo.write(
            tree,
            args.output_tree,
            "newick",
        )

    if args.output_metadata:
        # Collect attributes per node from the tree to export.
        records = []

        if args.attributes:
            attributes = args.attributes
        else:
            attributes = sorted(
                set(tree.root.node_attrs.keys()) |
                set(tree.root.branch_attrs.keys())
            )

        for node in tree.find_clades():
            if node.is_terminal() or args.include_internal_nodes:
                record = {
                    "name": node.name
                }

                for attribute in attributes:
                    if attribute in node.node_attrs:
                        value = node.node_attrs[attribute]
                    elif attribute in node.branch_attrs:
                        value = node.branch_attrs[attribute]
                    else:
                        print(f"Could not find attribute '{attribute}' for node '{node.name}'.", file=sys.stderr)
                        value = None

                    if value is not None:
                        if isinstance(value, dict) and "value" in value:
                            value = value["value"]

                    record[attribute] = value

                record["parent_name"] = getattr(node.parent, "name", "")
                records.append(record)

        # Convert records to a data frame and save as a tab-delimited file.
        df = pd.DataFrame(records)
        df.to_csv(
            args.output_metadata,
            sep="\t",
            header=True,
            index=False,
            columns=["name", "parent_name"] + list(attributes),
        )
