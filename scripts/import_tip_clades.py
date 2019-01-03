"""
Take clades.json file that lists:
{
 "nodes": {
  "A/AbuDhabi/16/2017": {
   "clade_membership": "A1b/135N"
  },
...
and creates a new file that has internal nodes 'clade_membership' set to 'unassigned'.
"""

import argparse
import Bio
import Bio.Phylo
import json


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Import clade membership",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick tree originally used to assign tips to clades")
    parser.add_argument("--clades", required=True, help="JSON file with clade memberships")
    parser.add_argument("--output", required=True, help="JSON file with scrubbed clade memberships")
    args = parser.parse_args()

    tree = Bio.Phylo.read(args.tree, 'newick')

    with open(args.clades) as infile:
        json_data = json.load(infile)

    scrubbed_json_data = {'nodes':{}}

    # Copy clade membership for tips to a new JSON and omit internal nodes from
    # the original tree that was used to assign tips to clades.
    for node in tree.find_clades():
        if node.is_terminal():
            clade_membership = json_data['nodes'][node.name]['clade_membership']
            scrubbed_json_data['nodes'][node.name] = {'clade_membership': clade_membership}

    with open(args.output, 'w') as outfile:
        json.dump(scrubbed_json_data, outfile, indent=1, sort_keys=True)
