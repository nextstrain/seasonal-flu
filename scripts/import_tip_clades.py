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
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Import clade membership",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--clades", required=True, help="JSON file with clade memberships")
    parser.add_argument("--output", required=True, help="JSON file with scrubbed clade memberships")
    args = parser.parse_args()

    with open(args.clades) as infile:
        json_data = json.load(infile)

    scrubbed_json_data = {'nodes':{}}

    for node, values in json_data['nodes'].items():
        clade_membership = values['clade_membership']
        if node[0:4] == 'NODE':
            clade_membership = 'unassigned'
        scrubbed_json_data['nodes'][node] = {'clade_membership': clade_membership}

    with open(args.output, 'w') as outfile:
        json.dump(scrubbed_json_data, outfile, indent=1, sort_keys=True)
