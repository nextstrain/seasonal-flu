#!/usr/bin/env python3
import argparse
import json


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--titers-sub")
    parser.add_argument("--titers-tree")
    parser.add_argument("--output-titers")
    parser.add_argument("--output-titers-sub")
    parser.add_argument("--output-titers-tree")

    args = parser.parse_args()

    with open(args.titers_sub) as fh:
        sub = json.load(fh)

    with open(args.output_titers_sub, 'wt') as sub_file:
        json.dump({'avidity': sub['avidity'],
                    'potency': sub['potency'],
                    'substitution': sub['substitution']},
                    sub_file, indent=1)

    with open(args.output_titers, 'wt') as raw_file:
        json.dump(sub['titers'], raw_file, indent=1)

    with open(args.titers_tree) as fh:
        tree = json.load(fh)

    with open(args.output_titers_tree, 'wt') as tree_file:
        json.dump({'avidity': tree['avidity'],
                    'potency': tree['potency'],
                    'dTiter': {k:v['dTiter'] for k,v in tree['nodes'].items()}},
                    tree_file, indent=1)
