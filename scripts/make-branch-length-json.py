import argparse
from treetime.arg import assign_all_mccs, get_MCC_dict, get_mcc_map
from Bio import Phylo
import numpy as np
import pandas as pd
import json
import random
from os import path

def to_float(x):
    try:
        return float(x)
    except:
        return None

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="turn newick trees into branch length jsons",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--timetree', type = str, help='timetree for clock length.', required=True)
    parser.add_argument('--divtree', type = str, help='divtree for branch length.', required=True)
    parser.add_argument('--mccs', type = str, help='mccs.', required=True)
    parser.add_argument('--dates', type=str, required=True)
    parser.add_argument('--output-tree', type=str, required=True)
    parser.add_argument('--output-node-data', help="name of the file to write node data to", required=True)
    args = parser.parse_args()

    dates = pd.read_csv(args.dates, sep='\t', index_col=0, skiprows=1)

    timetree = Phylo.read(args.timetree, 'nexus')
    timetree.root.up = None
    for n in timetree.get_nonterminals():
        for c in n:
            c.up = n

    divtree = Phylo.read(args.divtree, 'nexus')

    MCCs_dict = get_MCC_dict(args.mccs)
    focal_tree_name = path.splitext(path.basename(args.timetree))[0].replace("treetime_", "")

    other_tree_names = []
    MCCs_list = [] ##list of all MCC pairs of focal tree with other trees
    for key in MCCs_dict.keys():
        if focal_tree_name in key:
            other_tree_names.append(list(key.difference(frozenset([focal_tree_name])))[0])
            MCCs_list.append(MCCs_dict[key])

    leaf_to_MCCs = get_mcc_map(MCCs_list, shuffle=True)

    assign_all_mccs(timetree, leaf_to_MCCs, 0)

    node_data = {}
    for node in timetree.find_clades():
        node_name = node.name or node.confidence
        numdate = to_float(dates.loc[node_name,"numeric date"])
        node.confidence = None
        node.name = node_name
        node.comment = ''
        if node.is_terminal() and (numdate is None):
            print("pruning node:", node_name, dates.loc[node_name,"numeric date"])

        node_data[node_name] = {"branch_length":node.branch_length,
                                "clock_length":node.branch_length,
                                "date": dates.loc[node_name,"date"],
                                "numdate": numdate,
                                "num_date_confidence":
                                    [to_float(dates.loc[node_name,"lower bound"]),
                                     to_float(dates.loc[node_name,"upper bound"])]
                                }
    for i in range(len(other_tree_names)):
        other_tree = other_tree_names[i]
        node_data[node_name]["mcc_"+other_tree] = node.mcc[i]

    for node in divtree.find_clades():
        node_name = node.name or node.confidence
        if node_name in node_data:
            node_data[node_name]["mutation_length"] = node.branch_length

    Phylo.write(timetree, args.output_tree, 'newick')

    with open(args.output_node_data, 'w') as fh:
        json.dump({"nodes": node_data}, fh)
