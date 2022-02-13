import argparse
from Bio import Phylo
import numpy as np
import pandas as pd
import json

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

    divtree = Phylo.read(args.divtree, 'nexus')

    MCCs = []
    with open(args.mccs) as fh:
        for line in fh:
            if line.strip():
                MCCs.append(line.strip().split(','))

    leaf_to_MCC = {}
    for mi,mcc in enumerate(MCCs):
        for leaf in mcc:
            leaf_to_MCC[leaf] = mi


    node_data = {}
    for node in timetree.find_clades():
        numdate = to_float(dates.loc[node_name,"numeric date"])
        if node.is_terminals() and (numdate is None):
            timetree.prune(node)
            continue
        node_name = node.name or node.confidence
        node.confidence = None
        node.name = node_name
        node.comment = ''
        node_data[node_name] = {"branch_length":node.branch_length,
                                "clock_length":node.branch_length,
                                "date": dates.loc[node_name,"date"],
                                "numdate": numdate,
                                "num_date_confidence":
                                    [to_float(dates.loc[node_name,"lower bound"]),
                                     to_float(dates.loc[node_name,"upper bound"])],
                                "mcc": leaf_to_MCC.get(node_name, None)
                                }


    for node in divtree.find_clades():
        node_name = node.name or node.confidence
        if node_name in node_data:
            node_data[node_name]["mutation_length"] = node.branch_length

    Phylo.write(timetree, args.output_tree, 'newick')

    with open(args.output_node_data, 'w') as fh:
        json.dump({"nodes": node_data}, fh)