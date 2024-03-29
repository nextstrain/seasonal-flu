import argparse
from treetime.arg import assign_mccs
from Bio import Phylo
import numpy as np
import pandas as pd
import json
import random

def to_float(x):
    try:
        return float(x)
    except:
        return None

def get_clock(fname):
    with open(fname) as fh:
        for line in fh:
            if 'rate:' in line:
                rate = float(line.strip().split(':')[-1])

    return rate

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="turn newick trees into branch length jsons",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--timetree', type = str, help='timetree for clock length.', required=True)
    parser.add_argument('--divtree', type = str, help='divtree for branch length.', required=True)
    parser.add_argument('--molecular-clock', type = str, help='file with clock rate.', required=True)    
    parser.add_argument('--mccs', type = str, help='mccs.', required=True)
    parser.add_argument('--dates', type=str, required=True)
    parser.add_argument('--output-tree', type=str, required=True)
    parser.add_argument('--output-node-data', help="name of the file to write node data to", required=True)
    args = parser.parse_args()

    dates = pd.read_csv(args.dates, sep='\t', index_col=0, skiprows=1)
    rate = get_clock(args.molecular_clock)

    timetree = Phylo.read(args.timetree, 'nexus')
    timetree.root.up = None
    for n in timetree.get_nonterminals():
        for c in n:
            c.branch_length *=  rate
            c.up = n

    divtree = Phylo.read(args.divtree, 'nexus')

    MCCs = []
    n_large_mccs = 0
    with open(args.mccs) as fh:
        for line in fh:
            if line.strip():
                MCCs.append(line.strip().split(','))
                if len(MCCs[-1])>2:
                    n_large_mccs += 1


    MCCs.sort(key=lambda x:len(x), reverse=True)
    mcc_map = list(range(len(MCCs)))
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    random.seed(987)
    #random.shuffle(mcc_map)

    leaf_to_MCC = {}
    for mi,mcc in enumerate(MCCs):
        if mi<=n_large_mccs:
            label = letters[mi%len(letters)]
            if n_large_mccs>len(letters):
                label = letters[mi//len(letters)] + label
        else:
            label = '--'
        for leaf in mcc:
            leaf_to_MCC[leaf] = label

    assign_mccs(timetree, leaf_to_MCC, 0)

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
                                     to_float(dates.loc[node_name,"upper bound"])],
                                "mcc": node.mcc
                                }


    for node in divtree.find_clades():
        node_name = node.name or node.confidence
        if node_name in node_data:
            node_data[node_name]["mutation_length"] = node.branch_length

    Phylo.write(timetree, args.output_tree, 'newick')

    with open(args.output_node_data, 'w') as fh:
        json.dump({"nodes": node_data}, fh)
