'''
Cluster viruses based on full genome sequence by identifying connected components.
Measures Hamming distance between pairs of full genomes and draws a connection if
Hamming distance is less than some cutoff.
'''

import argparse
import os
import json
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def hamming(array1, array2):
    return np.sum(array1 != array2)

def shared_strains(files):
    '''
    Return list of strains shared accross multiple JSONs
    '''
    setlist = []
    for fname in files:
        strains = set()
        if os.path.isfile(fname):
            with open(fname) as jfile:
                json_data = json.load(jfile)
                for node in json_data["nodes"]:
                    if node[0:5] != "NODE_":
                        strains.add(node)
        setlist.append(strains)
    shared = set.intersection(*setlist)
    return shared

def sequence_mapping(files, strains):
    '''
    Return dictionary mapping of strains to concatenated sequence across input files
    Assumes that all strains are in all files
    '''
    mapping = {}
    for strain in strains:
        mapping[strain] = np.array([])
    for fname in files:
        if os.path.isfile(fname):
            with open(fname) as jfile:
                json_data = json.load(jfile)
                for node, values in json_data["nodes"].items():
                    if node in mapping:
                        seq = values["sequence"]
                        seq = seq.replace("A", "0")
                        seq = seq.replace("T", "1")
                        seq = seq.replace("G", "2")
                        seq = seq.replace("C", "3")
                        array = np.asarray(list(seq), dtype = int)
                        mapping[node] = np.concatenate((mapping[node], array), axis=0).astype(int)
    return mapping

def strains_to_adjacency_matrix(strains, mapping, cutoff):
    '''
    Return np array of 0/1 for connected edges between all pairs of strains
    Connected edges are edges where Hamming distance is less than cutoff
    '''
    counter = 0
    interval = 100
    length = len(mapping)
    print("progress")
    adj_matrix = np.zeros((length, length))
    for indexA, strainA in enumerate(strains):
        if counter % interval == 0:
            print("[", end = '')
            for x in range(int(counter/interval)):
                print("-", end = '')
            for x in range(int(length/interval) - int(counter/interval)):
                print(" ", end = '')
            print("]")
        for indexB, strainB in enumerate(strains):
            distance = hamming(mapping[strainA], mapping[strainB])
            if distance < cutoff:
                adj_matrix[indexA, indexB] = 1
            else:
                adj_matrix[indexA, indexB] = 0
            if indexA == indexB:
                adj_matrix[indexA, indexB] = 0
        counter += 1
    return adj_matrix

def adjacency_matrix_to_connected_components(adj_matrix):
    '''
    Return connected components from adjacency matrix
    '''
    graph = csr_matrix(adj_matrix)
    return connected_components(csgraph=graph, directed=False, return_labels=True)

def components_to_cluster_json(strains, labels):
    '''
    Return JSON compatible data structure representing mapping of nodes to component labels
    "nodes": {
        "A/Hyogo/1061/2017": {
            "cluster": 1
        },
        "A/Peru/27/2015": {
            "cluster": 2
        },
        ...
    '''
    data = {}
    data["nodes"] = {}
    for strain, label in zip(strains, labels):
        data["nodes"][strain] = {"cluster": str(label)}
    return data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Cluster viruses based on full genome sequence by identifying connected components",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--nt-muts', nargs='+', type=str, required=True, help="list of nt-muts JSON files")
    parser.add_argument('--cutoff', default=5, type = int, help = "Hamming distance cutoff to be considered connected")
    parser.add_argument('--output', required=True, help="name of the file to write JSON data to")
    args = parser.parse_args()

    # collect strains shared across segments
    strains = shared_strains(args.nt_muts)

    # mapping of strains to concatenated sequence
    mapping = sequence_mapping(args.nt_muts, strains)

    # adjacency matrix via Hamming distance matrix
    adj_matrix = strains_to_adjacency_matrix(strains, mapping, args.cutoff)

    # connected components via adjacency matrix
    n_components, labels = adjacency_matrix_to_connected_components(adj_matrix)

    print("Defining", n_components, "clusters based on connected components")

    cluster_json = components_to_cluster_json(strains, labels)

    with open(args.output, 'wt') as fh:
        json.dump(cluster_json, fh, indent=1)
