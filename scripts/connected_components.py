'''
Cluster viruses based on full genome sequence by identifying connected components.
Measures Hamming distance between pairs of full genomes and draws a connection if
Hamming distance is less than some cutoff.
'''

import argparse
import os
import json
import numpy as np

hamming_lookup = {}
def hamming(strain1, strain2, array1, array2):
    if (strain1, strain2) in hamming_lookup:
        return hamming_lookup[(strain1, strain2)]
    dist = np.sum(array1 != array2)
    hamming_lookup[(strain1, strain2)] = dist
    return dist

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
    Return dictionary mapping of strains to contenated sequence across input files
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
                        mapping[node] = np.concatenate((mapping[node], np.array(list(seq))), axis=0)
    return mapping

def attempt_merge(clusters, mapping, cutoff):
    for clusterA in clusters:
        for clusterB in clusters:
            if clusterA != clusterB:
                # compare all strains in clusterA to all strains in clusterB
                for strainA in clusters[clusterA]:
                    for strainB in clusters[clusterB]:
                        if strainA != strainB:
                            distance = hamming(strainA, strainB, mapping[strainA], mapping[strainB])
                            if distance <= cutoff:
                                clusters[clusterA] = [x for x in clusters[clusterA]] + [x for x in clusters[clusterB]]
                                clusters.pop(clusterB, None)
                                return True
    return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Cluster viruses based on full genome sequence by identifying connected components",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--files', nargs='+', type=str, required=True, help="list of nt-muts JSON files")
    parser.add_argument('--cutoff', default=5, type = int,  help = "Hamming distance cutoff to be considered connected")
    args = parser.parse_args()

    # collect strains shared across segments
    strains = shared_strains(args.files)

    # mapping of strains to concatenated sequence
    mapping = sequence_mapping(args.files, strains)

    # start with each strain in its own cluster
    # data structure is dict of cluster id -> list of strains
    clusters = {}
    for index, strain in enumerate(strains):
        clusters[index] = [strain]

    # iterate over clusters and attempt to merge
    test = True
    while test:
        test = attempt_merge(clusters, mapping, args.cutoff)

    print(clusters)
