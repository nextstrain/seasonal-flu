'''
Cluster viruses based on full genome sequence by identifying connected components.
Measures Hamming distance between pairs of full genomes and draws a connection if
Hamming distance is less than some cutoff.
'''

import argparse
import os
import json

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
        mapping[strain] = ""
    for fname in files:
        if os.path.isfile(fname):
            with open(fname) as jfile:
                json_data = json.load(jfile)
                for node, values in json_data["nodes"].items():
                    if node in mapping:
                        print(values["sequence"])
                        mapping[node] += values["sequence"]
    return(mapping)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Cluster viruses based on full genome sequence by identifying connected components",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--files', nargs='+', type=str, required=True, help="list of nt-muts JSON files")
    args = parser.parse_args()

    strains = shared_strains(args.files)
    print(strains)

    mapping = sequence_mapping(args.files, strains)
    print(mapping)
