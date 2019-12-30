import argparse, sys, os, glob, json
import numpy as np
from collections import defaultdict
from Bio import Phylo
from augur.utils import read_metadata, get_numerical_dates
from select_strains import read_strain_list, determine_time_interval, parse_metadata
from vaccination_coverage import read_all_vaccination_data

def calculate_average_on_tree(tree, func, min_clade_size=20):
    """Generic function that calculates clade averages on a tree

    Parameters
    ----------
    tree : Phylo.Clade.BaseTree
        phylogenetic tree
    func : callable
        a function that is called on terminal nodes and returns
        a tuple (value,1) when ther is data, otherwise (np.nan, 0)
    min_clade_size : int, optional
        smallest clade that is to get its own average. smaller clades
        inherit their parent averages

    Returns
    -------
    dict
        a dictionary with values for each node

    """
    for n in tree.find_clades(order='postorder'):
        if n.is_terminal():
            n.val, n.count = func(n)
        else:
            n.count = np.sum([c.count for c in n])
            n.val = np.sum([c.val for c in n if not np.isnan(c.val)])

    tree.root.val/=tree.root.count

    scores = {}
    scores[tree.root.name] = tree.root.val
    for n in tree.get_nonterminals(order='preorder'):
        for c in n:
            if c.count>min_clade_size:
                c.val /= c.count
            else:
                c.val = n.val
            scores[c.name] = c.val

    return scores


def calculate_average_age(tree, metadata, min_clade_size=20):
    """Function that assigns average host age to clades in the tree
    Parameters
    ----------
    tree : Phylo.Clade.BaseTree
        phylogenetic tree
    metadata : dict
        dictionary with meta data
    min_clade_size : int, optional
        smallest clade that is to get its own average. smaller clades
        inherit their parent averages

    Returns
    -------
    dict
        a dictionary with values for each node
    """
    def parse_age(n):
        if n.name in metadata and 'age' in metadata[n.name] and metadata[n.name]['age']!='unknown':
            return metadata[n.name]['age'], 1
        else:
            return np.nan, 0

    return calculate_average_on_tree(tree, parse_age, min_clade_size=min_clade_size)


def calc_average_vaccination_coverage(tree, metadata, min_clade_size=20):
    """calculate the average vaccination coverage for clades in the tree

    Parameters
    ----------
    tree : TYPE
        phylogenetic tree
    metadata : dict
        dictionary with metadata for each terminal node
    coverage : str
        file name of file that contains the vaccination coverage information
    min_clade_size : int, optional
        smallest clade that is to get its own average. smaller clades
        inherit their parent averages
    """
    vaccov = read_all_vaccination_data()
    def vaccination_coverage(n):
        if n.name in metadata and 'country' in metadata[n.name] \
            and metadata[n.name]['country'] in vaccov:
            return vaccov[metadata[n.name]['country']], 1
        else:
            return np.nan, 0

    scores = calculate_average_on_tree(tree, vaccination_coverage, min_clade_size=min_clade_size)
    if min_clade_size==0:
        for n in tree.get_terminals():
            if np.isnan(vaccination_coverage(n)[0]):
                scores[n.name] = 'unknown'

    return scores

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Annotate nodes with scores based on metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help="file with metadata associated with viral sequences, one for each segment")
    parser.add_argument('--tree', type=str, required=True, help="file inferred tree")
    parser.add_argument('--output', help="name of the file to write selected strains to")

    args = parser.parse_args()

    # read in meta data and tree
    metadata = parse_metadata(['s'], [args.metadata])['s']
    T = Phylo.read(args.tree, 'newick')

    # dictionary to hold calculated scores for terminal and internal nodes
    scores = dict()
    for k,v in calculate_average_age(T, metadata, min_clade_size=20).items():
        if k not in scores: scores[k] = {}
        scores[k]['avg_age'] = v

    for k,v in calc_average_vaccination_coverage(T, metadata, min_clade_size=0).items():
        if k not in scores: scores[k] = {}
        scores[k]['vaccov'] = v

    with open(args.output, 'w') as results:
        json.dump({"nodes":scores}, results, indent=1, sort_keys = True)
