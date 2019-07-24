import argparse
from Bio import Phylo
import numpy as np
from treetime.utils import numeric_date
from augur.utils import read_metadata, get_numerical_dates, read_node_data, write_json
from select_strains import parse_metadata

naming_alphabets = [list(map(str,range(1,10))),
                     list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'),
                     list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'.lower())]

separator = '.'

def next_clade(previous, alphabet):
    max_digits = len(previous)+1
    n = len(alphabet)
    digit_representation = [alphabet.index(x) for x in previous]
    digit_representation.reverse()
    c = int(np.sum([b*n**i for i,b in enumerate(digit_representation)]))+1
    new_digits = [alphabet[(c//(n**i))%n] for i in range(max_digits) if c//(n**i)]
    new_digits.reverse()
    return ''.join(new_digits)


def generate_clade_name(node, existing_names):
    '''
    to generate a name, we need to know how many sister clades at the same level
    exist and what the next 'digit' at this level will be.
    '''
    super_clade = node.clade
    prefix = super_clade+separator
    sister_clades = sorted([clade[len(prefix):].split(separator)[0] for clade in existing_names if clade.startswith(prefix)])

    # alphabet are simply cycled.
    alphabet = naming_alphabets[len(super_clade.split(separator))%len(naming_alphabets)]

    if len(sister_clades):
        print(sister_clades)
        new_clade = super_clade + separator + next_clade(sister_clades[-1], alphabet)
    else:
        new_clade = super_clade + separator + alphabet[0]
    print(super_clade, new_clade)
    return new_clade


def assign_new_clade_name(node, new_name):
    '''
    assign a new clade name to a node and reset names and distance of all descendent nodes
    '''
    for n in node.find_clades(order='preorder'):
        if n==node:
            node.clade = new_name
            node.distance_from_last_clade = 0
        else:
            n.clade = n.up.clade
            n.distance_from_last_clade = n.up.distance_from_last_clade + len(n.mutations)


def name_new_clades(tree, tree_frequency_index, frequency_threshold=0.1, distance_threshold=3):
    '''
    take stock of existing clades and then annotate new clades if they fulfill
    the frequencyd/divergence criteria in the time interval specified in tree_frequency_index
    '''
    assert hasattr(tree.root, 'clade')
    print("assigning clades with distance threshold {} and frequency threshold {}".format(distance_threshold, frequency_threshold))
    existing_clades = set()
    for n in tree.find_clades(order='preorder'):
        if n.up is None:
            existing_clades.add(n.clade)
            n.distance_from_last_clade = 0
            continue

        if not hasattr(n, "clade"):
            n.clade = n.up.clade

        if n.clade==n.up.clade:
            n.distance_from_last_clade = n.up.distance_from_last_clade + len(n.mutations)
        else:
            existing_clades.add(n.clade)
            n.distance_from_last_clade = 0

    for n in tree.find_clades(order='preorder'):
        if n.freq[tree_frequency_index]>frequency_threshold and \
           n.distance_from_last_clade>distance_threshold:
            new_name = generate_clade_name(n, existing_clades)
            assign_new_clade_name(n, new_name)
            existing_clades.add(new_name)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign clade names",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True,
                        help="newick file with the tree")
    parser.add_argument('--metadata', nargs='+', help="file with metadata associated with viral sequences, one for each segment")
    parser.add_argument('--mutations', nargs='+', help='JSON(s) containing ancestral and tip nucleotide and/or amino-acid mutations ')
    parser.add_argument('--frequency-threshold', type=float, default=0.2, help='minimal population fraction to name clade')
    parser.add_argument('--distance-threshold', type=int, default=3, help='minimal distance (number of mutations) from parent clade to name clade')
    parser.add_argument('--date-format', type=str, default="%Y-%m-%d", help="date format")
    parser.add_argument('--output', help="name of the file to write selected strains to")


    args = parser.parse_args()

    T = Phylo.read(args.tree, 'newick')
    # read in meta data, parse numeric dates
    metadata = parse_metadata(['ha'], args.metadata, date_format=args.date_format)['ha']
    node_data = read_node_data(args.mutations, args.tree)

    all_dates = []
    for n in T.get_terminals():
        if n.name in metadata:
            n.num_date = metadata[n.name]['num_date']
            all_dates.append(n.num_date)
        else:
            n.num_date = None

    # generate time bins to calculate frequencies in
    all_dates.sort()
    date_range = (all_dates[int(0.05*len(all_dates))], all_dates[-1])
    dt = 0.5
    date_bins = np.arange(date_range[0], date_range[1]+dt,dt)

    # count number of tips of each node in each time window => frequency
    for n in T.find_clades(order='postorder'):
        if n.name in node_data["nodes"]:
            n.mutations = node_data["nodes"][n.name]['muts']
        if n.is_terminal():
            n.bin_count = np.zeros_like(date_bins)
            n.bin_count[min(len(date_bins)-1, date_bins.searchsorted(n.num_date))] += 1
        else:
            n.bin_count = np.sum([c.bin_count for c in n], axis=0)

    # convert counts in time windows into frequencies and prep tree
    T.root.up = None
    for n in T.find_clades(order='preorder'):
        n.freq = n.bin_count/T.root.bin_count
        if not n.is_terminal():
            for c in n:
                c.up = n

    # name clades in the tree by going through all time slices consecutively
    T.root.clade = '1'
    for i,threshold in enumerate(date_bins):
        name_new_clades(T, i, frequency_threshold=args.frequency_threshold,
                        distance_threshold=args.distance_threshold)

    # collect clade names and generate augur clade json
    clades = {}
    for n in T.find_clades():
        if n.up:
            if n.up.clade!=n.clade:
                clades[n.name] = {'clade_annotation': n.clade, 'clade_membership':n.clade}
            else:
                clades[n.name] = {'clade_membership':n.clade}
        else:
            clades[n.name] = {'clade_annotation': n.clade, 'clade_membership':n.clade}

    write_json({'nodes': clades}, args.output)

