
'''
Creates full-genome fasta files for each clusters.
'''
import os
import json
import argparse

'''
Load clusters into usable python dictionary
'''
def clusters(files):
    diction = {}
    for fname in files:
        with open(fname) as jfile:
            json_data = json.load(jfile)
            for strain, cluster_dict in json_data["nodes"].items():
                cluster_number = cluster_dict['cluster']
                if cluster_number in diction:
                    diction[cluster_number][strain] = []
                else:
                    diction[cluster_number] = {strain : []}
    return diction


'''
Adds sequences into python dictionary by cluster
'''
def sequences(clusters, files):
    for fname in files:
        with open(fname) as jfile:
            json_data = json.load(jfile)
            for cluster_id, cluster_strains in clusters.items():
                for strain in cluster_strains.keys():
                    if strain in json_data["nodes"].keys():
                        seq = json_data["nodes"][strain]["sequence"]
                        clusters[cluster_id][strain].append(seq)
    for cluster_id, cluster_strains in clusters.items():
        for strain, sequence in cluster_strains.items():
            clusters[cluster_id][strain] = ''.join(sequence)
    return clusters


'''
Outputs python dictionary sequences into fasta files
'''
def genomes(seq_dict):
    for cluster, cluster_strains in seq_dict.items():
        with open ('results/genomes_{wildcards.lineage}_cluster'+cluster+'_{wildcards.resolution}.fasta', "w") as fasta: 
            for strain, seq in cluster_strains.items():
                fasta.write('>' + strain + ' | cluster'+cluster + '\n' + seq + '\n')
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create full-genome fasta files of clusters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--clusters', nargs='+', type = str, required=True, help = "list of cluster JSON files")
    parser.add_argument('--nt-muts', nargs='+', type =str, required=True, help = "list of nt-muts JSON files")
    args = parser.parse_args()

    #Create clusters dictionary
    cluster_dict = clusters(args.clusters)

    #Makes sequences dictionary
    cluster_seq_dict = sequences(cluster_dict, args.nt_muts)

    #Outputs fasta files
    genomes(cluster_seq_dict)

