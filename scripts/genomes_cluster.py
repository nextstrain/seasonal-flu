
'''
Creates full-genome fasta files of clusters.
'''
import os
import json


'''
Load .json clusters into python dictionary
'''
def clusters(file):
    diction = {}
    if os.path.isfile(file):
        with open(file) as jfile:
            json_data = json.load(jfile)
            for strain, cluster_dict in json_data["nodes"].items():
                for cluster_number in cluster_dict.values():
                    if cluster_number in diction:
                        diction[cluster_number].append(strain)
                    else:
                        diction[cluster_number] = [strain]
    for cluster_num, strains in diction.items():
        diction[cluster_num] = dict.fromkeys(strains, [])
    return diction


'''
Loads sequences into python dictionary by cluster
'''
def data(files):
    for fname in files:
        if os.path.isfile(fname):
            with open(fname) as jfile:
                json_data = json.load(jfile)
    return json_data


def sequences(clusters, file):
    data = clusters
    for cluster_id, cluster_strains in clusters.items():
        for strain, emptyset in cluster_strains.items():
            if strain in file["nodes"].keys():
                seq = file["nodes"][strain]["sequence"]
                #if strain == 'A/Louisiana/1/2019':
                 #   print(seq)
                data[cluster_id][strain].append(seq)
    return data

#Create clusters dictionary
y = clusters('clustering_h3n2_2y.json')

#Create nt_muts list
nt_muts = ['nt-muts_h3n2_ha_2y.json']

#Makes data to feed into sequences
z = data(nt_muts)
#print(z)

#There are not 335 of these strains with name "A/Calif..." in either z or y.
#set = []
#for strain in z["nodes"]:
 #   if "A/California/141/2019" in strain:
  #      set.append(1)
#print(set)

#print(z['nodes']['0622b49d-dada-40cb-9ba7-b7a75b0d9ab5'])

#There is only one sequence associated with z['nodes'][c]
for cluster_id, cluster_strains in y.items():
    for strain, emptyset in cluster_strains.items():
        if strain in z["nodes"].keys():
            print(strain, y[cluster_id][strain])

#Makes sequences dictionary
#x = sequences(y, z)
#print(x)
#print([(cluster, strain, len(seqs)) for cluster, cluster_strains in x.items() for strain, seqs in cluster_strains.items()])
#print([(cluster, len(cluster_strains)) for cluster, cluster_strains in x.items()]) 

'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create full-genome fasta files of clusters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--clusters', nargs='+', type = str, required=True, help = "list of cluster JSON files")
    args = parser.parse_args()

    # print out dict
    print(clusters(args.clusters))
'''
