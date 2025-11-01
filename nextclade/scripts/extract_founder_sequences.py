

def read_founder_nodes(clades_json_file):
    import json
    with open(clades_json_file) as fh:
        clades_data = json.load(fh)
    founder_nodes = dict()
    for node_id, node_data in clades_data.get('branches', {}).items():
        if 'labels' in node_data and 'clade' in node_data['labels']:
            founder_nodes[node_id] = node_data['labels']['clade']

    return founder_nodes


def extract_sequences_from_node_data(node_clade_map, ancestral_json):
    import json
    from Bio import SeqRecord, Seq

    sequences = dict()
    with open(ancestral_json) as fh:
        ancestral_json = json.load(fh)['nodes']

    for node_id, node_data in ancestral_json.items():
        if node_id in node_clade_map:
            clade = node_clade_map[node_id]
            sequences[clade] = SeqRecord.SeqRecord(Seq.Seq(node_data['sequence']), id=clade, description='')

    return sequences


if __name__=="__main__":
    import argparse
    from Bio import SeqIO

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--clades-json", required=True, help="JSON file with clade definitions")
    parser.add_argument("--ancestral-json", required=True, help="Ancestral sequence JSON file")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file with founder sequences")
    args = parser.parse_args()

    node_clade_map = read_founder_nodes(args.clades_json)
    sequences = extract_sequences_from_node_data(node_clade_map, args.ancestral_json)

    with open(args.output_fasta, 'w') as fh:
        SeqIO.write([sequences[clade] for clade in sorted(sequences.keys())], fh, 'fasta')
