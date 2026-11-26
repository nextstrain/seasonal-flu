import json

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Adjust branch lengths based on clade information")
    parser.add_argument("--branch-lengths", required=True, help="Branch lengths file in JSON format")
    parser.add_argument("--clade-json", required=True, help="Clade information in JSON format")
    parser.add_argument("--max-branch-length", type=float, required=True, help="Maximum allowed branch length")
    parser.add_argument("--output-node-data", required=True, help="Output node data file in JSON format")

    args = parser.parse_args()

    # Load input data
    with open(args.branch_lengths, 'r') as f:
        node_data = json.load(f)

    if args.max_branch_length <= 0:
        with open(args.output_node_data, 'w') as f:
            json.dump(node_data, f)

        exit(0)

    with open(args.clade_json, 'r') as f:
        clade_data = json.load(f)['nodes']

    # gather all unassigned nodes
    unassigned_nodes = []
    for node_id, clade_info in clade_data.items():
        if clade_info.get('clade_membership', 'unassigned') == 'unassigned':
            unassigned_nodes.append(node_id)

    unassigned_nodes = set(unassigned_nodes)

    # Adjust branch lengths
    for node_id, data in node_data['nodes'].items():
        if node_id in unassigned_nodes:
            data['branch_length'] = min(data.get('branch_length', 0), args.max_branch_length)

    # Save adjusted node data (currently unchanged)
    with open(args.output_node_data, 'w') as f:
        json.dump(node_data, f, indent=2)