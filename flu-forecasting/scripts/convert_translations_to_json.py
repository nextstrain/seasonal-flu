"""
Convert translation FASTA to a node data JSON
"""
import argparse
from augur.reconstruct_sequences import load_alignments
from augur.utils import write_json
import Bio.Phylo


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert translation FASTA to a node data JSON",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--tree", required=True, help="Newick file for the tree used to construct the given node data JSONs")
    parser.add_argument("--alignment", nargs="+", help="sequence(s) to be used, supplied as FASTA files", required=True)
    parser.add_argument('--gene-names', nargs="+", type=str, help="names of the sequences in the alignment, same order assumed", required=True)
    parser.add_argument("--output", help="JSON file with translated sequences by node", required=True)
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data associated with internal nodes in the output JSON")
    parser.add_argument("--attribute-name", default="aa_sequence", help="name of attribute to store the complete amino acid sequence of each node")
    args = parser.parse_args()

    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Load sequences.
    alignments = load_alignments(args.alignment, args.gene_names)

    # Concatenate translated sequences into a single sequence indexed by sample name.
    is_node_terminal = {node.name: node.is_terminal() for node in tree.find_clades()}

    translations = {}
    for gene in args.gene_names:
        alignment = alignments[gene]

        for record in alignment:
            if is_node_terminal[record.name] or args.include_internal_nodes:
                # Initialize new samples by name with an empty string.
                if record.name not in translations:
                    translations[record.name] = {args.attribute_name:  ""}

                # Append the current gene's amino acid sequence to the current
                # string for this sample.
                translations[record.name][args.attribute_name] += str(record.seq)

    # Write out the node annotations.
    write_json({"nodes": translations}, args.output)
