"""
Annotate derived haplotypes per node from annotated clades and store as node data JSON.
"""
import argparse
from collections import Counter

from augur.reconstruct_sequences import load_alignments
from augur.utils import annotate_parents_for_tree, read_node_data, read_tree, write_json


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Annotate derived haplotypes per node from annotated clades and store as node data JSON",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--tree", required=True, help="Newick file from augur refine")
    parser.add_argument("--alignment", help="aligned HA1 sequences with internal nodes in FASTA format", required=True)
    parser.add_argument("--clades", help="clade annotations in node data JSON format", required=True)
    parser.add_argument("--clade-label-attribute", help="name of the branch attribute for clade labels in the given clades JSON", default="clade")
    parser.add_argument("--clade-node-attribute", help="name of the node attribute for clade membership in the given clades JSON", default="clade_membership")
    parser.add_argument("--attribute-name", default="haplotype", help="name of attribute to store the complete amino acid sequence of each node")
    parser.add_argument("--min-tips", type=int, default=1, help="minimum number of tips with a derived haplotype to include in final output. Haplotypes with fewer than this number of tips will be annotated with the haplotype of their immediate parent in the tree.")
    parser.add_argument("--output-node-data", help="JSON file with translated sequences by node", required=True)
    args = parser.parse_args()

    # Load tree.
    tree = read_tree(args.tree)
    tree = annotate_parents_for_tree(tree)

    # Load sequences.
    gene_name = "HA1"
    alignment = load_alignments([args.alignment], [gene_name])

    # Index aligned HA1 amino acid sequences by node name.
    sequence_by_node = {
        sequence.id: str(sequence.seq)
        for sequence in alignment[gene_name]
    }

    # Load clades.
    clades = read_node_data(args.clades)

    # Check clades for "branches" annotations from Augur 22.0.0 and onward.
    # Otherwise, look for "clade_annotation" keys in "nodes" annotations.
    if "branches" in clades:
        # Index amino acid sequence by clade label.
        sequence_by_clade = {
            node_data["labels"][args.clade_label_attribute]: sequence_by_node[node]
            for node, node_data in clades["branches"].items()
        }
    else:
        # Index amino acid sequence by clade annotation.
        sequence_by_clade = {
            node_data["clade_annotation"]: sequence_by_node[node]
            for node, node_data in clades["nodes"].items()
            if "clade_annotation" in node_data
        }

    # Index clade membership by node name.
    clade_by_node = {
        node: node_data[args.clade_node_attribute]
        for node, node_data in clades["nodes"].items()
    }

    # Find clade and derived haplotypes from clade sequences per node.
    haplotypes = {}
    count_by_haplotype = Counter()
    for node in tree.find_clades():
        haplotypes[node.name] = {}
        clade = clade_by_node[node.name]

        if clade == "unassigned" or sequence_by_node[node.name] == sequence_by_clade[clade]:
            # If the current node's sequence is identical to its clade's
            # sequence, store the clade name.
            haplotype = clade
        else:
            # If the current node's sequence differs from its clade's sequence,
            # find the differences between the two.
            mutations = []
            for i in range(len(sequence_by_node[node.name])):
                if sequence_by_node[node.name][i] != sequence_by_clade[clade][i]:
                    # Store ancestral allele, 1-based mutation position, and derived allele.
                    mutations.append(f"{sequence_by_clade[clade][i]}{i + 1}{sequence_by_node[node.name][i]}")

            # Store the clade name plus a delimited list of derived mutations
            # present in the current node.
            haplotype = f"{clade}:{'-'.join(mutations)}"

        # Store the clade and haplotype values for this node.
        haplotypes[node.name][args.attribute_name] = haplotype

        # Count haplotypes for tips to collapse rare haplotypes later.
        if node.is_terminal():
            count_by_haplotype[haplotype] += 1

    # Collapse rare haplotypes into their parents, if requested.
    if args.min_tips > 1:
        for node in tree.find_clades():
            if (count_by_haplotype[haplotypes[node.name][args.attribute_name]] < args.min_tips) and node.parent is not None:
                haplotypes[node.name][args.attribute_name] = haplotypes[node.parent.name][args.attribute_name]

    # Write out the node annotations.
    write_json({"nodes": haplotypes}, args.output_node_data)
