from Bio import Phylo

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Root and sanitize a Newick tree")
    parser.add_argument("--input-tree", required=True, help="Input Newick tree file")
    parser.add_argument("--output-tree", required=True, help="Output rooted and sanitized Newick tree file")
    parser.add_argument("--root", required=False, help="Name of the reference sequence to root the tree on")
    parser.add_argument("--prune-length", type=float, default=0.03, help="Maximum branch length to retain (branches longer than this will be pruned)")
    args = parser.parse_args()

    # Load the tree
    T = Phylo.read(args.input_tree, "newick")

    # Root the tree
    if args.root:
        T.root_with_outgroup(args.root)

    to_prune = []
    for node in T.find_clades(order='postorder'):
        if node.is_terminal():
            node.total_length_per_child =  node.branch_length
            node.descendants = 1
        else:
            total_length = 0.0
            total_descendants = 0
            for child in node.clades:
                total_length += child.total_length_per_child
                total_descendants += child.descendants
                child.parent = node
            node.descendants = total_descendants
            node.total_length_per_child = total_length/total_descendants
        if node.total_length_per_child > args.prune_length:
            to_prune.append(node)

    for node in to_prune:
        if node == T.root:
            continue
        node.parent.clades = [c for c in node.parent.clades if c != node]

    T.ladderize()


    # Write the output tree
    Phylo.write(T, args.output_tree, "newick")