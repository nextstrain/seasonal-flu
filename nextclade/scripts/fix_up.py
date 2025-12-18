from Bio import Phylo



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Fix up tree")
    parser.add_argument("--build", required=True, help="Build identifier")


    args = parser.parse_args()

    if args.build == "h3n2/na/EPI1857215":
        tree = Phylo.read(f"build/{args.build}/tree_rooted.nwk", "newick")

        tree.root_with_outgroup("A/HongKong/1/1968")
        Phylo.write(tree, f"build/{args.build}/tree_rooted.nwk", "newick")

    if args.build == "vic/na/CY073894":
        tree = Phylo.read(f"build/{args.build}/tree_rooted.nwk", "newick")

        tree.root_with_outgroup("B/Memphis/6/1986")
        Phylo.write(tree, f"build/{args.build}/tree_rooted.nwk", "newick")


    if args.build == "vic/ha/KX058884":
        tree = Phylo.read(f"build/{args.build}/tree_rooted.nwk", "newick")

        tree.root_with_outgroup("B/Netherlands/1000/1977")
        Phylo.write(tree, f"build/{args.build}/tree_rooted.nwk", "newick")
