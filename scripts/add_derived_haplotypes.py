"""
Annotate derived haplotypes per node from annotated clades and store as node data JSON.
"""
import argparse
import json
import pandas as pd


def create_haplotype_for_record(record, clade_column, mutations_column, genes=None, strip_genes=False, sites_by_gene=None):
    """Create a haplotype string for the given record based on the values in its
    clade and mutations column. If a list of genes is given, filter mutations to
    only those in the requested genes.

    """
    clade = record[clade_column]

    if record[mutations_column] == "":
        return clade

    mutations = record[mutations_column].split(",")

    # Filter mutations to requested genes.
    if sites_by_gene is not None:
        filtered_mutations = []
        for mutation in mutations:
            # mutation looks like "HA1:N145K"
            gene, allele = mutation.split(":")
            position = allele[1:-1]
            if gene in sites_by_gene and position in sites_by_gene[gene]:
                filtered_mutations.append(mutation)

        mutations = filtered_mutations
    elif genes is not None:
        mutations = [
            mutation
            for mutation in mutations
            if mutation.split(":")[0] in genes
        ]

    mutations = "-".join(mutations).replace(":", "-")

    if mutations:
        if strip_genes and genes is not None:
            for gene in genes:
                mutations = mutations.replace(f"{gene}-", "")

        return f"{clade}:{mutations}"
    else:
        return clade


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Annotate derived haplotypes per record in Nextclade annotations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--nextclade", required=True, help="TSV file of Nextclade annotations with columns for clade and AA mutations derived from clade")
    parser.add_argument("--clade-column", help="name of the branch attribute for clade labels in the given Nextclade annotations", default="subclade")
    parser.add_argument("--mutations-column", help="name of the attribute for mutations relative to clades in the given Nextclade annotations", default="founderMuts['subclade'].aaSubstitutions")
    parser.add_argument("--genes", nargs="+", help="list of genes to filter mutations to. If not provided, all mutations will be used.")
    parser.add_argument("--distance-map", help="distance map JSON of genes and positions to include in haplotypes")
    parser.add_argument("--strip-genes", action="store_true", help="strip gene names from coordinates in output haplotypes")
    parser.add_argument("--attribute-name", default="haplotype", help="name of attribute to store the derived haplotype in the output file")
    parser.add_argument("--output", help="TSV file of Nextclade annotations with derived haplotype column added", required=True)
    args = parser.parse_args()

    # Load Nextclade annotations.
    df = pd.read_csv(
        args.nextclade,
        sep="\t",
        dtype={
            args.clade_column: "str",
            args.mutations_column: "str",
        },
        na_filter=False,
    )

    # Load distance map.
    sites_by_gene = None
    if args.distance_map:
        with open(args.distance_map, "r", encoding="utf-8") as fh:
            distance_map = json.load(fh)
            sites_by_gene = distance_map["map"]

    # Annotate derived haplotypes.
    df[args.attribute_name] = df.apply(
        lambda record: create_haplotype_for_record(
            record,
            args.clade_column,
            args.mutations_column,
            args.genes,
            args.strip_genes,
            sites_by_gene,
        ),
        axis=1
    )

    # Save updated Nextclade annotations
    df.to_csv(
        args.output,
        sep="\t",
        index=False,
    )
