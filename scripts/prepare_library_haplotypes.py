#!/usr/bin/env python3
import argparse
import json
import pandas as pd

from augur.reconstruct_sequences import load_alignments


def load_distance_map(distance_map_path):
    with open(distance_map_path, "r", encoding="utf-8") as fh:
        distance_map = json.load(fh)

    return distance_map["name"], set(distance_map["map"]["HA1"].keys())


def main(args):
    columns_to_keep = [
        "strain",
        "gisaid_epi_isl",
        "date",
        "region",
        "country",
        "derived_haplotype",
        "HA1_mutations",
        "glycosylation_mutations",
        "recurrent_substitutions",
        "HA1_coverage",
        "SigPep_sequence",
        "HA1_sequence",
        "HA2_sequence",
    ]

    metadata = pd.read_csv(
        args.metadata,
        sep="\t",
        dtype="str",
        na_filter=False,
    )
    alignments = load_alignments(
        args.alignments,
        args.gene_names,
    )

    alignment_by_gene_and_strain = {}
    for gene, alignment in alignments.items():
        alignment_by_gene_and_strain[gene] = {}

        for sequence in alignment:
            alignment_by_gene_and_strain[gene][sequence.name] = str(sequence.seq)

    # Parse distance map(s).
    distance_map_by_name = {}
    for distance_map_path in args.distance_maps:
        distance_map_name, distance_map = load_distance_map(distance_map_path)
        distance_map_by_name[distance_map_name] = distance_map

    metadata["all_ha1_substitutions"] = metadata["aaSubstitutions"].apply(
        lambda subs: ",".join([
            sub.removeprefix("HA1:")
            for sub in subs.split(",")
            if sub.startswith("HA1:")
        ])
    )

    metadata["founder_ha1_substitutions"] = metadata["founderMuts['clade'].aaSubstitutions"].apply(
        lambda subs: ",".join([
            sub.removeprefix("HA1:")
            for sub in subs.split(",")
            if sub.startswith("HA1:")
        ])
    )

    metadata["derived_haplotype"] = (metadata["clade"] + ":" + metadata["founder_ha1_substitutions"]).str.rstrip(":")

    metadata["HA1_mutations"] = metadata["all_ha1_substitutions"].apply(
        lambda subs: len(subs.split(","))
    )

    for distance_map_name, distance_map in distance_map_by_name.items():
        mutations_column = f"{distance_map_name}_mutations"
        columns_to_keep.append(mutations_column)

        metadata[mutations_column] = metadata["all_ha1_substitutions"].apply(
            lambda subs: sum(
                sub[1:-1] in distance_map
                for sub in subs.split(",")
            )
        )

    metadata["glycosylation_mutations"] = metadata["glycosylation"].apply(
        lambda glyc: len([site for site in glyc.split(";") if site.startswith("HA1:")])
    )

    # Score recurrent substitutions.
    metadata["recurrent_substitutions"] = 0
    if args.recurrent_substitutions_map:
        with open(args.recurrent_substitutions_map, "r", encoding="utf-8") as fh:
            recurrent_substitutions_by_clade = json.load(fh)

        for subclade, subclade_substitutions in recurrent_substitutions_by_clade.items():
            metadata.loc[metadata["clade"] == subclade, "recurrent_substitutions"] = metadata.loc[
                metadata["clade"] == subclade,
                "founder_ha1_substitutions"
            ].apply(
                lambda subs: sum(
                    sub in subclade_substitutions
                    for sub in subs.split(",")
                )
            )

    metadata["HA1_coverage"] = metadata["cdsCoverage"].apply(
        lambda coverage: float(
            dict(
                gene_coverage.split(":") for gene_coverage in coverage.split(",")
            ).get("HA1", 0)
        )
    )

    for gene in alignment_by_gene_and_strain:
        metadata[f"{gene}_sequence"] = metadata["strain"].map(
            alignment_by_gene_and_strain[gene]
        )

    metadata.to_csv(
        args.output,
        sep="\t",
        columns=columns_to_keep,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="TSV of Nextstrain-style metadata merged with Nextclade annotations")
    parser.add_argument("--distance-maps", nargs="+", help="JSON(s) in Augur's distance map format. Each map is used to annotate mutation counts in the output with count columns named by the distance map.")
    parser.add_argument("--recurrent-substitutions-map", help="JSON of substitutions (<ancestral allele><position><derived allele>) indexed by clade name for the given subtype")
    parser.add_argument("--alignments", nargs="+", required=True, help="FASTA amino acid alignment(s) per gene")
    parser.add_argument("--gene-names", nargs="+", required=True, help="list of genes to match the FASTA alignments")
    parser.add_argument("--output", required=True, help="TSV of metadata annotated by mutation scores and amino acid alignments per gene")

    args = parser.parse_args()
    main(args)
