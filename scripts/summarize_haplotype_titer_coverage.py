"""Summarize coverage of haplotypes by titers.
"""
import argparse
from augur.utils import read_node_data, write_json
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--haplotypes", required=True, help="node data JSON of derived haplotypes per strain")
    parser.add_argument("--antigenic-distances", required=True, help="TSV of antigenic distances between reference and test viruses with haplotype annotations")
    parser.add_argument("--frequencies", required=True, help="tip frequencies JSON")
    parser.add_argument("--attribute-name-for-haplotype-without-reference", default="haplotype_missing_reference_virus", help="name for attribute indicating whether haplotypes are missing a reference virus in the node data JSON output.")
    parser.add_argument("--output-table", required=True, help="TSV of haplotypes along with number of reference viruses, distinct reference viruses, number of test viruses, current frequency, and delta frequency from the last month.")
    parser.add_argument("--output-markdown-table", required=True, help="Markdown table of the TSV table above for use in narratives.")
    parser.add_argument("--output-node-data", required=True, help="node data JSON of non-zero haplotypes without reference viruses")

    args = parser.parse_args()

    # Load haplotypes.
    haplotypes = read_node_data(args.haplotypes)
    haplotypes = haplotypes["nodes"]

    # Load antigenic distances.
    distances = pd.read_csv(args.antigenic_distances, sep="\t")
    references_by_haplotype = distances.loc[
        :,
        ["reference_strain", "haplotype_reference"]
    ].drop_duplicates().groupby(
        "haplotype_reference",
    ).aggregate(
        total_references=("reference_strain", "count"),
        distinct_references=("reference_strain", "unique"),
    )

    # Load frequencies.
    with open(args.frequencies) as fh:
        frequencies_data = json.load(fh)

    # Get the current and previous frequency per strain.
    current_frequency_by_strain = {}
    previous_frequency_by_strain = {}
    for strain, strain_data in frequencies_data.items():
        if "frequencies" in strain_data:
            current_frequency_by_strain[strain] = strain_data["frequencies"][-1]
            previous_frequency_by_strain[strain] = strain_data["frequencies"][-2]

    # Calculate haplotype frequencies and delta frequencies.
    haplotype_frequencies = {}
    for strain, node_data in haplotypes.items():
        if strain not in current_frequency_by_strain:
            continue

        haplotype = node_data["haplotype"]

        if haplotype not in haplotype_frequencies:
            haplotype_frequencies[haplotype] = {
                "current_frequency": 0.0,
                "previous_frequency": 0.0,
            }

        haplotype_frequencies[haplotype]["current_frequency"] += current_frequency_by_strain[strain]
        haplotype_frequencies[haplotype]["previous_frequency"] += previous_frequency_by_strain[strain]

    # Build a data frame of haplotypes and their frequencies.
    haplotype_records = []
    for haplotype in haplotype_frequencies:
        haplotype_records.append({
            "haplotype": haplotype,
            "current_frequency": haplotype_frequencies[haplotype]["current_frequency"],
            "delta_frequency": (haplotype_frequencies[haplotype]["current_frequency"] - haplotype_frequencies[haplotype]["previous_frequency"]),
        })

    haplotypes_df = pd.DataFrame(haplotype_records).set_index("haplotype")

    # Join haplotypes with titer reference counts and names.
    annotated_haplotypes = haplotypes_df.join(
        references_by_haplotype,
        how="left",
    )
    annotated_haplotypes["total_references"] = annotated_haplotypes["total_references"].fillna(0).astype(int)
    annotated_haplotypes = annotated_haplotypes.query("current_frequency > 0").copy()
    annotated_haplotypes = annotated_haplotypes.sort_values("current_frequency", ascending=False)

    annotated_haplotypes.to_csv(
        args.output_table,
        sep="\t",
        header=True,
        index=True,
    )

    # Save non-zero frequency haplotypes with no references.
    nonzero_haplotypes_per_node = {}
    nonzero_haplotypes = set(annotated_haplotypes.loc[
        annotated_haplotypes["total_references"] == 0,
    ].index.values)

    for strain, node_data in haplotypes.items():
        if node_data["haplotype"] in nonzero_haplotypes:
            nonzero_haplotypes_per_node[strain] = {
                args.attribute_name_for_haplotype_without_reference: "true"
            }

    write_json({"nodes": nonzero_haplotypes_per_node}, args.output_node_data)

    # Create the Markdown version for display in narratives.
    annotated_haplotypes = annotated_haplotypes.reset_index()

    # Use spaces instead of commas, allowing Markdown to wrap
    # lines.
    annotated_haplotypes["haplotype"] = annotated_haplotypes["haplotype"].str.replace(",", " ")

    # Fill missing distinct references with empty string.
    annotated_haplotypes["distinct_references"] = annotated_haplotypes["distinct_references"].fillna("")

    # Round deminal values of frequencies.
    annotated_haplotypes["current_frequency"] = annotated_haplotypes["current_frequency"].round(2)

    # Keep non-zero frequencies.
    annotated_haplotypes = annotated_haplotypes.query("current_frequency > 0").copy()

    # Round change in frequency.
    annotated_haplotypes["delta_frequency"] = annotated_haplotypes["delta_frequency"].round(2)

    # Save Markdown table.
    with open(args.output_markdown_table, "w", encoding="utf-8") as oh:
        print(annotated_haplotypes.to_markdown(index=False), file=oh)
