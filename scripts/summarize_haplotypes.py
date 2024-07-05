"""Summarize coverage of haplotypes by titers.
"""
import argparse
import json
import pandas as pd

from augur.io import read_metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="metadata TSV with derived haplotypes per strain")
    parser.add_argument("--frequencies", required=True, help="tip frequencies JSON")
    parser.add_argument("--output-table", required=True, help="TSV of haplotypes along with number of reference viruses, distinct reference viruses, number of test viruses, current frequency, and delta frequency from the last month.")
    parser.add_argument("--output-markdown-table", required=True, help="Markdown table of the TSV table above for use in narratives.")

    args = parser.parse_args()

    # Load metadata with derived haplotypes.
    metadata = read_metadata(args.metadata)

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
    for strain, record in metadata.iterrows():
        if strain not in current_frequency_by_strain:
            continue

        haplotype = record["haplotype"]

        if haplotype not in haplotype_frequencies:
            haplotype_frequencies[haplotype] = {
                "current_frequency": 0.0,
                "previous_frequency": 0.0,
            }

        haplotype_frequencies[haplotype]["current_frequency"] += current_frequency_by_strain[strain]
        haplotype_frequencies[haplotype]["previous_frequency"] += previous_frequency_by_strain[strain]

    # Build a data frame of haplotypes and their frequencies.
    # TODO: Report haplotype counts associated with frequencies.
    haplotype_records = []
    for haplotype in haplotype_frequencies:
        haplotype_records.append({
            "haplotype": haplotype,
            "current_frequency": haplotype_frequencies[haplotype]["current_frequency"],
            "delta_frequency": (haplotype_frequencies[haplotype]["current_frequency"] - haplotype_frequencies[haplotype]["previous_frequency"]),
        })

    haplotypes_df = pd.DataFrame(haplotype_records).set_index("haplotype")

    # Join haplotypes with titer reference counts and names.
    annotated_haplotypes = haplotypes_df
    # TODO: Annotate haplotypes with number and list of titer references and number of titer test strains.
    annotated_haplotypes = annotated_haplotypes.query("current_frequency > 0").copy()
    annotated_haplotypes = annotated_haplotypes.sort_values("current_frequency", ascending=False)

    annotated_haplotypes.to_csv(
        args.output_table,
        sep="\t",
        header=True,
        index=True,
    )

    # Create the Markdown version for display in narratives.
    annotated_haplotypes = annotated_haplotypes.reset_index()

    # Use spaces instead of commas, allowing Markdown to wrap
    # lines.
    annotated_haplotypes["haplotype"] = annotated_haplotypes["haplotype"].str.replace(",", " ")

    # Keep haplotypes with at least 1% frequency.
    annotated_haplotypes = annotated_haplotypes.query("current_frequency >= 0.01").copy()

    # Round frequencies prior to writing out the markdown table.
    annotated_haplotypes["current_frequency"] = (annotated_haplotypes["current_frequency"] * 100).round(0).astype(int)
    annotated_haplotypes["delta_frequency"] = (annotated_haplotypes["delta_frequency"] * 100).round(0).astype(int)

    # Save Markdown table.
    with open(args.output_markdown_table, "w", encoding="utf-8") as oh:
        print(
            annotated_haplotypes.to_markdown(
                index=False,
            ),
            file=oh,
        )
