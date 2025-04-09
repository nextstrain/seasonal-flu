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
    parser.add_argument("--titers", nargs="+", help="titers TSV files with columns named 'virus_strain' and 'serum_strain' representing test and reference strains, respectively")
    parser.add_argument("--titer-names", nargs="+", help="names of the titer collections provided to the `--titers` argument to use in the table output")
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
    haplotype_by_strain = {}
    for strain, record in metadata.iterrows():
        haplotype = record["haplotype"]
        haplotype_by_strain[strain] = haplotype

        if strain not in current_frequency_by_strain:
            continue

        if haplotype not in haplotype_frequencies:
            haplotype_frequencies[haplotype] = {
                "frequency": 0.0,
                "previous_frequency": 0.0,
            }

        haplotype_frequencies[haplotype]["frequency"] += current_frequency_by_strain[strain]
        haplotype_frequencies[haplotype]["previous_frequency"] += previous_frequency_by_strain[strain]

    # Build a data frame of haplotypes and their frequencies.
    # TODO: Report haplotype counts associated with frequencies.
    haplotype_records = []
    for haplotype in haplotype_frequencies:
        haplotype_records.append({
            "haplotype": haplotype,
            "frequency": haplotype_frequencies[haplotype]["frequency"],
            "delta_frequency": (haplotype_frequencies[haplotype]["frequency"] - haplotype_frequencies[haplotype]["previous_frequency"]),
        })

    haplotypes_df = pd.DataFrame(haplotype_records)

    # Load titer collections and annotate haplotypes by the number and list of
    # available references.
    annotated_haplotypes = haplotypes_df
    for titer_name, titer_collection in zip(args.titer_names, args.titers):
        titers = pd.read_csv(
            titer_collection,
            sep="\t",
            usecols=["virus_strain", "serum_strain"],
        ).drop_duplicates()

        # Strip the "-egg" suffix from egg-passaged titer references, so we can
        # find those strains in the haplotype frequencies.
        titers["serum_strain"] = titers["serum_strain"].str.replace("-egg", "")

        reference_counts = titers.groupby("serum_strain")["virus_strain"].count().reset_index(name="count")
        reference_counts["reference"] = reference_counts.apply(
            lambda record: f"{record['serum_strain']} ({record['count']})",
            axis=1,
        )
        reference_counts["haplotype"] = reference_counts["serum_strain"].map(haplotype_by_strain)

        haplotype_references = reference_counts.groupby("derived_haplotype").aggregate(
            total_references=("reference", "count"),
            distinct_references=("reference", "unique"),
        ).rename(columns={
            "total_references": f"{titer_name} count",
            "distinct_references": f"{titer_name} names"
        })

        # Annotate haplotypes with number and list of titer references and
        # number of titer test strains.
        annotated_haplotypes = annotated_haplotypes.merge(
            haplotype_references,
            on="derived_haplotype",
            how="left",
        )

    annotated_haplotypes = annotated_haplotypes.set_index("derived_haplotype")
    annotated_haplotypes = annotated_haplotypes.query("frequency > 0").copy()
    annotated_haplotypes = annotated_haplotypes.sort_values("frequency", ascending=False)

    for column in annotated_haplotypes.columns:
        if column.endswith("count"):
            annotated_haplotypes[column] = annotated_haplotypes[column].apply(lambda value: "" if pd.isnull(value) else int(value))

        if column.endswith("names"):
            annotated_haplotypes[column] = annotated_haplotypes[column].apply(lambda value: "\n".join(value) if str(value) != "nan" else "")

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
    annotated_haplotypes["derived_haplotype"] = annotated_haplotypes["derived_haplotype"].str.replace(",", " ")

    # Round frequencies prior to writing out the markdown table.
    annotated_haplotypes["frequency"] = (annotated_haplotypes["frequency"] * 100).round(2).astype(float)
    annotated_haplotypes["delta_frequency"] = (annotated_haplotypes["delta_frequency"] * 100).round(2).astype(float)

    # Save Markdown table.
    with open(args.output_markdown_table, "w", encoding="utf-8") as oh:
        print(
            annotated_haplotypes.fillna("").to_markdown(
                index=False,
            ),
            file=oh,
        )
