import argparse
from augur.utils import write_json
import Bio.SeqIO
import datetime
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--titer-model", required=True, help="JSON with substitution effects from augur titers sub")
    parser.add_argument("--titers", required=True, help="TSV of raw titer measurements with reference strains in the `serum_strain` column")
    parser.add_argument("--forecasts", required=True, help="TSV with forecast frequencies per location and variant")
    parser.add_argument("--tip-attributes", required=True, help="TSV with tip name, date, region, and emerging haplotype")
    parser.add_argument("--ha1-sequences", required=True, help="FASTA with HA1 amino acid sequences per tip")
    parser.add_argument("--min-date", required=True, help="minimum date (YYYY-MM-DD) for tips to be included in analysis")
    parser.add_argument("--min-reference-year", type=int, required=True, help="minimum year for reference to be included in the analysis")
    parser.add_argument("--output-node-data", required=True, help="node data JSON of antigenic distances to the future by reference strain and boolean indicator `is_titer_reference` for filtering in Auspice")
    parser.add_argument("--output-table", required=True, help="TSV of antigenic distances to the future by reference strain")

    args = parser.parse_args()

    with open(args.titer_model, "r", encoding="utf-8") as fh:
        substitutions = json.load(fh)["substitution"]

    substitutions = {
        substitution.replace("HA1:", ""): weight
        for substitution, weight in substitutions.items()
    }
    print(f"Found {len(substitutions)} HA1 substitutions with antigenic weights")

    titers = pd.read_csv(
        args.titers,
        sep="\t",
    )

    date_range = set(range(args.min_reference_year, datetime.date.today().year + 1))
    references = [
        reference
        for reference in titers["serum_strain"].drop_duplicates().sort_values().values.tolist()
        if int(reference.split("/")[-1]) in date_range
    ]
    print(f"Found {len(references)} references")

    forecasts = pd.read_csv(
        args.forecasts,
        sep="\t",
    ).sort_values(
        "date",
        ascending=False,
    ).groupby([
        "location",
        "variant",
    ]).first().query(
        "median > 0"
    ).reset_index()
    forecasts["region_haplotype"] = forecasts.apply(lambda record: f"{record['location']}/{record['variant']}", axis=1)

    future_frequency_by_region_haplotype = dict(forecasts.loc[:, ["region_haplotype", "median"]].values)

    number_of_regions = forecasts["location"].drop_duplicates().shape[0]
    region_haplotypes = set(forecasts["region_haplotype"].drop_duplicates().values)
    print(f"Found {len(region_haplotypes)} region/haplotype combinations across {number_of_regions} regions")

    tips = pd.read_csv(
        args.tip_attributes,
        sep="\t",
        usecols=[
            "strain",
            "date",
            "region",
            "emerging_haplotype",
        ]
    ).query(
        f"date >= '{args.min_date}'"
    )
    tips.loc[tips["emerging_haplotype"] == "unassigned", "emerging_haplotype"] = "other"
    tips["region_haplotype"] = tips.apply(lambda record: f"{record['region']}/{record['emerging_haplotype']}", axis=1)

    tips_with_futures = tips[tips["region_haplotype"].isin(region_haplotypes)].copy()
    print(f"Found {tips_with_futures.shape[0]} tips")

    strains = set(tips_with_futures["strain"].values) | set(references)
    ha1_sequences_by_strain = {
        record.name: str(record.seq)
        for record in Bio.SeqIO.parse(args.ha1_sequences, "fasta")
        if record.name in strains
    }

    ha1_sequence_length = len(list(ha1_sequences_by_strain.values())[0])
    distance_by_reference = {}
    for reference in references:
        if reference not in ha1_sequences_by_strain:
            print(f"WARNING: Couldn't find a sequence for the reference '{reference}'")
            continue

        total_distance = 0
        reference_sequence = ha1_sequences_by_strain[reference]
        for region_haplotype in region_haplotypes:
            tips_for_region_haplotype = tips_with_futures.query(f"region_haplotype == '{region_haplotype}'")
            distance_for_region_haplotype = 0.0
            count_for_region_haplotype = tips_for_region_haplotype.shape[0]

            for tip in tips_for_region_haplotype["strain"]:
                if tip in ha1_sequences_by_strain:
                    tip_sequence = ha1_sequences_by_strain[tip]

                    for i in range(ha1_sequence_length):
                        if reference_sequence[i] != tip_sequence[i]:
                            distance_for_region_haplotype += substitutions.get(
                                f"{reference_sequence[i]}{i + 1}{tip_sequence[i]}",
                                0.0,
                            )

            # Calculate the average antigenic distance between the current reference
            # and all tips for this region and haplotype.
            average_distance = distance_for_region_haplotype / count_for_region_haplotype

            # Weight the average distance for this region and haplotype by the predicted
            # future frequency of this haplotpe in this region.
            weighted_average_distance = future_frequency_by_region_haplotype[region_haplotype] * average_distance
            total_distance += weighted_average_distance

        # Calculate the average weighted distance across all regions,
        # since the total weighted distance from above reflects all regions.
        distance_by_reference[reference] = total_distance / number_of_regions

    distances_df = pd.DataFrame([
        {
            "reference": reference,
            "distance": distance,
        }
        for reference, distance in distance_by_reference.items()
    ]).sort_values(
        "distance",
    )

    distances_df.to_csv(
        args.output_table,
        sep="\t",
        index=False,
    )

    node_data = {
        reference: {
            "antigenic_distance_to_future": distance,
            "is_titer_reference": True,
        }
        for reference, distance in distance_by_reference.items()
    }

    write_json({"nodes": node_data}, args.output_node_data)
