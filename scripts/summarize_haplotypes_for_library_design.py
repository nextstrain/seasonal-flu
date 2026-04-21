import argparse
import datetime
import hashlib
import json
import numpy as np
import pandas as pd


def summarize_haplotypes(data_frame, mutation_columns):
    strain_to_accession = dict(data_frame.loc[:, ["strain", "gisaid_epi_isl"]].values)

    aggregated_fields = {
        "count": ("derived_haplotype", "count"),
        "n_countries": ("country", "nunique"),
        "n_regions": ("region", "nunique"),
        "latest_sequence": ("date", "max"),
        "representative_strain": ("strain", "first"),
        "representative_strain_accession": ("gisaid_epi_isl", "first"),
        "representative_strain_ha1_sequence": ("HA1_sequence", "first"),
        "representative_strain_ha1_sequence_hash": ("HA1_sequence_hash", "first"),
    }

    if mutation_columns is not None:
        for mutation_column in mutation_columns:
            aggregated_fields[mutation_column] = (mutation_column, "mean")

    summary = data_frame.sort_values("HA1_coverage", ascending=False).groupby("derived_haplotype", sort=False).aggregate(
        **aggregated_fields
    ).sort_values(
        "count",
        ascending=False,
    )

    ha2_sequence_summary = data_frame.groupby(["derived_haplotype", "HA2_sequence"]).aggregate(
        HA2_sequence_count=("strain", "count"),
    ).reset_index().sort_values(
        "HA2_sequence_count",
        ascending=False,
    ).groupby("derived_haplotype").aggregate(
        representative_ha2_sequence=("HA2_sequence", "first"),
        representative_ha2_sequence_count=("HA2_sequence_count", "first"),
    )

    summary = summary.merge(
        ha2_sequence_summary,
        on="derived_haplotype",
    )

    total_countries = data_frame["country"].drop_duplicates().shape[0]
    total_regions = data_frame["region"].drop_duplicates().shape[0]
    summary["proportion_countries"] = summary["n_countries"] / total_countries
    summary["proportion_regions"] = summary["n_regions"] / total_regions

    summary["days_since_latest_sequence"] = (pd.to_datetime(datetime.date.today()) - pd.to_datetime(summary["latest_sequence"])).dt.days

    return summary


def main(args):
    # Load metadata with scores per strain and haplotype.
    data_frame_by_lineage = pd.read_csv(
        args.metadata,
        sep="\t",
        na_filter=False,
    )

    # We allow records with collection dates that have ambiguous days (e.g.,
    # 2020-01-XX), so we need to convert those ambiguous characters to valid
    # date strings. We might prefer haplotypes with more recent collection
    # dates, so use the first day of the month to avoid overestimating when a
    # sample was collected.
    data_frame_by_lineage["date"] = data_frame_by_lineage["date"].str.replace(
        "XX", "01"
    )

    # Hash the HA1 sequence to get a simpler identifier for the sequence. This
    # is useful for deduplication or communicating about specific sequences
    # without resharing the whole sequence.
    data_frame_by_lineage["HA1_sequence_hash"] = data_frame_by_lineage[
        "HA1_sequence"
    ].apply(lambda sequence: hashlib.sha256(sequence.encode()).hexdigest()[:7])

    # Create a map of mutation columns and their corresponding weights for the
    # scoring below.
    mutation_weights = {}
    mutation_columns = None
    if args.mutation_weights is not None:
        for mutation_weight in args.mutation_weights:
            column, weight = mutation_weight.split("=")
            weight = int(weight)
            mutation_weights[column] = weight

        mutation_columns = list(mutation_weights.keys())

    lineage_summary = summarize_haplotypes(
        data_frame_by_lineage,
        mutation_columns,
    )

    # Load fitness for lineage.
    fitnesses = pd.read_csv(args.fitnesses, sep="\t")
    hier_fitnesses = fitnesses.loc[
        fitnesses["location"] == "hierarchical", ["variant", "median"]
    ].rename(
        columns={
            "variant": "derived_haplotype",
            "median": "median_ga",
        }
    )
    hier_fitnesses["derived_haplotype"] = hier_fitnesses[
        "derived_haplotype"
    ].str.replace("-", ",")

    # Calculate the weighted population average GA.
    hier_fitnesses = hier_fitnesses.merge(
        lineage_summary.reset_index().loc[:, ["derived_haplotype", "count"]],
        on="derived_haplotype",
    )
    hier_fitnesses["frequency"] = (
        hier_fitnesses["count"] / hier_fitnesses["count"].sum()
    )
    weighted_mean_population_ga = (
        hier_fitnesses["frequency"] * hier_fitnesses["median_ga"]
    ).sum()

    hier_fitnesses["relative_fitness"] = np.log(
        hier_fitnesses["median_ga"]
    ) - np.log(weighted_mean_population_ga)
    hier_fitnesses = hier_fitnesses.drop(columns=["count"])

    # Annotate summary by median hierarchical GA.
    lineage_summary = lineage_summary.merge(
        hier_fitnesses,
        how="left",
        on="derived_haplotype",
    )
    lineage_summary["median_ga"] = lineage_summary["median_ga"].fillna(0.0)
    lineage_summary["relative_fitness"] = lineage_summary[
        "relative_fitness"
    ].fillna(0.0)

    # Score haplotypes following approach for clade assignment in:
    # https://github.com/influenza-clade-nomenclature/clade-suggestion-algorithm
    lineage_summary["score_mutations"] = 0.0
    for mutation_column, mutation_weight in mutation_weights.items():
        lineage_summary["score_mutations"] = lineage_summary["score_mutations"] + (lineage_summary[mutation_column] * mutation_weight)

    lineage_summary["score_mutations"] = lineage_summary["score_mutations"] / (
        lineage_summary["score_mutations"].median()
        + lineage_summary["score_mutations"]
    )

    lineage_summary["score_fitness"] = (
        lineage_summary["relative_fitness"]
        / lineage_summary["relative_fitness"].abs().max()
    )

    lineage_summary["score_total"] = (
        lineage_summary["score_mutations"]
        + lineage_summary["score_fitness"]
        + lineage_summary["proportion_countries"]
    )

    lineage_summary = lineage_summary.sort_values(
        "score_total",
        ascending=False,
    )
    lineage_summary.to_csv(
       args.output,
       sep="\t",
       index=False,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="TSV of metadata for potential library strains with derived haplotypes and scores annotated")
    parser.add_argument("--fitnesses", help="TSV of fitnesses per derived haplotype")
    parser.add_argument("--mutation-weights", nargs="+", help="key/value pairs of columns in the metadata to score and their corresponding weight (e.g., 'koel_mutations=3 wolf_mutations=2 HA1_mutations=1')")
    parser.add_argument("--output", required=True, help="TSV of derived haplotypes summaries with representative strains/sequences and scores")

    args = parser.parse_args()
    main(args)
