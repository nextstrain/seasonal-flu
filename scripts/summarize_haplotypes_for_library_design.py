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
        "emerging_haplotype": ("emerging_haplotype_ha", "first"),
        "n_countries": ("country", "nunique"),
        "n_regions": ("region", "nunique"),
        "latest_sequence": ("date", "max"),
        "representative_strain": ("strain", "first"),
        "representative_strain_accession": ("gisaid_epi_isl", "first"),
        "representative_strain_ha1_sequence": ("HA1_sequence", "first"),
        "representative_strain_ha1_sequence_hash": ("HA1_sequence_hash", "first"),
    }

    if "lcr_score" in data_frame.columns:
       aggregated_fields["lcr_score"] = ("lcr_score", "mean")

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

    # Apply a minimum coverage threshold to focus on haplotypes with more
    # complete HA1 sequences.
    data_frame_by_lineage = data_frame_by_lineage[
        data_frame_by_lineage["HA1_coverage"] >= args.min_ha1_coverage
    ].copy()

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

    def get_fitnesses_by_haplotype_resolution(haplotype_fitnesses, variant_column, lineage_summary):
        fitnesses = pd.read_csv(haplotype_fitnesses, sep="\t")
        hier_fitnesses = fitnesses.loc[
            fitnesses["location"] == "hierarchical", ["variant", "median"]
        ].rename(
            columns={
                "variant": variant_column,
                "median": "median_ga",
            }
        )
        hier_fitnesses[variant_column] = hier_fitnesses[variant_column].str.replace("-", ",")

        if variant_column != "derived_haplotype":
            counts = lineage_summary.groupby(variant_column).aggregate(
                count=("count", "sum"),
            ).reset_index()
        else:
            counts = lineage_summary.reset_index().loc[:, [variant_column, "count"]]

        # Calculate the weighted population average GA.
        hier_fitnesses = hier_fitnesses.merge(
            counts,
            on=variant_column,
        )
        hier_fitnesses["frequency"] = (
            hier_fitnesses["count"] / hier_fitnesses["count"].sum()
        )

        weighted_mean_population_ga = (
            hier_fitnesses["frequency"] * hier_fitnesses["median_ga"]
        ).sum()

        hier_fitnesses[f"relative_fitness_{variant_column}"] = np.log(hier_fitnesses["median_ga"]) - np.log(weighted_mean_population_ga)
        hier_fitnesses = hier_fitnesses.drop(columns=["count", "median_ga", "frequency"])

        lineage_summary = lineage_summary.merge(
            hier_fitnesses,
            how="left",
            on=variant_column,
        )
        lineage_summary[f"relative_fitness_{variant_column}"] = lineage_summary[f"relative_fitness_{variant_column}"].fillna(0.0)

        lineage_summary[f"score_fitness_{variant_column}"] = (
            lineage_summary[f"relative_fitness_{variant_column}"]
            / lineage_summary[f"relative_fitness_{variant_column}"].abs().max()
        )

        return lineage_summary.copy()

    lineage_summary = get_fitnesses_by_haplotype_resolution(args.aa_haplotype_fitnesses, "derived_haplotype", lineage_summary)
    lineage_summary = get_fitnesses_by_haplotype_resolution(args.emerging_haplotype_fitnesses, "emerging_haplotype", lineage_summary)

    # Annotate summary by relative fitness for different haplotype resolutions.

    # Score haplotypes following approach for clade assignment in:
    # https://github.com/influenza-clade-nomenclature/clade-suggestion-algorithm
    lineage_summary["score_mutations"] = 0.0
    for mutation_column, mutation_weight in mutation_weights.items():
        lineage_summary["score_mutations"] = lineage_summary["score_mutations"] + (lineage_summary[mutation_column] * mutation_weight)

    lineage_summary["score_mutations"] = lineage_summary["score_mutations"] / (
        lineage_summary["score_mutations"].median()
        + lineage_summary["score_mutations"]
    )

    lineage_summary["score_lcr"] = (
        lineage_summary["lcr_score"]
        / lineage_summary["lcr_score"].abs().max()
    )

    lineage_summary["score_total"] = (
        lineage_summary["score_mutations"]
        + lineage_summary["score_lcr"]
        + lineage_summary["score_fitness_derived_haplotype"]
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
       float_format="%.2f",
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="TSV of metadata for potential library strains with derived haplotypes and scores annotated")
    parser.add_argument("--emerging-haplotype-fitnesses", help="TSV of fitnesses per emerging haplotype (coarser haplotypes)")
    parser.add_argument("--aa-haplotype-fitnesses", help="TSV of fitnesses per amino acid haplotype (finer haplotypes)")
    parser.add_argument("--mutation-weights", nargs="+", help="key/value pairs of columns in the metadata to score and their corresponding weight (e.g., 'koel_mutations=3 wolf_mutations=2 HA1_mutations=1')")
    parser.add_argument("--min-ha1-coverage", type=float, default=0.75, help="minimum coverage for HA1 as a proportion of the full amino acid sequence length")
    parser.add_argument("--output", required=True, help="TSV of derived haplotypes summaries with representative strains/sequences and scores")

    args = parser.parse_args()
    main(args)
