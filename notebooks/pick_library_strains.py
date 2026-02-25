import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    import Bio.SeqIO
    import datetime
    import json
    import numpy as np
    import pandas as pd
    return Bio, datetime, json, mo, np, pd


@app.cell(hide_code=True)
def _(
    haplotypes_by_lineage,
    min_date,
    min_sequences,
    mo,
    nonsingleton_haplotypes_by_lineage,
):
    mo.md(rf"""
    ## Setup

    Selected all GISAID isolates with:

      - minimum collection date of {min_date}
      - an HA record for the isolate
      - a “good” overall QC status by Nextclade alignment

    Excluded records with:

      - ambiguous month in collection date (e.g., 2025-XX-XX)
      - excluded records with egg passaging

    ## Summary

    For H1N1pdm, found {haplotypes_by_lineage['h1n1pdm']} distinct haplotypes including {nonsingleton_haplotypes_by_lineage['h1n1pdm']} with at least {min_sequences} sequences.

    For H3N2, found {haplotypes_by_lineage['h3n2']} distinct haplotypes including {nonsingleton_haplotypes_by_lineage['h3n2']} with at least {min_sequences} sequences.
    """)
    return


@app.cell
def _():
    pivot_by_lineage = {
        "h1n1pdm": "D.3.1",
        "h3n2": "J.2.4",
    }
    return (pivot_by_lineage,)


@app.cell
def _(mo):
    min_date_picker = mo.ui.date(
        start="2025-05-01",
        value="2025-07-01",
        label="min date:",
    )
    min_date_picker
    return (min_date_picker,)


@app.cell
def _(min_date_picker):
    min_date = min_date_picker.value.strftime("%Y-%m-%d")
    return (min_date,)


@app.cell
def _(mo):
    min_sequences_picker = mo.ui.number(
        start=1,
        stop=20,
        value=1,
        label="min sequences:",
    )
    min_sequences_picker
    return (min_sequences_picker,)


@app.cell
def _(min_sequences_picker):
    min_sequences = min_sequences_picker.value
    return (min_sequences,)


@app.cell
def _():
    paths_by_lineage = {
        "h1n1pdm": {
            "distance_map": "config/distance_maps/h1n1pdm/ha/canton.json",
            "metadata": "builds/h1n1pdm/metadata.tsv",
            "HA": "builds/h1n1pdm/ha/translations/HA.fasta",
            "ga": "h1n1pdm_ga.tsv",
            "output": "h1n1pdm_haplotypes.tsv",
        },
        "h3n2": {
            "distance_map": "config/distance_maps/h3n2/ha/wolf.json",
            "rbs_distance_map": "config/distance_maps/h3n2/ha/koel.json",
            "metadata": "builds/h3n2/metadata.tsv",
            "HA": "builds/h3n2/ha/translations/HA.fasta",
            "ga": "h3n2_ga.tsv",
            "output": "h3n2_haplotypes.tsv",
        },
    }
    return (paths_by_lineage,)


@app.cell
def _(json):
    with open("substitutions_for_library_selection.json", "r", encoding="utf-8") as fh:
        recurrent_substitutions_lineage_and_subclade = json.load(fh)
    return (recurrent_substitutions_lineage_and_subclade,)


@app.cell
def _(recurrent_substitutions_lineage_and_subclade):
    recurrent_substitutions_lineage_and_subclade
    return


@app.cell
def _(json):
    def load_distance_map(distance_map_path):
        with open(distance_map_path, "r", encoding="utf-8") as fh:
            distance_map = json.load(fh)

        return set(distance_map["map"]["HA1"].keys())
    return (load_distance_map,)


@app.cell
def _(pd):
    def load_metadata(metadata_path, min_date=None, epitope_sites=None, rbs_sites=None, recurrent_substitutions_by_clade=None):
        if epitope_sites is None:
            epitope_sites = set()

        if rbs_sites is None:
            rbs_sites = set()

        df = pd.read_csv(
            metadata_path,
            sep="\t",
            dtype="str",
            na_filter=False,
        )

        df["date"] = df["date"].str.replace("XX", "01")

        if min_date:
            df = df[df["date"] >= min_date].copy()

        df["all_ha1_substitutions"] = df["aaSubstitutions"].apply(
            lambda subs: ",".join([
                sub.replace("HA1:", "")
                for sub in subs.split(",")
                if sub.startswith("HA1:")
            ])
        )

        df["founder_ha1_substitutions"] = df["founderMuts['clade'].aaSubstitutions"].apply(
            lambda subs: ",".join([
                sub.replace("HA1:", "")
                for sub in subs.split(",")
                if sub.startswith("HA1:")
            ])
        )

        df["derived_haplotype"] = (df["clade"] + ":" + df["founder_ha1_substitutions"]).str.rstrip(":")

        df["ha1_mutations_from_subclade"] = df["all_ha1_substitutions"].apply(
            lambda subs: len(subs.split(","))
        )

        df["epitope_mutations_from_subclade"] = df["all_ha1_substitutions"].apply(
            lambda subs: sum(
                sub[1:-1] in epitope_sites
                for sub in subs.split(",")
            )
        )

        df["rbs_mutations_from_subclade"] = df["all_ha1_substitutions"].apply(
            lambda subs: sum(
                sub[1:-1] in rbs_sites
                for sub in subs.split(",")
            )
        )

        df["glycosylation_count"] = df["glycosylation"].apply(
            lambda glyc: len([site for site in glyc.split(";") if site.startswith("HA1:")])
        )

        df["recurrent_substitutions"] = 0
        for subclade, subclade_substitutions in recurrent_substitutions_by_clade.items():
            df.loc[df["clade"] == subclade, "recurrent_substitutions"] = df.loc[
                df["clade"] == subclade,
                "founder_ha1_substitutions"
            ].apply(
                lambda subs: sum(
                    sub in subclade_substitutions
                    for sub in subs.split(",")
                )
            )

        return df
    return (load_metadata,)


@app.cell
def _(Bio, paths_by_lineage):
    alignment_by_lineage_by_strain = {}
    for _lineage, _lineage_paths in paths_by_lineage.items():
        alignment_by_lineage_by_strain[_lineage] = {
            sequence.name: str(sequence.seq)
            for sequence in Bio.SeqIO.parse(_lineage_paths["HA"], "fasta")
        }
    return (alignment_by_lineage_by_strain,)


@app.cell
def _(
    alignment_by_lineage_by_strain,
    load_distance_map,
    load_metadata,
    min_date,
    min_sequences,
    np,
    paths_by_lineage,
    recurrent_substitutions_lineage_and_subclade,
):
    data_frame_by_lineage = {}
    haplotypes_by_lineage = {}
    nonsingleton_haplotypes_by_lineage = {}

    for lineage, lineage_paths in paths_by_lineage.items():
        epitope_sites = load_distance_map(lineage_paths["distance_map"])
        rbs_sites = load_distance_map(lineage_paths["rbs_distance_map"]) if "rbs_distance_map" in lineage_paths else None

        data_frame_by_lineage[lineage] = load_metadata(
            lineage_paths["metadata"],
            min_date,
            epitope_sites,
            rbs_sites,
            recurrent_substitutions_lineage_and_subclade[lineage],
        )

        data_frame_by_lineage[lineage]["strain_ha1_gaps"] = data_frame_by_lineage[lineage]["strain"].apply(
            lambda strain: (alignment_by_lineage_by_strain[lineage][strain].count("X") + alignment_by_lineage_by_strain[lineage][strain].count("-")) if strain in alignment_by_lineage_by_strain[lineage] else np.inf
        )

        haplotypes_by_lineage[lineage] = (data_frame_by_lineage[lineage]["derived_haplotype"].value_counts() > 0).sum()
        nonsingleton_haplotypes_by_lineage[lineage] = (data_frame_by_lineage[lineage]["derived_haplotype"].value_counts() >= min_sequences).sum()
    return (
        data_frame_by_lineage,
        haplotypes_by_lineage,
        nonsingleton_haplotypes_by_lineage,
    )


@app.cell
def _(datetime, pd):
    def summarize_haplotypes(data_frame, alignment_by_strain):
        strain_to_accession = dict(data_frame.loc[:, ["strain", "gisaid_epi_isl"]].values)

        summary = data_frame.sort_values("strain_ha1_gaps").groupby("derived_haplotype", sort=False).aggregate(
            count=("derived_haplotype", "count"),
            mean_ha1_mutations=("ha1_mutations_from_subclade", "mean"),
            mean_epitope_mutations=("epitope_mutations_from_subclade", "mean"),
            mean_rbs_mutations=("rbs_mutations_from_subclade", "mean"),
            mean_glycosylation_count=("glycosylation_count", "mean"),
            mean_recurrent_mutations=("recurrent_substitutions", "mean"),
            n_countries=("country", "nunique"),
            n_regions=("region", "nunique"),
            latest_sequence=("date", "max"),
            representative_strain=("strain", "first"),
        ).sort_values(
            "count",
            ascending=False,
        )

        total_countries = data_frame["country"].drop_duplicates().shape[0]
        total_regions = data_frame["region"].drop_duplicates().shape[0]
        summary["proportion_countries"] = summary["n_countries"] / total_countries
        summary["proportion_regions"] = summary["n_regions"] / total_regions

        summary["representative_strain_accession"] = summary["representative_strain"].map(
            strain_to_accession
        )

        summary["representative_strain_ha_sequence"] = summary["representative_strain"].map(
            alignment_by_strain
        )

        summary["days_since_latest_sequence"] = (pd.to_datetime(datetime.date.today()) - pd.to_datetime(summary["latest_sequence"])).dt.days

        return summary
    return (summarize_haplotypes,)


@app.cell
def _(
    alignment_by_lineage_by_strain,
    data_frame_by_lineage,
    np,
    paths_by_lineage,
    pd,
    pivot_by_lineage,
    summarize_haplotypes,
):
    summary_by_lineage = {}
    for _lineage, _lineage_paths in paths_by_lineage.items():    
        lineage_summary = summarize_haplotypes(
            data_frame_by_lineage[_lineage],
            alignment_by_lineage_by_strain[_lineage],
        )

        # Load GAs for lineage.
        ga = pd.read_csv(_lineage_paths["ga"], sep="\t")
        hier_ga = ga.loc[ga["location"] == "hierarchical", ["variant", "median"]].rename(
            columns={
                "variant": "derived_haplotype",
                "median": "median_ga",
            }
        )

        hier_ga = pd.concat([
            hier_ga,
            pd.DataFrame([{
                "derived_haplotype": pivot_by_lineage[_lineage],
                "median_ga": 1.0,
            }])
        ])
        hier_ga["derived_haplotype"] = hier_ga["derived_haplotype"].str.replace("-", ",")

        # Calculate the weighted population average GA.
        hier_ga = hier_ga.merge(
            lineage_summary.reset_index().loc[:, ["derived_haplotype", "count"]],
            on="derived_haplotype",
        )
        hier_ga["frequency"] = hier_ga["count"] / hier_ga["count"].sum()
        weighted_mean_population_ga = (hier_ga["frequency"] * hier_ga["median_ga"]).sum()
        print(weighted_mean_population_ga)

        hier_ga["relative_fitness"] = np.log(hier_ga["median_ga"]) - np.log(weighted_mean_population_ga)
        hier_ga = hier_ga.drop(columns=["count"])

        # Annotate summary by median hierarchical GA.
        lineage_summary = lineage_summary.merge(
            hier_ga,
            how="left",
            on="derived_haplotype",
        )
        lineage_summary["median_ga"] = lineage_summary["median_ga"].fillna(0.0)
        lineage_summary["relative_fitness"] = lineage_summary["relative_fitness"].fillna(0.0)

        # Score haplotypes following approach for clade assignment in:
        # https://github.com/influenza-clade-nomenclature/clade-suggestion-algorithm
        lineage_summary["score_mutations"] = (
            (lineage_summary["mean_ha1_mutations"] * 1) +
            (lineage_summary["mean_epitope_mutations"] * 2) +
            (lineage_summary["mean_recurrent_mutations"] * 2) +
            (lineage_summary["mean_rbs_mutations"] * 3)
        )
        lineage_summary["score_mutations"] = lineage_summary["score_mutations"] / (lineage_summary["score_mutations"].median() + lineage_summary["score_mutations"])
        lineage_summary["score_fitness"] = lineage_summary["relative_fitness"] / lineage_summary["relative_fitness"].abs().max()

        lineage_summary["score_total"] = lineage_summary["score_mutations"] + lineage_summary["score_fitness"] + lineage_summary["proportion_countries"]

        summary_by_lineage[_lineage] = lineage_summary.sort_values(
            "score_total",
            ascending=False,
        )
        summary_by_lineage[_lineage].to_csv(
            _lineage_paths["output"],
            sep="\t",
            index=False,
        )
    return (summary_by_lineage,)


@app.cell
def _(summary_by_lineage):
    summary_by_lineage["h1n1pdm"]
    return


@app.cell
def _(summary_by_lineage):
    summary_by_lineage["h3n2"]
    return


@app.cell
def _(min_sequences, summary_by_lineage):
    ((~summary_by_lineage["h3n2"]["derived_haplotype"].str.startswith("J.2")) & (summary_by_lineage["h3n2"]["count"] >= min_sequences)).sum()
    return


@app.cell
def _(data_frame_by_lineage):
    data_frame_by_lineage["h1n1pdm"].loc[:, ["strain", "date", "region", "derived_haplotype"]].to_csv(
        "h1n1pdm_metadata_with_nextclade.tsv",
        sep="\t",
        index=False,
    )
    return


@app.cell
def _(data_frame_by_lineage):
    data_frame_by_lineage["h3n2"].loc[:, ["strain", "date", "region", "derived_haplotype"]].to_csv(
        "h3n2_metadata_with_nextclade.tsv",
        sep="\t",
        index=False,
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
