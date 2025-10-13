import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    import json
    import pandas as pd
    return json, mo, pd


@app.cell(hide_code=True)
def _(
    haplotypes_by_lineage,
    min_date,
    min_sequences,
    mo,
    nonsingleton_haplotypes_by_lineage,
):
    mo.md(
        rf"""
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
    """
    )
    return


@app.cell
def _():
    min_date = "2025-09-01"
    return (min_date,)


@app.cell
def _(mo):
    min_sequences_picker = mo.ui.number(
        start=1,
        stop=20,
        value=2,
        label="min sequences",
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
            "output": "h1n1pdm_haplotypes.tsv",
        },
        "h3n2": {
            "distance_map": "config/distance_maps/h3n2/ha/wolf.json",
            "rbs_distance_map": "config/distance_maps/h3n2/ha/koel.json",
            "metadata": "builds/h3n2/metadata.tsv",
            "output": "h3n2_haplotypes.tsv",
        },
    }
    return (paths_by_lineage,)


@app.cell
def _(json):
    def load_distance_map(distance_map_path):
        with open(distance_map_path, "r", encoding="utf-8") as fh:
            distance_map = json.load(fh)

        return set(distance_map["map"]["HA1"].keys())
    return (load_distance_map,)


@app.cell
def _(pd):
    def load_metadata(metadata_path, min_date=None, epitope_sites=None, rbs_sites=None):
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

        if min_date:
            df = df[df["date"] >= min_date].copy()

        df["founder_ha1_substitutions"] = df["founderMuts['subclade'].aaSubstitutions"].apply(
            lambda subs: ",".join([
                sub.replace("HA1:", "")
                for sub in subs.split(",")
                if sub.startswith("HA1:")
            ])
        )

        df["derived_haplotype"] = (df["subclade"] + ":" + df["founder_ha1_substitutions"]).str.rstrip(":")

        df["epitope_mutations_from_subclade"] = df["founder_ha1_substitutions"].apply(
            lambda subs: sum(
                sub[1:-1] in epitope_sites
                for sub in subs.split(",")
            )
        )

        df["rbs_mutations_from_subclade"] = df["founder_ha1_substitutions"].apply(
            lambda subs: sum(
                sub[1:-1] in rbs_sites
                for sub in subs.split(",")
            )
        )

        df["glycosylation_count"] = df["glycosylation"].apply(
            lambda glyc: len([site for site in glyc.split(";") if site.startswith("HA1:")])
        )

        return df
    return (load_metadata,)


@app.cell
def _(
    load_distance_map,
    load_metadata,
    min_date,
    min_sequences,
    paths_by_lineage,
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
        )

        haplotypes_by_lineage[lineage] = (data_frame_by_lineage[lineage]["derived_haplotype"].value_counts() > 0).sum()
        nonsingleton_haplotypes_by_lineage[lineage] = (data_frame_by_lineage[lineage]["derived_haplotype"].value_counts() >= min_sequences).sum()
    return (
        data_frame_by_lineage,
        haplotypes_by_lineage,
        lineage_paths,
        nonsingleton_haplotypes_by_lineage,
    )


@app.function
def summarize_haplotypes(data_frame):
    strain_to_accession = dict(data_frame.loc[:, ["strain", "gisaid_epi_isl"]].values)

    summary = data_frame.groupby("derived_haplotype").aggregate(
        count=("derived_haplotype", "count"),
        mean_epitope_mutations=("epitope_mutations_from_subclade", "mean"),
        mean_rbs_mutations=("rbs_mutations_from_subclade", "mean"),
        mean_glycosylation_count=("glycosylation_count", "mean"),
        latest_sequence=("date", "max"),
        representative_strain=("strain", "first"),
    ).sort_values(
        "count",
        ascending=False,
    )

    summary["representative_strain_accession"] = summary["representative_strain"].map(
        strain_to_accession
    )

    return summary


@app.cell
def _(data_frame_by_lineage, lineage_paths, paths_by_lineage):
    summary_by_lineage = {}
    for _lineage, _lineage_paths in paths_by_lineage.items():    
        summary_by_lineage[_lineage] = summarize_haplotypes(
            data_frame_by_lineage[_lineage],
        )

        summary_by_lineage[_lineage].to_csv(
            lineage_paths["output"],
            sep="\t",
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
def _():
    return


if __name__ == "__main__":
    app.run()
