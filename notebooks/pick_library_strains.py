import marimo

__generated_with = "0.14.17"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    import json
    import pandas as pd
    return json, mo, pd


@app.cell
def _(mo):
    lineage_picker = mo.ui.radio(options=["h1n1pdm", "h3n2"])
    lineage_picker
    return (lineage_picker,)


@app.cell
def _(lineage_picker):
    lineage = lineage_picker.value
    return (lineage,)


@app.cell
def _(lineage):
    if lineage == "h1n1pdm":
        distance_map_path = "config/distance_maps/h1n1pdm/ha/canton.json"
    elif lineage == "h3n2":
        distance_map_path = "config/distance_maps/h3n2/ha/wolf.json"
    return (distance_map_path,)


@app.cell
def _(distance_map_path, json):
    with open(distance_map_path, "r", encoding="utf-8") as _fh:
        distance_map = json.load(_fh)

    epitope_sites = set(distance_map["map"]["HA1"].keys())
    return (epitope_sites,)


@app.cell
def _(json, lineage):
    rbs_sites = set()
    if lineage == "h3n2":
        rbs_distance_map_path = "config/distance_maps/h3n2/ha/koel.json"
        with open(rbs_distance_map_path, "r", encoding="utf-8") as _fh:
            rbs_distance_map = json.load(_fh)

        rbs_sites = set(rbs_distance_map["map"]["HA1"].keys())
    return (rbs_sites,)


@app.cell
def _(epitope_sites, lineage, pd, rbs_sites):
    df = pd.read_csv(
        f"builds/{lineage}/metadata.tsv",
        sep="\t",
        dtype="str",
        na_filter=False,
    )

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
    return (df,)


@app.cell
def _(df):
    distinct_haplotypes = (df["derived_haplotype"].value_counts() > 0).sum()
    return (distinct_haplotypes,)


@app.cell
def _(df):
    non_singleton_distinct_haplotypes = (df["derived_haplotype"].value_counts() > 1).sum()
    return (non_singleton_distinct_haplotypes,)


@app.cell
def _(distinct_haplotypes, lineage, mo, non_singleton_distinct_haplotypes):
    mo.md(rf"""Found {distinct_haplotypes} distinct haplotypes for {lineage} including {non_singleton_distinct_haplotypes} with more than 1 sequence.""")
    return


@app.cell
def _(df):
    strain_to_accession = dict(df.loc[:, ["strain", "gisaid_epi_isl"]].values)
    return (strain_to_accession,)


@app.cell
def _(df, strain_to_accession):
    summary = df.groupby("derived_haplotype").aggregate(
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
    return (summary,)


@app.cell
def _(summary):
    summary
    return


@app.cell
def _(lineage, summary):
    summary.to_csv(
        f"{lineage}_haplotypes.tsv",
        sep="\t",
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
