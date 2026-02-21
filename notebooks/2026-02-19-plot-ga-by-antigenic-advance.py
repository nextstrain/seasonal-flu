import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    import altair as alt
    import json
    import pandas as pd
    return alt, json, pd


@app.cell
def _():
    lineage = "h3n2"
    return (lineage,)


@app.cell
def _():
    #human_advance_attribute = "kikawa_2025_cTiterSub"
    #human_advance_attribute = "vcm_SH_2025_cTiterSub"
    human_advance_attribute = "kikawa_2025_2026_cTiterSub"
    return (human_advance_attribute,)


@app.cell
def _(lineage):
    if lineage == "h1n1pdm":
        pivot = "D.3.1.1"
    elif lineage == "h3n2":
        pivot = "K"
    elif lineage == "vic":
        pivot = "C.5"
    return (pivot,)


@app.cell
def _(lineage):
    ga_path = f"builds/full-{lineage}/ha/mlr/fitnesses.tsv"
    return (ga_path,)


@app.cell
def _(lineage):
    antigenic_advance_path = f"full-{lineage}_ha.tsv"
    return (antigenic_advance_path,)


@app.cell
def _(lineage):
    auspice_config_path = f"config/{lineage}/ha/auspice_config.json"
    return (auspice_config_path,)


@app.cell
def _(lineage):
    output_path = f"{lineage}_ga_by_antigenic_advance.png"
    return (output_path,)


@app.cell
def _(auspice_config_path, json):
    with open(auspice_config_path, "r", encoding="utf-8") as fh:
        colorings = json.load(fh).get("colorings")
    return (colorings,)


@app.cell
def _(colorings):
    haplotype_colorings = [dict(coloring["scale"]) for coloring in colorings if coloring["key"] == "emerging_haplotype"][0]
    return (haplotype_colorings,)


@app.cell
def _(haplotype_colorings):
    haplotype_colorings_domain = list(haplotype_colorings.keys())
    haplotype_colorings_range = list(haplotype_colorings.values())
    return haplotype_colorings_domain, haplotype_colorings_range


@app.cell
def _(haplotype_colorings_domain):
    haplotype_colorings_domain
    return


@app.cell
def _(haplotype_colorings_range):
    haplotype_colorings_range
    return


@app.cell
def _(ga_path, pd, pivot):
    ga_without_pivot = pd.read_csv(
        ga_path,
        sep="\t",
    ).query(
        "location == 'hierarchical'"
    ).rename(
        columns={
            "variant": "emerging_haplotype",
            "median": "ga_median",
            "HDI_95_lower": "ga_HDI_95_lower",
            "HDI_95_upper": "ga_HDI_95_upper",
        }
    )

    ga_without_pivot.loc[ga_without_pivot["emerging_haplotype"] == "other", "emerging_haplotype"] = "unassigned"

    ga = pd.concat([
        ga_without_pivot,
        pd.DataFrame([
            {
                "location": "hierarchical",
                "emerging_haplotype": pivot,
                "ga_median": 1.0,
                "ga_HDI_95_lower": 1.0,
                "ga_HDI_95_upper": 1.0,
            }
        ]),
    ])
    return (ga,)


@app.cell
def _(ga):
    ga
    return


@app.cell
def _(antigenic_advance_path, pd):
    antigenic_advance = pd.read_csv(
        antigenic_advance_path,
        sep="\t",
    )
    return (antigenic_advance,)


@app.cell
def _(antigenic_advance):
    antigenic_advance
    return


@app.cell
def _(antigenic_advance, human_advance_attribute):
    antigenic_advance_per_haplotype = antigenic_advance.groupby("emerging_haplotype").aggregate(
        human_advance_median=(human_advance_attribute, "median"),
        human_advance_std=(human_advance_attribute, "std"),
        ferret_advance_median=("cell_hi_cTiterSub", "median"),
        ferret_advance_std=("cell_hi_cTiterSub", "std"),
    ).reset_index()

    antigenic_advance_per_haplotype["human_advance_lower"] = (
        antigenic_advance_per_haplotype["human_advance_median"] - antigenic_advance_per_haplotype["human_advance_std"]
    )
    antigenic_advance_per_haplotype["human_advance_upper"] = (
        antigenic_advance_per_haplotype["human_advance_median"] + antigenic_advance_per_haplotype["human_advance_std"]
    )

    antigenic_advance_per_haplotype["ferret_advance_lower"] = (
        antigenic_advance_per_haplotype["ferret_advance_median"] - antigenic_advance_per_haplotype["ferret_advance_std"]
    )
    antigenic_advance_per_haplotype["ferret_advance_upper"] = (
        antigenic_advance_per_haplotype["ferret_advance_median"] + antigenic_advance_per_haplotype["ferret_advance_std"]
    )
    return (antigenic_advance_per_haplotype,)


@app.cell
def _(antigenic_advance_per_haplotype):
    antigenic_advance_per_haplotype
    return


@app.cell
def _(antigenic_advance_per_haplotype, ga):
    ga_by_antigenic_advance = ga.merge(
        antigenic_advance_per_haplotype,
        on="emerging_haplotype",
    ).round(3)
    return (ga_by_antigenic_advance,)


@app.cell
def _(ga_by_antigenic_advance):
    ga_by_antigenic_advance
    return


@app.cell
def _(
    alt,
    ga_by_antigenic_advance,
    haplotype_colorings_domain,
    haplotype_colorings_range,
):
    base = alt.Chart(ga_by_antigenic_advance)
    tooltip = [
        "emerging_haplotype:N",
        "ga_HDI_95_lower:Q",
        "ga_median:Q",
        "ga_HDI_95_upper:Q",
        "human_advance_lower:Q",
        "human_advance_median:Q",
        "human_advance_upper:Q",
        "ferret_advance_lower:Q",
        "ferret_advance_median:Q",
        "ferret_advance_upper:Q",
    ]

    human_points = base.mark_circle(size=100).encode(
        x=alt.X("human_advance_median:Q", title="median antigenic advance by Kikawa et al. 2025 human sera").scale(zero=False),
        y=alt.Y("ga_median:Q", title="growth advantage").scale(zero=False),
        color=alt.Color("emerging_haplotype:N", title="emerging haplotype").scale(
            domain=haplotype_colorings_domain,
            range=haplotype_colorings_range,
        ),
        tooltip=tooltip,
    )

    human_text = base.mark_text(size=10).encode(
        x=alt.X("human_advance_median:Q", title="median antigenic advance by Kikawa et al. 2025 human sera").scale(zero=False),
        y=alt.Y("ga_median_text:Q", title="growth advantage").scale(zero=False),
        text=alt.Text("emerging_haplotype:N"),
        color=alt.Color("emerging_haplotype:N", title="emerging haplotype").scale(
            domain=haplotype_colorings_domain,
            range=haplotype_colorings_range,
        ),
        tooltip=tooltip,
    ).transform_calculate(
        ga_median_text="datum.ga_median + 0.01",
    )

    human_ga_errors = base.mark_rule().encode(
        x=alt.X("human_advance_median:Q"),
        y=alt.Y("ga_HDI_95_lower:Q"),
        y2=alt.Y2("ga_HDI_95_upper:Q"),
        color=alt.Color("emerging_haplotype"),
        tooltip=tooltip,
    )

    human_advance_errors = base.mark_rule().encode(
        x=alt.X("human_advance_lower:Q"),
        x2=alt.X2("human_advance_upper:Q"),
        y=alt.Y("ga_median:Q"),
        color=alt.Color("emerging_haplotype"),
        tooltip=tooltip,
    )

    human_plot = (human_points + human_ga_errors + human_advance_errors + human_text).properties(
        width=350,
        height=350,
    )

    ferret_points = base.mark_circle(size=100).encode(
        x=alt.X("ferret_advance_median:Q", title="median antigenic advance by ferret cell HIs").scale(zero=False),
        y=alt.Y("ga_median:Q", title="growth advantage").scale(zero=False),
        color=alt.Color("emerging_haplotype:N", title="emerging haplotype").scale(
            domain=haplotype_colorings_domain,
            range=haplotype_colorings_range,
        ),
        tooltip=tooltip,
    )

    ferret_text = base.mark_text(size=10).encode(
        x=alt.X("ferret_advance_median:Q", title="median antigenic advance by ferret cell HIs").scale(zero=False),
        y=alt.Y("ga_median_text:Q", title="growth advantage").scale(zero=False),
        text=alt.Text("emerging_haplotype:N"),
        color=alt.Color("emerging_haplotype:N", title="emerging haplotype").scale(
            domain=haplotype_colorings_domain,
            range=haplotype_colorings_range,
        ),
        tooltip=tooltip,
    ).transform_calculate(
        ga_median_text="datum.ga_median + 0.01",
    )

    ferret_ga_errors = base.mark_rule().encode(
        x=alt.X("ferret_advance_median:Q"),
        y=alt.Y("ga_HDI_95_lower:Q"),
        y2=alt.Y2("ga_HDI_95_upper:Q"),
        color=alt.Color("emerging_haplotype"),
        tooltip=tooltip,
    )

    ferret_advance_errors = base.mark_rule().encode(
        x=alt.X("ferret_advance_lower:Q"),
        x2=alt.X2("ferret_advance_upper:Q"),
        y=alt.Y("ga_median:Q"),
        color=alt.Color("emerging_haplotype"),
        tooltip=tooltip,
    )

    ferret_plot = (ferret_points + ferret_ga_errors + ferret_advance_errors + ferret_text).properties(
        width=350,
        height=350,
    )

    plot = human_plot | ferret_plot

    plot = plot.configure_axis(
        grid=False,
        labelFontSize=12,
        titleFontSize=12,
    ).configure_legend(
        disable=True,
    ).configure_view(
        stroke=None,
    )

    plot
    return (plot,)


@app.cell
def _(output_path, plot):
    plot.save(
        output_path,
        ppi=200,
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
