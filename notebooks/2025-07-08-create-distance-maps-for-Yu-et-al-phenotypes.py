import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Create distance maps for Yu et al. H3N2 phenotypes

    Tim Yu's project measured [multiple H3N2 phenotypes including serum escape, cell entry, and pH stability](https://github.com/dms-vep/Flu_H3_Massachusetts2022_DMS/blob/main/results/summaries/Phenotypes.csv).
    This notebook converts these phenotype to [Augur distance maps](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/distance.html) which allows us to calculate cumulative phenotypic effects for each node in a phylogeny based on the mutations present on all branches leading to each node.
    Because the phenotypes reflect the effect of individual amino acid substitutions from a specific wildtype allele, the distance maps define effects relative to the wildtype allele at each position as shown in the following example:

    ```json
    {
        "default": 0.0,
        "map": {
           "HA1": {
               "112": [
                   {
                       "from": "V",
                       "to": "I",
                       "weight": 1.192
                   },
                   {
                       "from": "V",
                       "to": "F",
                       "weight": 0.002
                   }
               ]
           }
       }
    }
    ```
    """
    )
    return


@app.cell
def _():
    import marimo as mo

    import json
    import pandas as pd
    return json, mo, pd


@app.cell
def _(pd):
    phenotypes = pd.read_csv("https://github.com/dms-vep/Flu_H3_Massachusetts2022_DMS/raw/refs/heads/main/results/summaries/Phenotypes.csv")
    return (phenotypes,)


@app.cell
def _(phenotypes):
    phenotypes
    return


@app.cell
def _(phenotypes):
    ha1_max_site = phenotypes.loc[phenotypes["region"] == "HA1", "site"].max()
    return (ha1_max_site,)


@app.cell
def _(phenotypes):
    phenotypes.loc[phenotypes["region"] == "HA2", "site"].min()
    return


@app.cell
def _(phenotypes):
    phenotypes.head()
    return


@app.cell
def _(ha1_max_site, pd, phenotypes):
    cell_entry_map = {
        "name": "Yu et al. 2025 cell entry",
        "default": 0.0,
        "map": {
            "HA1": {},
            "HA2": {},
        }
    }

    ph_stability_map = {
        "name": "Yu et al. 2025 pH stability",
        "default": 0.0,
        "map": {
            "HA1": {},
            "HA2": {},
        }
    }

    for record_index, record in phenotypes.iterrows():
        if record["wildtype"] == record["mutant"]:
            continue

        site = record["site"]
        # Adjust coordinates of HA2 sites to start at 1 relative to the end of HA1.
        if record["region"] == "HA2":
            site -= ha1_max_site
            region = "HA2"
        else:
            region = "HA1"

        site = str(site)

        if not pd.isnull(record["MDCKSIAT1 cell entry"]):
            if site not in cell_entry_map["map"][region]:
                cell_entry_map["map"][region][site] = []

            cell_entry_map["map"][region][str(site)].append({
                "from": record["wildtype"],
                "to": record["mutant"],
                "weight": record["MDCKSIAT1 cell entry"],
            })

        if not pd.isnull(record["pH stability"]):
            if site not in ph_stability_map["map"][region]:
                ph_stability_map["map"][region][site] = []

            ph_stability_map["map"][region][str(site)].append({
                "from": record["wildtype"],
                "to": record["mutant"],
                "weight": record["pH stability"],
            })
    return cell_entry_map, ph_stability_map


@app.cell
def _(cell_entry_map, json):
    with open("Yu_et_al_2025_cell_entry.json", "w") as _oh:
        json.dump(cell_entry_map, _oh)
    return


@app.cell
def _(json, ph_stability_map):
    with open("Yu_et_al_2025_ph_stability.json", "w") as _oh:
        json.dump(ph_stability_map, _oh)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
