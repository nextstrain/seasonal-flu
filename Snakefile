import datetime
import pandas as pd
from treetime.utils import numeric_date


wildcard_constraints:
    lineage = r'h1n1pdm|h3n2|vic|yam',
    segment = r'pb2|pb1|pa|ha|np|na|mp|ns',
    center = r'who|cdc|crick|niid|crick|vidrl',
    passage = r'cell|egg',
    assay = r'fra|hi',
    host = r'ferret|human|mouse'

# Load distance map configuration for lineages and segments.
distance_map_config = pd.read_table("config/distance_maps.tsv")

lineage_name_by_abbreviation = {
    "h3n2": "H3N2",
    "h1n1pdm": "H1N1pdm",
    "vic": "Vic",
    "yam": "Yam",
}

root_lineage = {
    "h3n2": {
        "ha":{
            "6m": "J.2",
            "1y": "J.2",
            "2y": "J",
            "3y": "G",
            "6y": "B",
            "12y": "A"
        },
        "na":{
            "6m": "B.4",
            "1y": "B.4",
            "2y": "B.4",
            "3y": "B.2",
            "6y": "A.2",
            "12y": "A"
        }
    },
    "h1n1pdm": {
        "ha":{
            "6m": "C.1",
            "1y": "C.1",
            "2y": "C.1",
            "3y": "C.1",
            "6y": "C",
            "12y": "A"
        },
        "na":{
            "6m": "C.5.3",
            "1y": "C.5.3",
            "2y": "C.5.3",
            "3y": "C",
            "6y": "A",
            "12y": "A"
        }
    },
    "vic": {
        "ha":{
            "6m": "C",
            "1y": "C",
            "2y": "C",
            "3y": "A.3",
            "6y": "A",
            "12y": "A"
        },
        "na":{
            "6m": "B",
            "1y": "B",
            "2y": "B",
            "3y": "A.1",
            "6y": "A",
            "12y": "A"
        }
    },
}

# Repositories of the hierarchical lineage nomenclature, which provide the clade and
# subclade definitions, the subclade guide tree and the lineage founder sequences.
nomenclature_repo_by_lineage_and_segment = {
    "h1n1pdm": {
        "ha": "seasonal_A-H1N1pdm_HA",
        "na": "seasonal_A-H1N1pdm_NA",
    },
    "h3n2": {
        "ha": "seasonal_A-H3N2_HA",
        "na": "seasonal_A-H3N2_NA",
    },
    "vic": {
        "ha": "seasonal_B-Vic_HA",
        "na": "seasonal_B-Vic_NA",
    }
}

# The founders are only on the `founders` branch until it is merged into main.
nomenclature_branch = config.get("nomenclature", {}).get("branch", "founders")

def nomenclature_url(lineage, segment, filename):
    """URL of a generated file of the nomenclature, or None if the segment has none."""
    repo = nomenclature_repo_by_lineage_and_segment.get(lineage, {}).get(segment)
    if repo is None:
        return None
    return f"https://raw.githubusercontent.com/influenza-clade-nomenclature/{repo}/{nomenclature_branch}/.auto-generated/{filename}"

include: "workflow/snakemake_rules/common.smk"

if "data_source" in config and config["data_source"]=='fauna':
    include: "workflow/snakemake_rules/download_from_fauna.smk"

include: "workflow/snakemake_rules/array_builds.smk"

include: "workflow/snakemake_rules/select_strains.smk"

include:  "workflow/snakemake_rules/core.smk"

include:  "workflow/snakemake_rules/export.smk"

include:  "workflow/snakemake_rules/titer_models.smk"

include:  "workflow/snakemake_rules/fitness.smk"

def _get_build_outputs():
    outputs = []
    for build_name, build_params in config["builds"].items():
        for segment in config["segments"]:
            outputs.extend([
                f"auspice/{build_name}_{segment}.json",
                f"auspice/{build_name}_{segment}_tip-frequencies.json",
            ])

            if segment == "ha" and build_params.get("enable_measurements", False):
                outputs.append(f"auspice/{build_name}_{segment}_measurements.json")

        if build_params.get("enable_forecasts", False):
            for model in config["fitness_model"]["models"]:
               outputs.append(f"auspice/{build_name}_ha_{model}_forecast-tip-frequencies.json")

    return outputs

if "custom_rules" in config:
    for rule_file in config["custom_rules"]:
        include: rule_file

rule all:
    input:
        _get_build_outputs()
