import datetime
import pandas as pd
from treetime.utils import numeric_date


wildcard_constraints:
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

clade_url_by_lineage_and_segment = {
    "h1n1pdm": {
        "ha": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_A-H1N1pdm_HA/main/.auto-generated/clades.tsv",
    },
    "h3n2": {
        "ha": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/main/.auto-generated/clades.tsv",
    },
    "vic": {
        "ha": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_B-Vic_HA/main/.auto-generated/clades.tsv",
    }
}

subclade_url_by_lineage_and_segment = {
    "h1n1pdm": {
        "ha": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_A-H1N1pdm_HA/main/.auto-generated/subclades.tsv",
        "na": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_A-H1N1pdm_NA/main/.auto-generated/subclades.tsv",
    },
    "h3n2": {
        "ha": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/main/.auto-generated/subclades.tsv",
        "na": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_A-H3N2_NA/main/.auto-generated/subclades.tsv",
    },
    "vic": {
        "ha": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_B-Vic_HA/main/.auto-generated/subclades.tsv",
        "na": "https://raw.githubusercontent.com/influenza-clade-nomenclature/seasonal_B-Vic_NA/main/.auto-generated/subclades.tsv",
    }
}

if "data_source" in config and config["data_source"]=='fauna':
    include: "workflow/snakemake_rules/download_from_fauna.smk"

include: "workflow/snakemake_rules/array_builds.smk"

include: "workflow/snakemake_rules/select_strains.smk"

include:  "workflow/snakemake_rules/core.smk"

include:  "workflow/snakemake_rules/export.smk"

include:  "workflow/snakemake_rules/titer_models.smk"

include:  "workflow/snakemake_rules/fitness.smk"

include:  "workflow/snakemake_rules/report.smk"

def _get_build_outputs():
    outputs = []
    for build_name, build_params in config["builds"].items():
        for segment in config["segments"]:
            outputs.extend([
                f"auspice/{build_name}_{segment}.json",
                f"auspice/{build_name}_{segment}_root-sequence.json",
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
