import datetime
import pandas as pd
from treetime.utils import numeric_date
from snakemake.utils import min_version


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

# Optionally support inputs to keep workflow backwards compatible
if config.get("inputs"):
    # Minimum Snakemake version needed for the storage plugins used in remote_files.smk
    min_version("8.0.0")

    from packaging import version
    from augur.__version__ import __version__ as augur_version
    import sys

    # Minimum Augur version needed for simple merge rules
    min_augur_version = "34.0.0"
    if version.parse(augur_version) < version.parse(min_augur_version):
      print("This pipeline needs a newer version of augur than you currently have...")
      print(f"Current augur version: {augur_version}. Minimum required: {min_augur_version}")
      sys.exit(1)

    include: "shared/vendored/snakemake/config.smk"
    include: "shared/vendored/snakemake/remote_files.smk"
    include: "workflow/snakemake_rules/merge_inputs.smk"

include: "workflow/snakemake_rules/common.smk"

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
