import datetime
import pandas as pd
from treetime.utils import numeric_date
from snakemake.utils import min_version


wildcard_constraints:
    lineage = r'h1n1pdm|h3n2|vic|yam|h1n1',
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

# Nextclade datasets that provide the reference sequence and genome annotation of the
# builds. Where there are two datasets for a segment, the one whose reference has the
# length of the reference of the nomenclature repository, so that positions match the
# nomenclature: A/Wisconsin/67/2005 for H3N2 HA rather than A/Darwin/6/2021; for H1N1pdm
# HA both have that length and A/California/7/2009 is the reference of the nomenclature.
# A build can name another dataset with `nextclade_dataset`, or give a `reference` and
# `annotation` for a segment without a dataset.
INTERNAL_SEGMENTS = ["pb2", "pb1", "pa", "np", "mp", "ns"]
nextclade_dataset_by_lineage_and_segment = {
    "h3n2": {
        "ha": "nextstrain/flu/h3n2/ha/CY163680",
        "na": "nextstrain/flu/h3n2/na/EPI1857215",
        **{segment: f"nextstrain/flu/h3n2/{segment}" for segment in INTERNAL_SEGMENTS},
    },
    "h1n1pdm": {
        "ha": "nextstrain/flu/h1n1pdm/ha/CY121680",
        "na": "nextstrain/flu/h1n1pdm/na/MW626056",
        **{segment: f"nextstrain/flu/h1n1pdm/{segment}" for segment in INTERNAL_SEGMENTS},
    },
    "vic": {
        "ha": "nextstrain/flu/vic/ha/KX058884",
        "na": "nextstrain/flu/vic/na/CY073894",
        **{segment: f"nextstrain/flu/vic/{segment}" for segment in INTERNAL_SEGMENTS},
    },
    # the other segments of yam from the datasets of all influenza B
    "yam": {
        "ha": "nextstrain/flu/yam/ha/JN993010",
        "na": "nextstrain/flu/b/na/CY073894",
        **{segment: f"nextstrain/flu/b/{segment}" for segment in INTERNAL_SEGMENTS},
    },
    "h1n1": {
        segment: f"nextstrain/flu/h1n1/{segment}" for segment in ["ha", "na", *INTERNAL_SEGMENTS]
    },
}

lcr_url_by_lineage_and_segment = {
    "h1n1pdm": {
        "ha": {
            "weights": "https://raw.githubusercontent.com/SamT123/nextstrain-flu-convergence/main/results/H1pdm09/lcr_weights.tsv",
            "windows": "https://raw.githubusercontent.com/SamT123/nextstrain-flu-convergence/main/results/H1pdm09/lcr_windows.tsv",
        },
    },
    "h3n2": {
        "ha": {
            "weights": "https://raw.githubusercontent.com/SamT123/nextstrain-flu-convergence/main/results/H3/lcr_weights.tsv",
            "windows": "https://raw.githubusercontent.com/SamT123/nextstrain-flu-convergence/main/results/H3/lcr_windows.tsv",
        },
    },
    "vic": {
        "ha": {
            "weights": "https://raw.githubusercontent.com/SamT123/nextstrain-flu-convergence/main/results/BVic/lcr_weights.tsv",
            "windows": "https://raw.githubusercontent.com/SamT123/nextstrain-flu-convergence/main/results/BVic/lcr_windows.tsv",
        },
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
