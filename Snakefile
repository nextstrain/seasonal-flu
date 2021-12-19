

if "data_source" in config and config["data_source"]=='fauna':
    include: "workflow/snakemake_rules/download_from_fauna.smk"

if "array-builds" in config:
    include: "workflow/snakemake_rules/array_builds.smk"

include: "workflow/snakemake_rules/select_strains.smk"

include:  "workflow/snakemake_rules/core.smk"

rule all:
    input:
        [f"builds/{b}/ha/tree_raw.nwk" for b in config["builds"]]
