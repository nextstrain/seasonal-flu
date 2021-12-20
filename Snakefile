

if "data_source" in config and config["data_source"]=='fauna':
    include: "workflow/snakemake_rules/download_from_fauna.smk"

if "array-builds" in config:
    include: "workflow/snakemake_rules/array_builds.smk"

include: "workflow/snakemake_rules/select_strains.smk"

include:  "workflow/snakemake_rules/core.smk"

include:  "workflow/snakemake_rules/export.smk"

rule all:
    input:
        [f"auspice/{b}/ha.json" for b in config[build_dir + ""]] + \
        [f"auspice/{b}/ha_tip-frequencies.json" for b in config[build_dir + ""]] + \
        [f"auspice/{b}/na.json" for b in config[build_dir + ""]] + \
        [f"auspice/{b}/na_tip-frequencies.json" for b in config[build_dir + ""]]
        #        [fbuild_dir + "/{b}/ha/aa_muts.json" for b in config[build_dir + ""]]
