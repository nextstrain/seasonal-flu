import datetime

wildcard_constraints:
    segment = r'pb2|pb1|pa|ha|np|na|ma'

if "data_source" in config and config["data_source"]=='fauna':
    include: "workflow/snakemake_rules/download_from_fauna.smk"

if "array-builds" in config:
    include: "workflow/snakemake_rules/array_builds.smk"

include: "workflow/snakemake_rules/select_strains.smk"

include:  "workflow/snakemake_rules/core.smk"

include:  "workflow/snakemake_rules/export.smk"

include:  "workflow/snakemake_rules/titer_models.smk"

include:  "workflow/snakemake_rules/fitness.smk"

rule all:
    input:
        [f"auspice/{b}_ha.json" for b in config[build_dir + ""]] + \
        [f"auspice/{b}_ha_tip-frequencies.json" for b in config[build_dir + ""]] + \
        [f"auspice/{b}_na.json" for b in config[build_dir + ""]] + \
        [f"auspice/{b}_na_tip-frequencies.json" for b in config[build_dir + ""]]
