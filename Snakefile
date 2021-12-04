

if "data_source" in config and config["data_source"]=='fauna':
    include: "workflow/snakemake_rules/download_from_fauna.smk"

