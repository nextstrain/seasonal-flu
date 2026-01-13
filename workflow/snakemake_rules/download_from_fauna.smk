'''
This file contains rules to interact with the fauna titer data base.
It produces files in the directory `data` and requires no input files.
the endpoints are

sequences = "data/{lineage}/{segment}.fasta"
metadata = "data/{lineage}/metadata.tsv"
titers = "data/{lineage}/{center}_{passage}_{assay}_titers.tsv"

'''

import json

# Limit the number of concurrent fauna connections so that we are less likely
# to overwhelm the rethinkdb server.
# Inspired by the ncov's limit on concurrent deploy jobs
# <https://github.com/nextstrain/ncov/blob/20f5fc3c7032f4575a99745cee3238ecbeebb6e0/workflow/snakemake_rules/export_for_nextstrain.smk#L340-L362>
workflow.global_resources.setdefault("concurrent_fauna", 2)

# fields that will be canonicized by augur parse (upper/lower casing etc)

path_to_fauna = '../fauna'
localrules: download_sequences, download_titers, parse
#
# Define titer data sets to be used.
#
def _get_tdb_databases(wildcards):
    if wildcards.center in ['cdc', 'crick', 'niid', 'vidrl']:
        return wildcards.center + "_tdb tdb"
    else:
        return "cdc_tdb crick_tdb niid_tdb vidrl_tdb tdb"


def _get_tdb_assays(wildcards):
    if wildcards.assay == 'fra':
        return 'fra,hint'
    if wildcards.assay == 'hi':
        return 'hi,hi_oseltamivir'
    return wildcards.assay

def _get_virus_passage_category(wildcards):
    # Exclude titer measurements for egg-passaged test viruses from
    # cell-passaged titer data. This filter prevents mutations associated with
    # egg-passaging from appearing in the titer substitition model for
    # cell-passaged titers. For egg-passaged titer models, we accept all other
    # passage types, too.
    if wildcards.passage == "cell":
        return "virus_passage_category:unpassaged,cell,undetermined,N/A"
    else:
        return ""

def _get_prioritized_seqs_file(wildcards):
    prioritized_seqs_file = []
    for build_name, build_params in config["builds"].items():
        if build_params["lineage"] == wildcards.lineage:
            prioritized_seqs_file = build_params.get('prioritized_seqs_file', prioritized_seqs_file)
            break
    return prioritized_seqs_file

rule download_sequences:
    input:
        prioritized_seqs_file = _get_prioritized_seqs_file,
    output:
        sequences = "data/{lineage}/raw_{segment}.fasta"
    params:
        fasta_fields = config["fauna_fasta_fields"],
        prioritized_seqs_file = lambda wildcards, input:
            f"--prioritized_seqs_file {input.prioritized_seqs_file!r}"
            if input.prioritized_seqs_file
            else ""
    resources:
        concurrent_fauna = 1
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/download_sequences_{lineage}_{segment}.txt"
    log:
        "logs/download_sequences_{lineage}_{segment}.txt"
    shell:
        r"""
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus flu \
            --fasta_fields {params.fasta_fields} \
            --resolve_method split_passage \
            --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
            {params.prioritized_seqs_file} \
            --path data \
            --fstem {wildcards.lineage}/raw_{wildcards.segment} 2>&1 | tee {log}
        """

S3_PATH = config.get("s3_path", "s3://nextstrain-data-private/files/workflows/seasonal-flu")

# This is almost identical to `download_parsed_metadata` in `download_from_s3.smk` however
# we include it here so that fauna (titer) workflows can access it. If we import the
# rules file we get a number of rule name & file name clashes which I (James) don't
# want to resolve! Note that this rule only works (conceptually) since we have switched
# our "source of truth" metadata from fauna to the new curation pipeline.
rule download_curated_metadata:
    output:
        metadata="data/{lineage}/curated-metadata-for-titer-matching.tsv.xz",
    params:
        s3_path=S3_PATH + "/{lineage}/metadata.tsv.xz"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - > {output.metadata}
        """

# The fauna strain map is a TSV linking (fauna) strain name to EPI ISL.
# It was manually created from a fauna download on 2026-01-12, after
# which time no more sequence data will be added.
rule download_fauna_strain_map:
    output:
        metadata="data/{lineage}/fauna-strain-map.tsv.xz",
    params:
        s3_path=S3_PATH + "/{lineage}/fauna-strain-map.tsv.xz"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - > {output.metadata}
        """

rule download_titers:
    output:
        titers = "data/{lineage}/{center}_{passage}_{assay}_titers.tsv"
    params:
        dbs = _get_tdb_databases,
        assays = _get_tdb_assays,
        virus_passage_category=_get_virus_passage_category,
    resources:
        concurrent_fauna = 1
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/download_titers_{lineage}_{center}_{passage}_{assay}.txt"
    log:
        "logs/download_titers_{lineage}_{center}_{passage}_{assay}.txt"
    shell:
        """
        python3 {path_to_fauna}/tdb/download.py \
            --database {params.dbs} \
            --virus flu \
            --subtype {wildcards.lineage} \
            --select assay_type:{params.assays} {params.virus_passage_category} serum_passage_category:{wildcards.passage} \
            --path data \
            --fstem {wildcards.lineage}/{wildcards.center}_{wildcards.passage}_{wildcards.assay} 2>&1 | tee {log}
        """

rule remap_titer_strain_names:
    input:
        titers = "data/{lineage}/{center}_{passage}_{assay}_titers.tsv",
        fauna_strain_map = "data/{lineage}/fauna-strain-map.tsv.xz",
        curated_metadata = "data/{lineage}/curated-metadata-for-titer-matching.tsv.xz",
    output:
        titers = "data/{lineage}/{center}_{passage}_{assay}_titers-strains-remapped.tsv",
        stats = "data/{lineage}/{center}_{passage}_{assay}_matching-stats.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/remap_titer_strain_names_{lineage}_{center}_{passage}_{assay}.txt"
    log:
        "logs/remap_titer_strain_names_{lineage}_{center}_{passage}_{assay}.txt"
    params:
        stats_metadata = lambda w: json.dumps({'lineage': w.lineage, 'center': w.center, 'passage': w.passage, 'assay': w.assay})
    shell:
        """
        scripts/remap-titer-strain-names.py \
            --fauna-strain-map {input.fauna_strain_map} \
            --metadata {input.curated_metadata} \
            --titers {input.titers} \
            --output {output.titers} \
            --stats-metadata {params.stats_metadata:q} \
            --stats {output.stats} 2> {log}
        """

def get_host_query(wildcards):
    query_by_host = {
        "ferret": "--istr-not-in-fld serum_id:mouse --istr-not-in-fld serum_id:human",
        "human": "--istr-in-fld serum_id:human",
        "mouse": "--istr-in-fld serum_id:mouse",
    }
    return query_by_host[wildcards.host]

rule select_titers_by_host:
    input:
        titers = "data/{lineage}/{center}_{passage}_{assay}_titers-strains-remapped.tsv",
    output:
        titers = "data/{lineage}/{center}_{host}_{passage}_{assay}_titers.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/select_titers_by_host_{lineage}_{center}_{host}_{passage}_{assay}.txt"
    log:
        "logs/select_titers_by_host_{lineage}_{center}_{host}_{passage}_{assay}.txt"
    params:
        host_query=get_host_query,
    shell:
        """
        tsv-filter -H {params.host_query} {input.titers} > {output.titers} 2> {log}
        """
