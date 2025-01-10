'''
This file contains rules to interact with the fauna titer data base.
It produces files in the directory `data` and requires no input files.
the endpoints are

sequences = "data/{lineage}/{segment}.fasta"
metadata = "data/{lineage}/metadata.tsv"
titers = "data/{lineage}/{center}_{passage}_{assay}_titers.tsv"

'''

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

rule download_sequences:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/{lineage}/raw_{segment}.fasta"
    params:
        fasta_fields = config["fauna_fasta_fields"],
    resources:
        concurrent_fauna = 1
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/download_sequences_{lineage}_{segment}.txt"
    log:
        "logs/download_sequences_{lineage}_{segment}.txt"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus flu \
            --fasta_fields {params.fasta_fields} \
            --resolve_method split_passage \
            --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
            --path data \
            --fstem {wildcards.lineage}/raw_{wildcards.segment} 2>&1 | tee {log}
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

def get_host_query(wildcards):
    query_by_host = {
        "ferret": "--istr-not-in-fld serum_id:mouse --istr-not-in-fld serum_id:human",
        "human": "--istr-in-fld serum_id:human",
        "mouse": "--istr-in-fld serum_id:mouse",
    }
    return query_by_host[wildcards.host]

rule select_titers_by_host:
    input:
        titers = "data/{lineage}/{center}_{passage}_{assay}_titers.tsv",
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
