'''
This file contains rules to interact with the fauna titer data base.
It produces files in the directory `data` and requires no input files.
the endpoints are

sequences = "data/{lineage}/{segment}.fasta"
metadata = "data/{lineage}/metadata.tsv"
titers = "data/{lineage}/{center}_{passage}_{assay}_titers.tsv"

'''

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

rule download_sequences:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/{lineage}/raw_{segment}.fasta"
    params:
        fasta_fields = config["fauna_fasta_fields"],
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
        assays = _get_tdb_assays
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
            --select assay_type:{params.assays} serum_passage_category:{wildcards.passage} \
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
