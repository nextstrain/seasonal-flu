'''
This file contains rules to interact with the fauna titer data base.
It produces files in the directory `data` and requires no input files.
the endpoints are

sequences = "data/{lineage}/{segment}.fasta"
metadata = "data/{lineage}/metadata.tsv"
titers = "data/{lineage}/{center}_{passage}_{assay}_titers.tsv"

'''

fasta_fields = ["strain", "virus", "accession", "collection_date", "virus_inclusion_date",
                "region",  "country", "division", "location", "passage_category",
                "originating_lab", "submitting_lab", "age", "gender"]

output_fasta_fields = ["strain", "virus", "accession", "date", "virus_inclusion_date",
                "region",  "country", "division", "location", "passage_category",
                "originating_lab", "submitting_lab", "age", "gender"]

# fields that will be canonicized by augur parse (upper/lower casing etc)
prettify_fields = ["region","country","division","location","originating_lab","submitting_lab"]

path_to_fauna = '../fauna'
localrules: download_sequences, download_titers, parse, metadata
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
        fasta_fields = " ".join(fasta_fields)
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus flu \
            --fasta_fields {params.fasta_fields} \
            --resolve_method split_passage \
            --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
            --path data \
            --fstem {wildcards.lineage}/raw_{wildcards.segment}
        """

rule download_titers:
    message: "Downloading titers from fauna: {wildcards.lineage}, {wildcards.assay}, {wildcards.center}"
    output:
        titers = "data/{lineage}/{center}_{passage}_{assay}_titers.tsv"
    params:
        dbs = _get_tdb_databases,
        assays = _get_tdb_assays
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 {path_to_fauna}/tdb/download.py \
            --database {params.dbs} \
            --virus flu \
            --subtype {wildcards.lineage} \
            --select assay_type:{params.assays} serum_passage_category:{wildcards.passage} \
            --path data \
            --fstem {wildcards.lineage}/{wildcards.center}_{wildcards.passage}_{wildcards.assay}
        """

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.download_sequences.output.sequences
    output:
        sequences = "data/{lineage}/{segment}.fasta",
        metadata = "data/{lineage}/metadata_{segment}.tsv"
    params:
        fasta_fields =  " ".join(output_fasta_fields),
        prettify_fields = " ".join(prettify_fields)
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields} \
            --prettify-fields {params.prettify_fields}
        """

rule join_metadata:
    input:
        segment_metadata=lambda w: [f"data/{w.lineage}/metadata_{segment}.tsv" for segment in config['segments']],
    output:
        metadata="data/{lineage}/metadata.tsv",
    conda: "../envs/nextstrain.yaml"
    params:
        segments=lambda w: config['segments'],
        segment_columns=["accession"],
        how="outer",
    shell:
        """
        python3 scripts/join_metadata.py \
            --metadata {input.segment_metadata:q} \
            --segments {params.segments:q} \
            --segment-columns {params.segment_columns:q} \
            --how {params.how:q} \
            --output {output.metadata:q}
        """
