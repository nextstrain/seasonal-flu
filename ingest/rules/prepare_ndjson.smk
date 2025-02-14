"""
This part of the workflow handles preparing the GISAID NDJSON for the
curation pipeline.

INPUTS:

    metadata    = data/gisaid_epiflu_isolates.xls
    sequences   = data/gisaid_epiflu_sequence.fasta

OUTPUTS:

    ndjson      = data/gisaid.ndjson

"""

# TODO: Semi-automated workflow
# 1. rule to download the gisaid.ndjson.zst file from S3
# 2. rule to download raw GISAID download files that have not been processed from S3 (gisaid/downloads/unprocessed/*).
# 3. rule with wildcards to link multiple GISAID Excel/FASTA pairs -> NDJSON
# 4. rule to concat all NDJSON, latest downloads _first_
# 5. rule to dedup concatenated gisaid.ndjson by GISAID EPI ISL (Isolate_Id),
#    keep the first record because the latest downloads are _first_


# ------------- For local single GISAID download only ------------------- #

rule link_gisaid_metadata_and_fasta:
    input:
        metadata="data/gisaid_epiflu_isolates.xls",
        sequences="data/gisaid_epiflu_sequence.fasta",
    output:
        ndjson="data/gisaid.ndjson",
    log: "logs/link_gisaid_metadata_and_fasta.txt"
    shell:
        r"""
        ./scripts/link-gisaid-metadata-and-fasta \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            > {output.ndjson:q} \
            2> {log}
        """
