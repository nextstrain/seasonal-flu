"""
This part of the workflow handles preparing the GISAID NDJSON for the
curation pipeline.

INPUTS:

    metadata    = data/<YYYY-MM-DD-N>-metadata.xls
    sequences   = data/<YYYY-MM-DD-N>-sequences.fasta

<YYYY-MM-DD> is the date the files were downloaded from GISAID
<N> is the number of the download since GISAID limits the number of records per download

OUTPUTS:

    ndjson      = data/gisaid.ndjson

"""


rule link_gisaid_metadata_and_fasta:
    input:
        metadata="data/{gisaid_pair}-metadata.xls",
        sequences="data/{gisaid_pair}-sequences.fasta",
    output:
        ndjson="data/{gisaid_pair}.ndjson",
    log: "logs/link_gisaid_metadata_and_fasta/{gisaid_pair}.txt"
    shell:
        r"""
        ./scripts/link-gisaid-metadata-and-fasta \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            > {output.ndjson:q} \
            2> {log}
        """


def aggregate_gisaid_ndjsons(wildcards):
    """
    Input function for rule concatenate_gisaid_ndjsons to check which
    GISAID pairs to include the output.
    """
    if len(config.get("gisaid_pairs", [])):
        GISAID_PAIRS = config["gisaid_pairs"]
    elif config.get('s3_src') and hasattr(checkpoints, "fetch_unprocessed_files"):
        # Use checkpoint for the Nextstrain automation
        checkpoint_output = checkpoints.fetch_unprocessed_files.get(**wildcards).output[0]
        GISAID_PAIRS, = glob_wildcards(os.path.join(checkpoint_output, "{gisaid_pair}-metadata.xls.zst"))
        # Reverse sort to list latest downloads first
        GISAID_PAIRS.sort(reverse=True)
        # Add the GISAID cache last to prioritize the latest downloads
        GISAID_PAIRS.append("gisaid_cache")
    else:
        # Create wildcards for pairs of GISAID downloads
        GISAID_PAIRS, = glob_wildcards("data/{gisaid_pair}-metadata.xls")
        # Reverse sort to list latest downloads first
        GISAID_PAIRS.sort(reverse=True)

    assert len(GISAID_PAIRS), "No GISAID metadata and sequences inputs were found"

    return expand("data/{gisaid_pair}.ndjson", gisaid_pair=GISAID_PAIRS)


rule concatenate_gisaid_ndjsons:
    input:
        ndjsons=aggregate_gisaid_ndjsons,
    output:
        ndjson="data/gisaid.ndjson",
    params:
        gisaid_id_field=config["gisaid_id_field"],
    log: "logs/concatenate_gisaid_ndjsons.txt"
    shell:
        r"""
        (cat {input.ndjsons:q} \
            | ./scripts/dedup-by-gisaid-id \
                --id-field {params.gisaid_id_field:q} \
            > {output.ndjson:q}) 2>> {log:q}
        """
