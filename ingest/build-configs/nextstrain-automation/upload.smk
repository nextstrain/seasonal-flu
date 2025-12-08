"""
This part of the workflow handles uploading files to AWS S3.
"""


def all_processed_gisaid_pairs(wildcards):
    """
    Check which unprocessed files were fetched so that we only move the
    `gisaid_pairs` that were processed in this run of the workflow.
    """
    checkpoint_output = checkpoints.fetch_unprocessed_files.get(**wildcards).output[0]
    return expand(
        "results/mv-processed/{gisaid_pair}.done",
        gisaid_pair=glob_wildcards(os.path.join(checkpoint_output, "{gisaid_pair}-metadata.xls.zst")).gisaid_pair
    )


rule upload_all:
    input:
        ndjson="results/upload/gisaid.ndjson.upload",
        metadata=expand("results/upload/{dataset}/metadata.tsv.upload",
                        dataset=list(config['filtering'].keys())),
        sequences=expand("results/upload/{dataset}/{segment}.fasta.upload",
                         dataset=list(config['filtering'].keys()),
                         segment=config["segments"]),
        mv_processed=all_processed_gisaid_pairs,


rule upload_gisaid_ndjson:
    input:
        ndjson="data/gisaid.ndjson",
    output:
        flag="results/upload/gisaid.ndjson.upload",
    params:
        s3_dst=config["s3_dst"],
    shell:
        r"""
        ./vendored/upload-to-s3 \
            --quiet \
            {input.ndjson:q} \
            {params.s3_dst:q}/gisaid.ndjson.zst \
            2>&1 | tee {output.flag:q}
        """


rule mv_processed_gisaid_pair:
    """
    Move the processed gisaid_pair from /unprocessed to /processed on AWS S3
    so that they are not reprocessed in the next run of the automated ingest.

    The records in the processed gisaid_pair should be included in the cached
    gisaid.ndjson, so only move it after the gisaid.ndjson was successfully
    uploaded to S3.
    """
    input:
        ndjson_flag="results/upload/gisaid.ndjson.upload",
    output:
        flag=touch("results/mv-processed/{gisaid_pair}.done")
    params:
        s3_dst=f"{config['s3_dst']}/gisaid-downloads/processed/",
        s3_src=f"{config['s3_src']}/gisaid-downloads/unprocessed/",
    shell:
        r"""
        aws s3 mv \
            {params.s3_src:q} \
            {params.s3_dst:q} \
            --recursive \
            --exclude "*" \
            --include "{wildcards.gisaid_pair}-metadata.xls.zst" \
            --include "{wildcards.gisaid_pair}-sequences.fasta.zst"
        """


rule upload_metadata:
    input:
        metadata="results/{dataset}/metadata.tsv",
    output:
        flag="results/upload/{dataset}/metadata.tsv.upload",
    params:
        s3_dst=config["s3_dst"],
    shell:
        r"""
        ./vendored/upload-to-s3 \
            --quiet \
            {input.metadata:q} \
            {params.s3_dst:q}/{wildcards.dataset}/metadata.tsv.xz \
            2>&1 | tee {output.flag:q}
        """


rule upload_sequences:
    input:
        sequences="results/{dataset}/{segment}.fasta",
    output:
        flag="results/upload/{dataset}/{segment}.fasta.upload",
    params:
        s3_dst=config["s3_dst"],
    shell:
        r"""
        ./vendored/upload-to-s3 \
            --quiet \
            {input.sequences:q} \
            {params.s3_dst:q}/{wildcards.dataset}/{wildcards.segment}/sequences.fasta.xz \
            2>&1 | tee {output.flag:q}
        """
