"""
This part of the workflow handles uploading files to AWS S3.
"""


rule upload_all:
    input:
        ndjson="results/upload/gisaid.ndjson.upload",
        metadata=expand("results/upload/{lineage}/metadata.tsv.upload",
                        lineage=config["lineages"]),
        sequences=expand("results/upload/{lineage}/{segment}.fasta.upload",
                         lineage=config["lineages"],
                         segment=config["segments"])


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


rule upload_metadata:
    input:
        metadata="results/{lineage}/metadata.tsv",
    output:
        flag="results/upload/{lineage}/metadata.tsv.upload",
    params:
        s3_dst=config["s3_dst"],
    shell:
        r"""
        ./vendored/upload-to-s3 \
            --quiet \
            {input.metadata:q} \
            {params.s3_dst:q}/{wildcards.lineage}/metadata.tsv.xz \
            2>&1 | tee {output.flag:q}
        """


rule upload_sequences:
    input:
        sequences="results/{lineage}/{segment}.fasta",
    output:
        flag="results/upload/{lineage}/{segment}.fasta.upload",
    params:
        s3_dst=config["s3_dst"],
    shell:
        r"""
        ./vendored/upload-to-s3 \
            --quiet \
            {input.sequences:q} \
            {params.s3_dst:q}/{wildcards.lineage}/{wildcards.segment}/sequences.fasta.xz \
            2>&1 | tee {output.flag:q}
        """
