"""
This part of the workflow handles uploading files to AWS S3.

The rule `upload_all` can be used as a target to upload all files.
"""
import os


rule upload_all:
    input:
        metadata=expand("results/upload/{lineage}/metadata.tsv.upload",
                        lineage=config['lineages']),
        sequences=expand("results/upload/{lineage}/{segment}.fasta.upload",
                         lineage=config['lineages'],
                         segment=config["segments"]),
        nextclade=[
            f"results/upload/{lineage}/{segment}/nextclade.tsv.upload"
            for lineage in (set(config["lineages"]) & set(config["nextclade"].keys()))
            for segment in (set(config["segments"]) & set(config["nextclade"][lineage].keys()))
        ],


rule upload_metadata:
    input:
        metadata="results/{lineage}/metadata.tsv",
    output:
        flag="results/upload/{lineage}/metadata.tsv.upload",
    params:
        s3_dst=config["s3_dst"],
        cloudfront_domain=config["cloudfront_domain"],
        vendored_scripts=f"{workflow.current_basedir}/../../../shared/vendored/scripts",
    benchmark:
        "benchmarks/{lineage}/upload_metadata.txt"
    shell:
        r"""
        exec &> >(tee {output.flag:q})

        {params.vendored_scripts:q}/upload-to-s3 \
            --quiet \
            {input.metadata:q} \
            {params.s3_dst:q}/{wildcards.lineage:q}/metadata.tsv.zst \
            {params.cloudfront_domain:q}
        """


rule upload_sequences:
    input:
        sequences="results/{lineage}/{segment}.fasta",
    output:
        flag="results/upload/{lineage}/{segment}.fasta.upload",
    params:
        s3_dst=config["s3_dst"],
        cloudfront_domain=config["cloudfront_domain"],
        vendored_scripts=f"{workflow.current_basedir}/../../../shared/vendored/scripts",
    benchmark:
        "benchmarks/{lineage}/{segment}/upload_sequences.txt"
    shell:
        r"""
        exec &> >(tee {output.flag:q})

        {params.vendored_scripts:q}/upload-to-s3 \
            --quiet \
            {input.sequences:q} \
            {params.s3_dst:q}/{wildcards.lineage}/{wildcards.segment}/sequences.fasta.zst \
            {params.cloudfront_domain:q}
        """


rule upload_nextclade:
    input:
        nextclade=_get_nextclade,
    output:
        flag="results/upload/{lineage}/{segment}/nextclade.tsv.upload",
    params:
        s3_dst=config["s3_dst"],
        cloudfront_domain=config["cloudfront_domain"],
        vendored_scripts=f"{workflow.current_basedir}/../../../shared/vendored/scripts",
    benchmark:
        "benchmarks/{lineage}/{segment}/upload_nextclade.txt"
    shell:
        r"""
        exec &> >(tee {output.flag:q})

        {params.vendored_scripts:q}/upload-to-s3 \
            --quiet \
            {input.nextclade:q} \
            {params.s3_dst:q}/{wildcards.lineage}/{wildcards.segment}/nextclade.tsv.zst \
            {params.cloudfront_domain:q}
        """
