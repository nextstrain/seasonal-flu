"""
This part of the workflow handles downloading files from AWS S3.
"""


rule fetch_gisaid_ndjson:
    """
    Fetch previously uploaded gisaid.ndjson if it exists.
    This is a cache of the raw data from previous GISAID ingest(s).
    If it doesn't exist, then just create an empty file.
    """
    output:
        ndjson=temp("data/gisaid_cache.ndjson.zst"),
    params:
        s3_file=f"{config['s3_src']}/gisaid.ndjson.zst",
        vendored_scripts=f"{workflow.current_basedir}/../../../shared/vendored/scripts",
    shell:
        r"""
        if $({params.vendored_scripts:q}/s3-object-exists {params.s3_file:q}); then
            aws s3 cp {params.s3_file:q} {output.ndjson:q}
        else
            echo "{params.s3_file:q} does not exist, creating empty file."
            touch {output.ndjson:q}
        fi
        """


checkpoint fetch_unprocessed_files:
    """
    Fetch unprocessed GISAID files.
    These are pairs of metadata.xls.zst and sequences.fasta.zst files.

    This is a checkpoint because the DAG needs to be re-evaluated to determine
    which `gisaid_pair` need to be processed.
    """
    output:
        directory("data/unprocessed-gisaid-downloads/"),
    params:
        s3_prefix=f"{config['s3_src']}/gisaid-downloads/unprocessed/"
    shell:
        r"""
        aws s3 cp {params.s3_prefix:q} {output:q} \
            --recursive \
            --exclude "*" \
            --include "*-metadata.xls.zst" \
            --include "*-sequences.fasta.zst"
        """


rule decompress_unprocessed_files:
    input:
        metadata="data/unprocessed-gisaid-downloads/{gisaid_pair}-metadata.xls.zst",
        sequences="data/unprocessed-gisaid-downloads/{gisaid_pair}-sequences.fasta.zst",
    output:
        metadata=temp("data/{gisaid_pair}-metadata.xls"),
        sequences=temp("data/{gisaid_pair}-sequences.fasta"),
    shell:
        r"""
        zstd --decompress --stdout {input.metadata:q} > {output.metadata:q}
        zstd --decompress --stdout {input.sequences:q} > {output.sequences:q}
        """
