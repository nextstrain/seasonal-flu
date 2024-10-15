rule upload_all_nextclade_files:
    input:
        files=lambda wildcards: [
            "data/upload/s3/{filetype}_{lineage}_{segment}.done".format(filetype=filetype, lineage=build["lineage"], segment=segment)
            for filetype in ("alignment", "nextclade")
            for build in config["builds"].values()
            for segment in build.get("segments", config["segments"])
        ]

rule get_nextclade_dataset_for_lineage_and_segment:
    output:
        nextclade_dir=directory("nextclade_dataset/{lineage}_{segment}/"),
    shell:
        """
        nextclade3 dataset get \
            -n flu_{wildcards.lineage}_{wildcards.segment} \
            --output-dir {output.nextclade_dir}
        """

rule run_nextclade:
    input:
        nextclade_dir="nextclade_dataset/{lineage}_{segment}/",
        sequences="data/{lineage}/{segment}.fasta",
    output:
        alignment="data/upload/s3/{lineage}/{segment}/aligned.fasta",
        annotations="data/upload/s3/{lineage}/{segment}/nextclade.tsv",
    log:
        "logs/run_nextclade_{lineage}_{segment}.txt"
    threads: 8
    shell:
        """
        nextclade3 run \
            -j {threads} \
            -D {input.nextclade_dir} \
            --output-fasta {output.alignment} \
            --output-tsv {output.annotations} \
            {input.sequences}
        """

rule upload_alignment:
    input:
        alignment="data/upload/s3/{lineage}/{segment}/aligned.fasta",
    output:
        flag="data/upload/s3/alignment_{lineage}_{segment}.done",
    params:
        s3_dst=config["s3_dst"],
    log:
        "logs/upload_alignment_{lineage}_{segment}.txt"
    shell:
        """
        ./scripts/upload-to-s3 \
            --quiet \
            {input.alignment:q} \
            {params.s3_dst:q}/{wildcards.lineage}/{wildcards.segment}/aligned.fasta.xz 2>&1 | tee {output.flag}
        """

rule upload_nextclade_annotations:
    input:
        annotations="data/upload/s3/{lineage}/{segment}/nextclade.tsv",
    output:
        flag="data/upload/s3/nextclade_{lineage}_{segment}.done",
    params:
        s3_dst=config["s3_dst"],
    log:
        "logs/upload_nextclade_annotations_{lineage}_{segment}.txt"
    shell:
        """
        ./scripts/upload-to-s3 \
            --quiet \
            {input.annotations:q} \
            {params.s3_dst:q}/{wildcards.lineage}/{wildcards.segment}/nextclade.tsv.xz 2>&1 | tee {output.flag}
        """
