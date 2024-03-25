nextclade_dataset_by_lineage_and_segment = {
    "h1n1pdm": {
        "ha": "nextstrain/flu/h1n1pdm/ha/california-7-2009",
    },
    "h3n2": {
        "ha": "nextstrain/flu/h3n2/ha/wisconsin-67-2005",
    },
    "vic": {
        "ha": "nextstrain/flu/vic/ha/brisbane-60-2008",
    },
}

rule upload_all_nextclade_files:
    input:
        files=lambda wildcards: [
            "data/upload/s3/{filetype}_{lineage}_{segment}.done".format(filetype=filetype, lineage=lineage, segment=segment)
            for filetype in ("alignment", "nextclade")
            for lineage in nextclade_dataset_by_lineage_and_segment.keys()
            for segment in nextclade_dataset_by_lineage_and_segment[lineage].keys()
        ]

rule get_nextclade_dataset_for_lineage_and_segment:
    output:
        nextclade_dir=directory("nextclade_dataset/{lineage}_{segment}/"),
    params:
        dataset_name=lambda wildcards: nextclade_dataset_by_lineage_and_segment.get(wildcards.lineage, {}).get(wildcards.segment),
    shell:
        """
        nextclade3 dataset get \
            -n {params.dataset_name:q} \
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
