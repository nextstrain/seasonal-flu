rule upload_all_raw_sequences:
    input:
        sequences=lambda wildcards: [
            "data/upload/s3/raw_sequences_{lineage}_{segment}.done".format(lineage=build_params["lineage"], segment=segment)
            for build_name, build_params in config["builds"].items()
            for segment in config["segments"]
        ]

rule upload_all_sequences:
    input:
        sequences=lambda wildcards: [
            "data/upload/s3/sequences_{lineage}_{segment}.done".format(lineage=build_params["lineage"], segment=segment)
            for build_name, build_params in config["builds"].items()
            for segment in config["segments"]
        ]

rule upload_all_metadata:
    input:
        metadata=lambda wildcards: [
            "data/upload/s3/metadata_{lineage}.done".format(lineage=build_params["lineage"])
            for build_name, build_params in config["builds"].items()
        ]

rule upload_all_titers:
    input:
        titers=lambda wildcards: [
            "data/upload/s3/titers/{build_name}/{titer_collection}.done".format(build_name=build_name, titer_collection=titer_collection["name"])
            for build_name, build_params in config["builds"].items()
            for titer_collection in build_params["titer_collections"]
        ]

rule upload_raw_sequences:
    input:
        sequences="data/{lineage}/raw_{segment}.fasta",
    output:
        flag="data/upload/s3/raw_sequences_{lineage}_{segment}.done",
    params:
        s3_dst=config["s3_dst"],
        vendored_scripts=f"{workflow.current_basedir}/../../shared/vendored/scripts",
    log:
        "logs/upload_raw_sequences_{lineage}_{segment}.txt"
    shell:
        """
        {params.vendored_scripts:q}/upload-to-s3 \
            --quiet \
            {input.sequences:q} \
            {params.s3_dst:q}/{wildcards.lineage}/{wildcards.segment}/raw_sequences.fasta.xz 2>&1 | tee {output.flag}
        """

rule upload_sequences:
    input:
        sequences="data/{lineage}/{segment}.fasta",
    output:
        flag="data/upload/s3/sequences_{lineage}_{segment}.done",
    params:
        s3_dst=config["s3_dst"],
        vendored_scripts=f"{workflow.current_basedir}/../../shared/vendored/scripts",
    log:
        "logs/upload_sequences_{lineage}_{segment}.txt"
    shell:
        """
        {params.vendored_scripts:q}/upload-to-s3 \
            --quiet \
            {input.sequences:q} \
            {params.s3_dst:q}/{wildcards.lineage}/{wildcards.segment}/sequences.fasta.xz 2>&1 | tee {output.flag}
        """

rule upload_metadata:
    input:
        metadata="data/{lineage}/metadata.tsv",
    output:
        flag="data/upload/s3/metadata_{lineage}.done",
    params:
        s3_dst=config["s3_dst"],
        vendored_scripts=f"{workflow.current_basedir}/../../shared/vendored/scripts",
    log:
        "logs/upload_metadata_{lineage}.txt"
    shell:
        """
        {params.vendored_scripts:q}/upload-to-s3 \
            --quiet \
            {input.metadata:q} \
            {params.s3_dst:q}/{wildcards.lineage}/metadata.tsv.xz 2>&1 | tee {output.flag}
        """

def get_titer_collection_data(wildcards):
    for collection in config["builds"][wildcards.build_name]["titer_collections"]:
        if collection["name"] == wildcards.titer_collection:
            return collection["data"]

rule upload_titers:
    input:
        titers=get_titer_collection_data,
    output:
        flag="data/upload/s3/titers/{build_name}/{titer_collection}.done",
    params:
        s3_dst=config["s3_dst"],
        lineage=lambda wildcards: config["builds"][wildcards.build_name]["lineage"],
        vendored_scripts=f"{workflow.current_basedir}/../../shared/vendored/scripts",
    log:
        "logs/upload_titers_{build_name}_{titer_collection}.txt"
    shell:
        """
        {params.vendored_scripts:q}/upload-to-s3 \
            --quiet \
            {input.titers:q} \
            {params.s3_dst:q}/{params.lineage}/{wildcards.titer_collection}_titers.tsv.gz 2>&1 | tee {output.flag}
        """
