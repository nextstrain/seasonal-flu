rule upload_all_titers:
    input:
        titers=lambda wildcards: [
            "data/upload/s3/titers/{build_name}/{titer_collection}.done".format(build_name=build_name, titer_collection=titer_collection["name"])
            for build_name, build_params in config["builds"].items()
            for titer_collection in build_params["titer_collections"]
        ]

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
    log:
        "logs/upload_titers_{build_name}_{titer_collection}.txt"
    shell:
        """
        ./ingest/vendored/upload-to-s3 \
            --quiet \
            {input.titers:q} \
            {params.s3_dst:q}/{params.lineage}/{wildcards.titer_collection}_titers.tsv.gz 2>&1 | tee {output.flag}
        """
