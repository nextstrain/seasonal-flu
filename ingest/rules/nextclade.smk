"""
This part of the workflow handles running Nextclade on the curated metadata
and sequences.
"""


def _get_nextclade_config(wildcards):
    return config["nextclade"][wildcards.dataset][wildcards.segment]


rule get_nextclade_dataset:
    """Download Nextclade dataset"""
    output:
        dataset="data/nextclade/{dataset}/{segment}.zip",
    params:
        nextclade_dataset=lambda w: _get_nextclade_config(w)["dataset_name"],
    benchmark:
        "benchmarks/get_nextclade_dataset/{dataset}/{segment}.txt"
    log:
        "logs/get_nextclade_dataset/{dataset}/{segment}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade3 dataset get \
            --name={params.nextclade_dataset:q} \
            --output-zip={output.dataset:q} \
            --verbose
        """


rule run_nextclade:
    input:
        dataset="data/nextclade/{dataset}/{segment}.zip",
        sequences="results/{dataset}/{segment}.fasta",
    output:
        nextclade="results/{dataset}/{segment}/nextclade.tsv",
    benchmark:
        "benchmarks/run_nextclade/{dataset}/{segment}.txt"
    log:
        "logs/run_nextclade/{dataset}/{segment}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade3 run \
            {input.sequences:q} \
            --input-dataset {input.dataset:q} \
            --output-tsv {output.nextclade:q}
        """


rule nextclade_metadata:
    input:
        nextclade="results/{dataset}/{segment}/nextclade.tsv",
    output:
        nextclade_metadata=temp("results/{dataset}/{segment}/nextclade_metadata.tsv"),
    params:
        nextclade_id_field=lambda w: _get_nextclade_config(w)["id_field"],
        nextclade_field_map=lambda w: [f"{old}={new}" for old, new in _get_nextclade_config(w)["field_map"].items()],
        nextclade_fields=lambda w: ",".join(_get_nextclade_config(w)["field_map"].values()),
    benchmark:
        "benchmarks/{dataset}/{segment}/nextclade_metadata.txt"
    log:
        "logs/{dataset}/{segment}/nextclade_metadata.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur curate rename \
            --metadata {input.nextclade:q} \
            --id-column {params.nextclade_id_field:q} \
            --field-map {params.nextclade_field_map:q} \
            --output-metadata - \
        | csvtk cut -t --fields {params.nextclade_fields:q} \
        > {output.nextclade_metadata:q}
        """


def _get_nextclade_metadata(wildcards):
    return {
        segment: f"results/{wildcards.dataset}/{segment}/nextclade_metadata.tsv"
        for segment in config["segments"]
        if segment in config["nextclade"][wildcards.dataset]
    }


rule join_metadata_and_nextclade:
    input:
        unpack(_get_nextclade_metadata),
        metadata="data/{dataset}/subset_metadata.tsv",
    output:
        metadata="data/{dataset}/metadata_with_nextclade.tsv",
    params:
        metadata_id_field=config["curate"]["output_id_field"],
        # augur merge requires named inputs
        named_nextclade_metadata=lambda w, input: (
            f"{segment}={path}"
            for segment, path in input.items()
            if segment != "metadata"
        ),
        nextclade_id_fields=lambda w, input: (
            f"{segment}={config['nextclade'][w.dataset][segment]['id_field']}"
            for segment in input.keys()
            if segment != "metadata"
        )
    benchmark:
        "benchmarks/{dataset}/join_metadata_and_nextclade.txt"
    log:
        "logs/{dataset}/join_metadata_and_nextclade.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur merge \
            --metadata \
                metadata={input.metadata:q} \
                {params.named_nextclade_metadata} \
            --metadata-id-columns \
                metadata={params.metadata_id_field:q} \
                {params.nextclade_id_fields} \
            --output-metadata {output.metadata:q} \
            --no-source-columns
        """


rule final_metadata:
    input:
        metadata=lambda w: (
            "data/{dataset}/metadata_with_nextclade.tsv"
            if config["nextclade"].get(w.dataset)
            else "data/{dataset}/subset_metadata.tsv"),
    output:
        metadata="results/{dataset}/metadata.tsv"
    shell:
        r"""
        cp {input.metadata:q} {output.metadata:q}
        """
