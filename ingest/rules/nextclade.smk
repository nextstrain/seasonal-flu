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
    threads: 1
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade3 run \
            {input.sequences:q} \
            --jobs {threads} \
            --input-dataset {input.dataset:q} \
            --output-tsv {output.nextclade:q}
        """


rule add_derived_haplotypes:
    input:
        nextclade="results/{dataset}/{segment}/nextclade.tsv",
    output:
        haplotypes="results/{dataset}/{segment}/nextclade_with_derived_haplotypes.tsv",
    params:
        genes=lambda w: _get_nextclade_config(w)["derived_haplotypes"]["genes"],
        clade_column=lambda w: _get_nextclade_config(w)["derived_haplotypes"]["clade_column"],
        mutations_column=lambda w: _get_nextclade_config(w)["derived_haplotypes"]["mutations_column"],
        derived_haplotype_column=lambda w: _get_nextclade_config(w)["derived_haplotypes"]["derived_haplotype_column"],
    benchmark:
        "benchmarks/{dataset}/{segment}/add_derived_haplotypes.txt"
    log:
        "logs/{dataset}/{segment}/add_derived_haplotypes.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 ../scripts/add_derived_haplotypes.py \
            --nextclade {input.nextclade:q} \
            --genes {params.genes:q} \
            --strip-genes \
            --clade-column {params.clade_column:q} \
            --mutations-column {params.mutations_column:q} \
            --attribute-name {params.derived_haplotype_column:q} \
            --output {output.haplotypes:q}
        """


rule add_emerging_haplotypes:
    input:
        # Adds to derived haplotypes file if it was also requested
        nextclade=lambda w: (
            "results/{dataset}/{segment}/nextclade_with_derived_haplotypes.tsv"
            if _get_nextclade_config(w).get("derived_haplotypes")
            else "results/{dataset}/{segment}/nextclade.tsv"),
        haplotypes=lambda w: _get_nextclade_config(w)["emerging_haplotypes"]["haplotypes"],
    output:
        haplotypes="results/{dataset}/{segment}/nextclade_with_emerging_haplotypes.tsv",
    params:
        clade_column=lambda w: _get_nextclade_config(w)["emerging_haplotypes"]["clade_column"],
        emerging_haplotype_column=lambda w: _get_nextclade_config(w)["emerging_haplotypes"]["emerging_haplotype_column"],
    benchmark:
        "benchmarks/{dataset}/{segment}/add_emerging_haplotypes.txt"
    log:
        "logs/{dataset}/{segment}/add_emerging_haplotypes.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python ../scripts/assign_haplotypes.py \
            --substitutions {input.nextclade:q} \
            --haplotypes {input.haplotypes:q} \
            --clade-column {params.clade_column:q} \
            --haplotype-column-name {params.emerging_haplotype_column:q} \
            --output-table {output.haplotypes:q}
        """


def _get_nextclade(wildcards):
    """
    The Nextclade file used depends on whether haplotypes were requested.
    Emerging haplotypes get added to the derived haplotypes file if both were
    requested, so just use the emerging haplotypes file.
    """
    if _get_nextclade_config(wildcards).get("emerging_haplotypes"):
        return "results/{dataset}/{segment}/nextclade_with_emerging_haplotypes.tsv"
    elif _get_nextclade_config(wildcards).get("derived_haplotypes"):
        return "results/{dataset}/{segment}/nextclade_with_derived_haplotypes.tsv"
    else:
        return "results/{dataset}/{segment}/nextclade.tsv"


def _get_nextclade_field_map(wildcards):
    """
    Automatically adds the emerging_haplotype column and derived haplotype
    column to the field map if they were provided in the config and were not
    already defined in the field map.
    """
    nextclade_config = _get_nextclade_config(wildcards)
    field_map = nextclade_config["field_map"].copy()

    if emerging_haplotype_field := nextclade_config.get("emerging_haplotypes", {}).get("emerging_haplotype_column"):
        if emerging_haplotype_field not in field_map:
            field_map[emerging_haplotype_field] = emerging_haplotype_field

    if derived_haplotype_field := nextclade_config.get("derived_haplotypes", {}).get("derived_haplotype_column"):
        if derived_haplotype_field not in field_map:
            field_map[derived_haplotype_field] = derived_haplotype_field

    return field_map


rule nextclade_metadata:
    input:
        nextclade=_get_nextclade,
    output:
        nextclade_metadata=temp("results/{dataset}/{segment}/nextclade_metadata.tsv"),
    params:
        nextclade_id_field=lambda w: _get_nextclade_config(w)["id_field"],
        nextclade_field_map=lambda w: [f"{old}={new}" for old, new in _get_nextclade_field_map(w).items()],
        nextclade_fields=lambda w: ",".join(_get_nextclade_field_map(w).values()),
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
