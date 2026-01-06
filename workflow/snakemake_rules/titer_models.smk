'''
This file contains rules infer titer models

input:
 - builds/{build_name}/metadata.tsv
 - builds/{build_name}/titers/{titer_collection}.tsv
 - builds/{build_name}/{segment}/tree.nwk

output:

'''
build_dir = config.get("build_dir", "builds")

def get_titer_collection_attribute_prefix(wildcards):
    for collection in config["builds"][wildcards.build_name]["titer_collections"]:
        if collection["name"] == wildcards.titer_collection:
            return collection.get("prefix", "")

def get_titer_collection_attribute_prefix_argument(wildcards):
    collection_prefix = get_titer_collection_attribute_prefix(wildcards)
    if collection_prefix:
        return f"--attribute-prefix {collection_prefix}"
    else:
        return ""

def get_titer_collection_genes(wildcards):
    for collection in config["builds"][wildcards.build_name]["titer_collections"]:
        if collection["name"] == wildcards.titer_collection:
            return collection.get("genes", GENES[wildcards.segment])

rule titers_sub:
    input:
        titers = build_dir +"/{build_name}/titers/{titer_collection}.tsv",
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    params:
        genes = get_titer_collection_genes,
        translations = lambda wildcards: [f"{build_dir}/{wildcards.build_name}/{wildcards.segment}/translations/{gene}_withInternalNodes.fasta" for gene in get_titer_collection_genes(wildcards)],
        attribute_prefix_argument = get_titer_collection_attribute_prefix_argument,
    output:
        titers_model = build_dir + "/{build_name}/{segment}/titers-sub-model/{titer_collection}.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/titers_sub_{build_name}_{segment}_{titer_collection}.txt",
    log:
        "logs/titers_sub_{build_name}_{segment}_{titer_collection}.txt",
    resources:
        mem_mb=8000,
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {params.translations} \
            --gene-names {params.genes} \
            --tree {input.tree} \
            --allow-empty-model \
            {params.attribute_prefix_argument} \
            --output {output.titers_model} 2>&1 | tee {log}
        """

rule titers_tree:
    input:
        titers = "builds/{build_name}/titers/{titer_collection}.tsv",
        tree = rules.refine.output.tree,
    params:
        attribute_prefix_argument = get_titer_collection_attribute_prefix_argument,
    output:
        titers_model = "builds/{build_name}/{segment}/titers-tree-model/{titer_collection}.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/titers_tree_{build_name}_{segment}_{titer_collection}.txt",
    log:
        "logs/titers_tree_{build_name}_{segment}_{titer_collection}.txt",
    resources:
        mem_mb=8000,
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --allow-empty-model \
            {params.attribute_prefix_argument} \
            --output {output.titers_model} 2>&1 | tee {log}
        """

rule antigenic_distances_between_strains:
    input:
        tree="builds/{build_name}/{segment}/tree.nwk",
        clades="builds/{build_name}/{segment}/clades.json",
        subclades="builds/{build_name}/{segment}/subclades.json",
        emerging_haplotypes="builds/{build_name}/{segment}/emerging_haplotypes.json",
        derived_haplotypes="builds/{build_name}/{segment}/derived_haplotypes.json",
        titer_model="builds/{build_name}/{segment}/titers-sub-model/{titer_collection}.json",
        titers="builds/{build_name}/titers/{titer_collection}.tsv",
        branch_lengths="builds/{build_name}/{segment}/branch-lengths.json",
        frequencies="builds/{build_name}/{segment}/tip-frequencies.json",
    output:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains/{titer_collection}.tsv",
    benchmark:
        "benchmarks/antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}.txt"
    log:
        "logs/antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}.txt"
    conda: "../envs/nextstrain.yaml"
    params:
        lineage=lambda wildcards: config["builds"][wildcards.build_name].get("lineage", ""),
    shell:
        """
        python3 scripts/get_antigenic_distances_between_strains.py \
            --tree {input.tree} \
            --clades {input.clades} \
            --subclades {input.subclades} \
            --emerging-haplotypes {input.emerging_haplotypes} \
            --derived-haplotypes {input.derived_haplotypes} \
            --titer-model {input.titer_model} \
            --titers {input.titers} \
            --branch-lengths {input.branch_lengths} \
            --frequencies {input.frequencies} \
            --annotations lineage={params.lineage} \
            --output {output.distances} &> {log}
        """

rule generate_collection_config_json:
    input:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains/{titer_collection}.tsv",
        tree="builds/{build_name}/{segment}/tree.nwk",
    output:
        config_json="builds/{build_name}/{segment}/measurements_collection_config/{titer_collection}.json",
    conda: "../envs/nextstrain.yaml"
    params:
        groupings=[
            "reference_strain",
            "reference_strain_source",
            "subclade_reference",
            "emerging_haplotype_reference",
            "derived_haplotype_reference",
            "source",
            "serum"
        ],
        fields=[
            "strain",
            "reference_strain",
            "reference_strain_source",
            "serum",
            "value",
            "raw_titer",
            "source",
            "test_date",
            "reference_date",
            "subclade_test",
            "subclade_reference",
            "emerging_haplotype_test",
            "emerging_haplotype_reference",
            "derived_haplotype_test",
            "derived_haplotype_reference",
        ],
    log:
        "logs/generate_collection_config_json_{build_name}_{segment}_{titer_collection}.txt"
    shell:
        """
        python3 scripts/generate_collection_config_json.py \
            --tree {input.tree} \
            --collection {input.distances} \
            --groupings {params.groupings:q} \
            --fields {params.fields:q} \
            --output {output.config_json} &> {log}
        """

def get_titer_collection_title(wildcards):
    for collection in config["builds"][wildcards.build_name]["titer_collections"]:
        if collection["name"] == wildcards.titer_collection:
            return collection.get("title", collection["name"])

rule export_measurements:
    input:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains/{titer_collection}.tsv",
        configuration="builds/{build_name}/{segment}/measurements_collection_config/{titer_collection}.json",
    output:
        measurements="builds/{build_name}/{segment}/measurements/{titer_collection}.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/export_measurements_{build_name}_{segment}_{titer_collection}.txt"
    log:
        "logs/export_measurements_{build_name}_{segment}_{titer_collection}.txt"
    params:
        strain_column="test_strain",
        value_column="log2_titer",
        title=get_titer_collection_title,
        x_axis_label="normalized log2 titer",
        thresholds=[0.0, 2.0],
        filters=[
            "reference_strain",
            "reference_strain_source",
            "subclade_reference",
            "emerging_haplotype_reference",
            "derived_haplotype_reference",
            "source",
            "serum"
        ],
        include_columns=[
            "test_strain",
            "reference_strain",
            "reference_strain_source",
            "serum",
            "log2_titer",
            "raw_titer",
            "source",
            "test_date",
            "reference_date",
            "subclade_test",
            "subclade_reference",
            "emerging_haplotype_test",
            "emerging_haplotype_reference",
            "derived_haplotype_test",
            "derived_haplotype_reference",
        ],
    shell:
        """
        augur measurements export \
            --collection {input.distances} \
            --collection-config {input.configuration} \
            --include-columns {params.include_columns:q} \
            --strain-column {params.strain_column} \
            --value-column {params.value_column} \
            --key {wildcards.titer_collection} \
            --title {params.title:q} \
            --x-axis-label {params.x_axis_label:q} \
            --thresholds {params.thresholds} \
            --filters {params.filters} \
            --show-threshold \
            --hide-overall-mean \
            --minify-json \
            --output-json {output.measurements} 2>&1 | tee {log}
        """

checkpoint get_titers_per_reference:
    input:
        titers="builds/{build_name}/titers/{titer_collection}.tsv",
    output:
        references="builds/{build_name}/titer_references/{titer_collection}.txt",
        reference_titers_directory=directory("builds/{build_name}/reference_titers/{titer_collection}/"),
    conda: "../envs/nextstrain.yaml"
    shell:
        r"""
        mkdir -p {output.reference_titers_directory};

        python scripts/get_titers_per_reference.py \
            --titers {input.titers} \
            --output-references {output.references} \
            --output-titers-directory {output.reference_titers_directory}
        """

rule reference_model_titers_sub:
    input:
        titers = build_dir +"/{build_name}/reference_titers/{titer_collection}/{reference}.tsv",
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    params:
        genes = get_titer_collection_genes,
        translations = lambda wildcards: [f"{build_dir}/{wildcards.build_name}/{wildcards.segment}/translations/{gene}_withInternalNodes.fasta" for gene in get_titer_collection_genes(wildcards)],
        attribute_prefix_argument = get_titer_collection_attribute_prefix_argument,
    output:
        titers_model = build_dir + "/{build_name}/{segment}/reference-titers-sub-model/{titer_collection}/{reference}.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/titers_sub_{build_name}_{segment}_{titer_collection}_{reference}.txt",
    log:
        "logs/titers_sub_{build_name}_{segment}_{titer_collection}_{reference}.txt",
    resources:
        mem_mb=8000,
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {params.translations} \
            --gene-names {params.genes} \
            --tree {input.tree} \
            --allow-empty-model \
            {params.attribute_prefix_argument} \
            --output {output.titers_model} 2>&1 | tee {log}
        """

rule reference_model_antigenic_distances_between_strains:
    input:
        titer_model="builds/{build_name}/{segment}/reference-titers-sub-model/{titer_collection}/{reference}.json",
        titers="builds/{build_name}/reference_titers/{titer_collection}/{reference}.tsv",
    output:
        distances="builds/{build_name}/{segment}/reference_model_antigenic_distances_between_strains/{titer_collection}/{reference}.tsv",
    benchmark:
        "benchmarks/reference_model_antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}_{reference}.txt"
    log:
        "logs/reference_model_antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}_{reference}.txt"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/get_antigenic_distances_for_reference_model.py \
            --titer-model {input.titer_model} \
            --titers {input.titers} \
            --output {output.distances} &> {log}
        """

def aggregate_reference_model_distances_input(wildcards):
    with checkpoints.get_titers_per_reference.get(**wildcards).output["references"].open() as fh:
        distances = [
            f"builds/{wildcards.build_name}/{wildcards.segment}/reference_model_antigenic_distances_between_strains/{wildcards.titer_collection}/{reference.strip()}.tsv"
            for reference in fh
        ]

    return distances

rule aggregate_reference_model_distances:
    input:
        distances=aggregate_reference_model_distances_input,
    output:
        distances="builds/{build_name}/{segment}/reference_model_antigenic_distances_between_strains/{titer_collection}.tsv",
    conda: "../envs/nextstrain.yaml"
    shell:
        r"""
        tsv-append -H {input.distances} > {output.distances}
        """

rule generate_reference_model_collection_config_json:
    input:
        distances="builds/{build_name}/{segment}/reference_model_antigenic_distances_between_strains/{titer_collection}.tsv",
        tree="builds/{build_name}/{segment}/tree.nwk",
    output:
        config_json="builds/{build_name}/{segment}/reference_model_measurements_collection_config/{titer_collection}.json",
    conda: "../envs/nextstrain.yaml"
    params:
        groupings=[
            "reference_strain",
        ],
        fields=[
            "strain",
            "reference_strain",
            "value",
        ],
    log:
        "logs/generate_reference_model_collection_config_json_{build_name}_{segment}_{titer_collection}.txt"
    shell:
        """
        python3 scripts/generate_collection_config_json.py \
            --tree {input.tree} \
            --collection {input.distances} \
            --groupings {params.groupings:q} \
            --fields {params.fields:q} \
            --output {output.config_json} &> {log}
        """

rule export_reference_model_measurements:
    input:
        distances="builds/{build_name}/{segment}/reference_model_antigenic_distances_between_strains/{titer_collection}.tsv",
        configuration="builds/{build_name}/{segment}/reference_model_measurements_collection_config/{titer_collection}.json",
    output:
        measurements="builds/{build_name}/{segment}/reference_model_measurements/{titer_collection}.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/export_reference_model_measurements_{build_name}_{segment}_{titer_collection}.txt"
    log:
        "logs/export_reference_model_measurements_{build_name}_{segment}_{titer_collection}.txt"
    params:
        strain_column="test_strain",
        value_column="log2_titer",
        title=lambda wildcards: get_titer_collection_title(wildcards) + " (inferred)",
        x_axis_label="inferred log2 titer",
        thresholds=[0.0, 2.0],
        filters=[
            "reference_strain",
        ],
        include_columns=[
            "reference_strain",
        ],
    shell:
        """
        augur measurements export \
            --collection {input.distances} \
            --collection-config {input.configuration} \
            --include-columns {params.include_columns:q} \
            --strain-column {params.strain_column} \
            --value-column {params.value_column} \
            --key {wildcards.titer_collection}_inferred \
            --title {params.title:q} \
            --x-axis-label {params.x_axis_label:q} \
            --thresholds {params.thresholds} \
            --filters {params.filters} \
            --show-threshold \
            --hide-overall-mean \
            --minify-json \
            --output-json {output.measurements} 2>&1 | tee {log}
        """

def get_titer_collections(wildcards):
    files = []
    for collection in config["builds"][wildcards.build_name]["titer_collections"]:
        files.append(f"builds/{wildcards.build_name}/{wildcards.segment}/measurements/{collection['name']}.json")

        if collection.get("run_reference_models"):
            files.append(f"builds/{wildcards.build_name}/{wildcards.segment}/reference_model_measurements/{collection['name']}.json")

    return files

rule concat_measurements:
    input:
        measurements=get_titer_collections,
    output:
        measurements="auspice/{build_name}_{segment}_measurements.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/concat_measurements_{build_name}_{segment}.txt"
    log:
        "logs/concat_measurements_{build_name}_{segment}.txt"
    shell:
        """
        augur measurements concat \
            --jsons {input.measurements} \
            --minify-json \
            --output-json {output.measurements} 2>&1 | tee {log}
        """
