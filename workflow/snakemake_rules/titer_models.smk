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

rule titers_sub:
    input:
        titers = build_dir +"/{build_name}/titers/{titer_collection}.tsv",
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    params:
        translations = lambda w: [f"{build_dir}/{w.build_name}/{w.segment}/nextalign/masked.gene.{gene}_withInternalNodes.fasta" for gene in GENES[w.segment]],
        genes = lambda w: GENES[w.segment],
        attribute_prefix_argument = get_titer_collection_attribute_prefix_argument,
    output:
        titers_model = build_dir + "/{build_name}/{segment}/titers-sub-model/{titer_collection}.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/titers_sub_{build_name}_{segment}_{titer_collection}.txt",
    log:
        "logs/titers_sub_{build_name}_{segment}_{titer_collection}.txt",
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
        haplotypes="builds/{build_name}/{segment}/haplotypes.json",
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
            --haplotypes {input.haplotypes} \
            --titer-model {input.titer_model} \
            --titers {input.titers} \
            --branch-lengths {input.branch_lengths} \
            --frequencies {input.frequencies} \
            --annotations lineage={params.lineage} \
            --output {output.distances} &> {log}
        """

def get_titer_collection_title(wildcards):
    for collection in config["builds"][wildcards.build_name]["titer_collections"]:
        if collection["name"] == wildcards.titer_collection:
            return collection.get("title", collection["name"])

rule export_measurements:
    input:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains/{titer_collection}.tsv",
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
        grouping_column=["reference_strain", "clade_reference", "haplotype_reference", "source", "serum"],
        title=get_titer_collection_title,
        x_axis_label="normalized log2 titer",
        threshold=2.0,
        filters=["reference_strain", "clade_reference", "haplotype_reference", "source", "serum"],
    shell:
        """
        augur measurements export \
            --collection {input.distances} \
            --strain-column {params.strain_column} \
            --value-column {params.value_column} \
            --grouping-column {params.grouping_column} \
            --key {wildcards.titer_collection} \
            --title {params.title:q} \
            --x-axis-label {params.x_axis_label:q} \
            --threshold {params.threshold} \
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
