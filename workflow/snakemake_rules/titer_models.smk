'''
This file contains rules infer titer models

input:
 - builds/{build_name}/metadata.tsv
 - builds/{build_name}/titers.tsv
 - builds/{build_name}/{segment}/tree.nwk

output:

'''
build_dir = config.get("build_dir", "builds")

rule titers_sub:
    input:
        titers = build_dir +"/{build_name}/titers.tsv",
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    params:
        translations = lambda w: [f"{build_dir}/{w.build_name}/{w.segment}/nextalign/masked.gene.{gene}_withInternalNodes.fasta" for gene in GENES[w.segment]],
        genes = lambda w: GENES[w.segment]
    output:
        titers_model = build_dir + "/{build_name}/{segment}/titers-sub-model.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/titers_sub_{build_name}_{segment}.txt",
    log:
        "logs/titers_sub_{build_name}_{segment}.txt",
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {params.translations} \
            --gene-names {params.genes} \
            --tree {input.tree} \
            --allow-empty-model \
            --output {output.titers_model} 2>&1 | tee {log}
        """

rule titers_tree:
    input:
        titers = build_dir +"/{build_name}/titers.tsv",
        tree = rules.refine.output.tree
    output:
        titers_model = build_dir + "/{build_name}/{segment}/titers-tree-model.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/titers_tree_{build_name}_{segment}.txt",
    log:
        "logs/titers_tree_{build_name}_{segment}.txt",
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --allow-empty-model \
            --output {output.titers_model} 2>&1 | tee {log}
        """

rule antigenic_distances_between_strains:
    input:
        tree="builds/{build_name}/{segment}/tree.nwk",
        clades="builds/{build_name}/{segment}/clades.json",
        haplotypes="builds/{build_name}/{segment}/haplotypes.json",
        titer_model="builds/{build_name}/{segment}/titers-sub-model.json",
        titers="builds/{build_name}/titers.tsv",
        branch_lengths="builds/{build_name}/{segment}/branch-lengths.json",
        frequencies="builds/{build_name}/{segment}/tip-frequencies.json",
    output:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains.tsv",
    benchmark:
        "benchmarks/antigenic_distances_between_strains_{build_name}_{segment}.txt"
    log:
        "logs/antigenic_distances_between_strains_{build_name}_{segment}.txt"
    conda: "../envs/nextstrain.yaml"
    params:
        lineage=lambda wildcards: config["builds"][wildcards.build_name].get("lineage", ""),
        passage=lambda wildcards: config["builds"][wildcards.build_name].get("passage", ""),
        assay=lambda wildcards: config["builds"][wildcards.build_name].get("assay", ""),
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
            --annotations lineage={params.lineage} passage={params.passage} assay={params.assay} \
            --output {output.distances} &> {log}
        """

rule export_measurements:
    input:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains.tsv",
    output:
        measurements="auspice/{build_name}_{segment}_measurements.json",
    conda: "../envs/nextstrain.yaml"
    params:
        strain_column="test_strain",
        value_column="log2_titer",
        grouping_column=["reference_strain", "clade_reference", "haplotype_reference", "source", "serum"],
        key=lambda wildcards: f"{config['builds'][wildcards.build_name]['lineage']}_{wildcards.segment}_{config['builds'][wildcards.build_name]['passage']}_{config['builds'][wildcards.build_name]['assay']}",
        title=lambda wildcards: f"{lineage_name_by_abbreviation[config['builds'][wildcards.build_name]['lineage']]} {config['builds'][wildcards.build_name]['passage']}-passaged {config['builds'][wildcards.build_name]['assay'].upper()} measurements",
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
            --key {params.key} \
            --title {params.title:q} \
            --x-axis-label {params.x_axis_label:q} \
            --threshold {params.threshold} \
            --filters {params.filters} \
            --show-threshold \
            --hide-overall-mean \
            --minify-json \
            --output-json {output.measurements}
        """
