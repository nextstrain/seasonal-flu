ruleorder: export_private > export

def get_antigenic_plot_paths(wildcards):
    paths = []
    for build_name in config["builds"].keys():
        for collection in config["builds"][build_name]["titer_collections"]:
            if "ferret" in collection["data"]:
                paths.append(f"builds/{build_name}/ha/plots/antigenic_distances_between_strains_{build_name}_{collection['name']}.png")

    return paths

rule all_antigenic_plots:
    input:
        get_antigenic_plot_paths,

rule plot_antigenic_distances_between_strains:
    input:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains/{titer_collection}.tsv",
        clades=lambda wildcards: f"config/subclades_for_titer_plots_{config['builds'][wildcards.build_name]['lineage']}.txt",
        references=lambda wildcards: f"config/references_for_titer_plots_{config['builds'][wildcards.build_name]['lineage']}.txt",
        colors=lambda wildcards: f"config/colors_for_titer_plots_{config['builds'][wildcards.build_name]['lineage']}.tsv",
    output:
        plot="builds/{build_name}/{segment}/plots/antigenic_distances_between_strains_{build_name}_{titer_collection}.png",
    benchmark:
        "benchmarks/plot_antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}.txt"
    log:
        "logs/plot_antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}.txt"
    params:
        min_test_date=2022.5,
        title=get_titer_collection_title,
        clade_color_field="subclade_test",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/plot_antigenic_distances_between_strains.py \
            --antigenic-distances {input.distances} \
            --min-test-date {params.min_test_date} \
            --clades {input.clades} \
            --clade-color-field {params.clade_color_field} \
            --references {input.references} \
            --colors {input.colors} \
            --title {params.title:q} \
            --output {output.plot} 2>&1 | tee {log}
        """

rule annotate_titer_counts_for_reference_viruses:
    input:
        tree="builds/{build_name}/{segment}/tree.nwk",
        titers="builds/{build_name}/titers/{titer_collection}.tsv",
    output:
        titer_counts="builds/{build_name}/{segment}/titers_for_reference_viruses/{titer_collection}.json",
    benchmark:
        "benchmarks/annotate_titer_counts_for_reference_viruses_{build_name}_{segment}_{titer_collection}.txt"
    log:
        "logs/annotate_titer_counts_for_reference_viruses_{build_name}_{segment}_{titer_collection}.txt"
    params:
        attribute_name=lambda wildcards: f"titers_for_reference_viruses_{wildcards.titer_collection}",
        attribute_name_for_is_reference_virus=lambda wildcards: f"is_titer_reference_virus_in_{wildcards.titer_collection}",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/annotate_titers_per_node.py \
            --tree {input.tree} \
            --titers {input.titers} \
            --attribute-name {params.attribute_name} \
            --attribute-name-for-is-reference-virus {params.attribute_name_for_is_reference_virus} \
            --use-references \
            --output {output.titer_counts} 2>&1 | tee {log}
        """

rule summarize_haplotype_titer_coverage:
    input:
        haplotypes="builds/{build_name}/{segment}/haplotypes.json",
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains/{titer_collection}.tsv",
        frequencies="builds/{build_name}/{segment}/tip-frequencies.json",
    output:
        table="builds/{build_name}/{segment}/haplotype_summary/{titer_collection}.tsv",
        node_data="builds/{build_name}/{segment}/haplotypes_without_references/{titer_collection}.json",
    benchmark:
        "benchmarks/summarize_haplotype_titer_coverage_{build_name}_{segment}_{titer_collection}.txt"
    log:
        "logs/summarize_haplotype_titer_coverage_{build_name}_{segment}_{titer_collection}.txt"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/summarize_haplotype_titer_coverage.py \
            --haplotypes {input.haplotypes} \
            --antigenic-distances {input.distances} \
            --frequencies {input.frequencies} \
            --attribute-name-for-haplotype-without-reference "haplotype_missing_reference_virus_in_{wildcards.titer_collection}" \
            --output-table {output.table} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

def get_private_node_data(wildcards):
    node_data = []

    # Only try to annotate titer collections for HA.
    if wildcards.segment == "ha":
        for collection in config["builds"][wildcards.build_name]["titer_collections"]:
            node_data.append(f"builds/{wildcards.build_name}/{wildcards.segment}/titers_for_reference_viruses/{collection['name']}.json")
            node_data.append(f"builds/{wildcards.build_name}/{wildcards.segment}/haplotypes_without_references/{collection['name']}.json")

    return node_data

rule export_private:
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
        node_data = _get_node_data_by_wildcards,
        private_node_data = get_private_node_data,
        auspice_config = lambda w: config['builds'][w.build_name]['auspice_config'],
        lat_longs = config['lat-longs']
    output:
        auspice_json = "auspice/{build_name}_{segment}.json",
        root_sequence_json = "auspice/{build_name}_{segment}_root-sequence.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/export_{build_name}_{segment}.txt"
    log:
        "logs/export_{build_name}_{segment}.txt"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} {input.private_node_data} \
            --include-root-sequence \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --minify-json \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
