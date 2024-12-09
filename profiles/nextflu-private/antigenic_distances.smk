ruleorder: export_private > export

def get_antigenic_plot_paths(wildcards):
    paths = []
    for build_name in config["builds"].keys():
        if "titers" in build_name:
            for collection in config["builds"][build_name]["titer_collections"]:
                if "ferret" in collection["data"]:
                    paths.append(f"figures/antigenic_distances_between_strains_{build_name}_ha_{collection['name']}.png")

    return paths

rule all_antigenic_plots:
    input:
        get_antigenic_plot_paths,

rule plot_antigenic_distances_between_strains:
    input:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains/{titer_collection}.tsv",
        clades=lambda wildcards: f"config/subclades_for_titer_plots_{config['builds'][wildcards.build_name]['lineage']}.txt",
        references=lambda wildcards: f"config/references_for_titer_plots/{config['builds'][wildcards.build_name]['lineage']}/{wildcards.titer_collection}.txt",
        auspice_config=lambda wildcards: f"profiles/nextflu-private/{config['builds'][wildcards.build_name]['lineage']}/{wildcards.segment}/auspice_config.json",
    output:
        plot="figures/antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}.png",
    benchmark:
        "benchmarks/plot_antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}.txt"
    log:
        "logs/plot_antigenic_distances_between_strains_{build_name}_{segment}_{titer_collection}.txt"
    params:
        min_test_date=2024.4167,
        title=get_titer_collection_title,
        clade_color_field="proposed_subclade_test",
        auspice_config_color_field="proposed_subclade",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/plot_antigenic_distances_between_strains.py \
            --antigenic-distances {input.distances} \
            --min-test-date {params.min_test_date} \
            --clades {input.clades} \
            --clade-color-field {params.clade_color_field} \
            --references {input.references} \
            --auspice-config {input.auspice_config} \
            --auspice-config-color-field {params.auspice_config_color_field} \
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
        markdown_table="builds/{build_name}/{segment}/haplotype_summary/{titer_collection}.md",
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
            --output-markdown-table {output.markdown_table} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

rule scores:
    input:
        metadata = "builds/{build_name}/metadata.tsv",
        tree = "builds/{build_name}/{segment}/tree.nwk",
    output:
        node_data = "builds/{build_name}/{segment}/scores.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/scores.py \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --output {output}
        """

rule welsh_epitope_distances:
    input:
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done",
        distance_maps = [
            "config/distance_maps/h3n2/ha/welsh_epitope_sites.json",
            "config/distance_maps/h3n2/ha/welsh_escape_by_site_and_amino_acid.json",
        ],
    output:
        distances = "builds/{build_name}/{segment}/welsh_epitope_distances.json",
    params:
        alignments = lambda w: [f"{build_dir}/{w.build_name}/{w.segment}/translations/{gene}_withInternalNodes.fasta" for gene in GENES[w.segment]],
        genes = lambda w: GENES[w.segment],
        comparisons = ["root", "root"],
        attribute_names = ["welsh_ep", "welsh_escape"],
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/welsh_epitope_distances_{build_name}_{segment}.txt"
    log:
        "logs/welsh_epitope_distances_{build_name}_{segment}.txt"
    resources:
        mem_mb=8000,
        time="00:30:00",
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {params.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --output {output.distances} 2>&1 | tee {log}
        """

def get_private_node_data(wildcards):
    node_data = [
        "builds/{build_name}/{segment}/scores.json",
    ]

    # Only try to annotate titer collections for HA.
    if wildcards.segment == "ha" and config["builds"][wildcards.build_name].get("enable_titer_models", False):
        for collection in config["builds"][wildcards.build_name]["titer_collections"]:
            node_data.append(f"builds/{wildcards.build_name}/{wildcards.segment}/titers_for_reference_viruses/{collection['name']}.json")
            node_data.append(f"builds/{wildcards.build_name}/{wildcards.segment}/haplotypes_without_references/{collection['name']}.json")
            node_data.append(rules.titer_tree_cross_immunities.output.cross_immunities.format(titer_collection=collection["name"], **wildcards))

    # Only annotate Welsh epitope distances for H3N2 HA builds.
    if "h3n2" in wildcards.build_name and wildcards.segment == "ha":
        node_data.append(f"builds/{wildcards.build_name}/{wildcards.segment}/welsh_epitope_distances.json")

    return node_data

rule export_private:
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
        node_data = _get_node_data_by_wildcards,
        private_node_data = get_private_node_data,
        auspice_config = lambda w: config['builds'][w.build_name]['auspice_config'],
        description = lambda w: config['builds'][w.build_name].get("description", "config/description.md"),
        lat_longs = config['lat-longs']
    output:
        auspice_json = "auspice/{build_name}_{segment}.json"
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
            --include-root-sequence-inline \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --minify-json \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
