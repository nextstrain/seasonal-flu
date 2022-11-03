ruleorder: export_private > export

rule all_antigenic_plots:
    input:
        expand("builds/{build_name}/{segment}/antigenic_distances_between_strains.pdf", build_name=list(config["builds"].keys()), segment=["ha"])

rule plot_antigenic_distances_between_strains:
    input:
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains.tsv",
        clades=lambda wildcards: f"config/clades_for_titer_plots_{config['builds'][wildcards.build_name]['lineage']}.txt",
        references=lambda wildcards: f"config/references_for_titer_plots_{config['builds'][wildcards.build_name]['lineage']}.txt",
        colors="config/colors_for_titer_plots.tsv",
    output:
        plot="builds/{build_name}/{segment}/antigenic_distances_between_strains.pdf",
    benchmark:
        "benchmarks/plot_antigenic_distances_between_strains_{build_name}_{segment}.txt"
    log:
        "logs/plot_antigenic_distances_between_strains_{build_name}_{segment}.txt"
    params:
        min_test_date=2021.0,
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/plot_antigenic_distances_between_strains.py \
            --antigenic-distances {input.distances} \
            --min-test-date {params.min_test_date} \
            --clades {input.clades} \
            --references {input.references} \
            --colors {input.colors} \
            --plot-raw-data \
            --output {output.plot} 2>&1 | tee {log}
        """

rule annotate_titer_counts_for_test_viruses:
    input:
        tree="builds/{build_name}/{segment}/tree.nwk",
        titers="builds/{build_name}/titers.tsv",
    output:
        titer_counts="builds/{build_name}/{segment}/titers_for_test_viruses.json",
    benchmark:
        "benchmarks/annotate_titer_counts_for_test_viruses_{build_name}_{segment}.txt"
    log:
        "logs/annotate_titer_counts_for_test_viruses_{build_name}_{segment}.txt"
    params:
        attribute_name="titers_for_test_viruses",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/annotate_titers_per_node.py \
            --tree {input.tree} \
            --titers {input.titers} \
            --attribute-name {params.attribute_name} \
            --include-internal-nodes \
            --use-categorical-ranges \
            --output {output.titer_counts} 2>&1 | tee {log}
        """

rule annotate_titer_counts_for_reference_viruses:
    input:
        tree="builds/{build_name}/{segment}/tree.nwk",
        titers="builds/{build_name}/titers.tsv",
    output:
        titer_counts="builds/{build_name}/{segment}/titers_for_reference_viruses.json",
    benchmark:
        "benchmarks/annotate_titer_counts_for_reference_viruses_{build_name}_{segment}.txt"
    log:
        "logs/annotate_titer_counts_for_reference_viruses_{build_name}_{segment}.txt"
    params:
        attribute_name="titers_for_reference_viruses",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/annotate_titers_per_node.py \
            --tree {input.tree} \
            --titers {input.titers} \
            --attribute-name {params.attribute_name} \
            --use-references \
            --output {output.titer_counts} 2>&1 | tee {log}
        """

rule summarize_haplotype_titer_coverage:
    input:
        haplotypes="builds/{build_name}/{segment}/haplotypes.json",
        distances="builds/{build_name}/{segment}/antigenic_distances_between_strains.tsv",
        frequencies="builds/{build_name}/{segment}/tip-frequencies.json",
    output:
        table="builds/{build_name}/{segment}/haplotype_summary.tsv",
        node_data="builds/{build_name}/{segment}/haplotypes_without_references.json",
    benchmark:
        "benchmarks/summarize_haplotype_titer_coverage_{build_name}_{segment}.txt"
    log:
        "logs/summarize_haplotype_titer_coverage_{build_name}_{segment}.txt"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/summarize_haplotype_titer_coverage.py \
            --haplotypes {input.haplotypes} \
            --antigenic-distances {input.distances} \
            --frequencies {input.frequencies} \
            --output-table {output.table} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

rule export_private:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
        node_data = _get_node_data_by_wildcards,
        private_node_data = [
            "builds/{build_name}/{segment}/titers_for_test_viruses.json",
            "builds/{build_name}/{segment}/titers_for_reference_viruses.json",
            "builds/{build_name}/{segment}/haplotypes_without_references.json",
        ],
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
        export AUGUR_RECURSION_LIMIT=10000;
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} {input.private_node_data} \
            --include-root-sequence \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
