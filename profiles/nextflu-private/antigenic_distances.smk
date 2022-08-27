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
