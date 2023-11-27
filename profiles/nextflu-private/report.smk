rule plot_lineage_counts:
    input:
        h1n1pdm_metadata="data/h1n1pdm/metadata.tsv",
        h3n2_metadata="data/h3n2/metadata.tsv",
        vic_metadata="data/vic/metadata.tsv",
    output:
        total_sample_count_by_lineage="figures/total-sample-count-by-lineage.png",
        total_sample_count_h1n1pdm="figures/total-sample-count_h1n1pdm.png",
        total_sample_count_h3n2="figures/total-sample-count_h3n2.png",
        total_sample_count_vic="figures/total-sample-count_vic.png",
    conda: "../../workflow/envs/notebook.yaml"
    log:
        "logs/plot-counts-per-lineage.ipynb"
    notebook: "../../notebooks/plot-counts-per-lineage.py.ipynb"

rule all_counts_of_recent_tips_by_clade:
    input:
        counts=expand("builds/{build_name}/counts_of_recent_tips_by_clade.md", build_name=list(config["builds"].keys()))

rule count_recent_tips_by_clade:
    input:
        recency="builds/{build_name}/recency.json",
        clades="builds/{build_name}/ha/subclades.json",
    output:
        counts="builds/{build_name}/counts_of_recent_tips_by_clade.md",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/count_recent_tips_by_clade.py \
            --recency {input.recency} \
            --clades {input.clades} \
            --output {output.counts}
        """
