rule all_report_outputs:
    input:
        derived_haplotypes=expand("tables/{lineage}/derived_haplotypes.md", lineage=["h1n1pdm", "h3n2", "vic"]),
        counts_by_clade=expand("tables/{lineage}/counts_of_recent_sequences_by_clade.md", lineage=["h1n1pdm", "h3n2", "vic"]),
        total_sample_count_by_lineage="figures/total-sample-count-by-lineage.png",

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
    params:
        lineages=["H1N1pdm", "H3N2", "Vic"],
        min_date="2023-04-01",
    shell:
        """
        python3 scripts/plot_counts_per_lineage.py \
            --metadata {input.h1n1pdm_metadata} {input.h3n2_metadata} {input.vic_metadata} \
            --lineages {params.lineages:q} \
            --min-date {params.min_date} \
            --output-count-by-lineage {output.total_sample_count_by_lineage} \
            --output-h1n1pdm-count {output.total_sample_count_h1n1pdm} \
            --output-h3n2-count {output.total_sample_count_h3n2} \
            --output-vic-count {output.total_sample_count_vic}
        """

rule annotate_recency_of_all_submissions:
    input:
        metadata = "data/{lineage}/metadata.tsv",
    output:
        node_data = "tables/{lineage}/recency.json",
    params:
        submission_date_field=config.get("submission_date_field"),
        date_bins=config.get("recency", {}).get("date_bins"),
        date_bin_labels=config.get("recency", {}).get("date_bin_labels"),
        upper_bin_label=config.get("recency", {}).get("upper_bin_label"),
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/all_recency_{lineage}.txt"
    log:
        "logs/all_recency_{lineage}.txt"
    shell:
        """
        python3 scripts/construct-recency-from-submission-date.py \
            --metadata {input.metadata} \
            --submission-date-field {params.submission_date_field} \
            --date-bins {params.date_bins} \
            --date-bin-labels {params.date_bin_labels:q} \
            --upper-bin-label {params.upper_bin_label} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule download_nextclade:
    output:
        nextclade="data/{lineage}/{segment}/nextclade.tsv.xz",
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{segment}/nextclade.tsv.xz"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} {output.nextclade}
        """

rule filter_nextclade_by_qc:
    input:
        nextclade="data/{lineage}/{segment}/nextclade.tsv.xz",
    output:
        nextclade="data/{lineage}/{segment}/nextclade_without_bad_qc.tsv",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        xz -c -d {input.nextclade} \
            | tsv-filter -H --str-ne "qc.overallStatus:bad" > {output.nextclade}
        """

rule count_recent_tips_by_clade:
    input:
        recency="tables/{lineage}/recency.json",
        clades="data/{lineage}/ha/nextclade_without_bad_qc.tsv",
    output:
        counts="tables/{lineage}/counts_of_recent_sequences_by_clade.md",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/count_recent_tips_by_clade.py \
            --recency {input.recency} \
            --clades {input.clades} \
            --output {output.counts}
        """

rule get_derived_haplotypes:
    input:
        nextclade="data/{lineage}/ha/nextclade_without_bad_qc.tsv",
    output:
        haplotypes="data/{lineage}/nextclade_with_derived_haplotypes.tsv",
    conda: "../../workflow/envs/nextstrain.yaml"
    params:
        genes=["HA1"],
    shell:
        """
        python3 scripts/add_derived_haplotypes.py \
            --nextclade {input.nextclade} \
            --genes {params.genes:q} \
            --strip-genes \
            --output {output.haplotypes}
        """

rule join_metadata_and_nextclade:
    input:
        metadata="data/{lineage}/metadata.tsv",
        nextclade="data/{lineage}/nextclade_with_derived_haplotypes.tsv",
    output:
        metadata="data/{lineage}/metadata_with_derived_haplotypes.tsv",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        tsv-join -H -f {input.nextclade} -a haplotype -k seqName -d strain {input.metadata} > {output.metadata}
        """

rule estimate_derived_haplotype_frequencies:
    input:
        metadata="data/{lineage}/metadata_with_derived_haplotypes.tsv",
    output:
        frequencies="tables/{lineage}/derived_haplotype_frequencies.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    params:
        narrow_bandwidth=1 / 12.0,
        min_date="16W",
        max_date=config.get("build_date", "4W"),
    shell:
        """
        python3 scripts/estimate_frequencies_from_metadata.py \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output.frequencies}
        """

rule summarize_derived_haplotypes:
    input:
        metadata="data/{lineage}/metadata_with_derived_haplotypes.tsv",
        frequencies="tables/{lineage}/derived_haplotype_frequencies.json",
        titers=lambda wildcards: [
            collection["data"]
            for collection in config["builds"][f"{wildcards.lineage}_2y_titers"]["titer_collections"]
            if "ferret" in collection["data"]
        ],
    output:
        table="tables/{lineage}/derived_haplotypes.tsv",
        markdown_table="tables/{lineage}/derived_haplotypes.md",
    conda: "../../workflow/envs/nextstrain.yaml"
    params:
        titer_names=lambda wildcards: [
            collection["name"]
            for collection in config["builds"][f"{wildcards.lineage}_2y_titers"]["titer_collections"]
            if "ferret" in collection["data"]
        ],
    shell:
        """
        python3 scripts/summarize_haplotypes.py \
            --metadata {input.metadata} \
            --frequencies {input.frequencies} \
            --titers {input.titers:q} \
            --titer-names {params.titer_names:q} \
            --output-table {output.table} \
            --output-markdown-table {output.markdown_table}
        """
