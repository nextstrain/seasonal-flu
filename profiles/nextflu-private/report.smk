rule all_report_outputs:
    input:
        derived_haplotypes=expand("tables/{lineage}/derived_haplotypes.md", lineage=["h1n1pdm", "h3n2", "vic"]),
        counts_by_clade=expand("tables/{lineage}/counts_of_recent_sequences_by_clade.md", lineage=["h1n1pdm", "h3n2", "vic"]),
        total_sample_count_by_lineage="figures/total-sample-count-by-lineage.png",
        narrative="narratives/nextstrain-cdc.md",

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
        min_date="2024-09-01",
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

rule annotate_haplotypes_for_all_nextclade_data:
    input:
        nextclade="data/{lineage}/ha/nextclade_without_bad_qc.tsv",
        haplotypes="config/{lineage}/ha/emerging_haplotypes.tsv",
    output:
        haplotypes="data/{lineage}/ha/nextclade_with_haplotypes.tsv",
    params:
        clade_column="subclade",
        membership_name="emerging_haplotype",
    conda: "../../workflow/envs/nextstrain.yaml"
    log:
        "logs/annotate_haplotypes_for_all_nextclade_data_{lineage}_ha.txt"
    shell:
        r"""
        python scripts/assign_haplotypes.py \
            --substitutions {input.nextclade:q} \
            --haplotypes {input.haplotypes:q} \
            --clade-column {params.clade_column:q} \
            --haplotype-column-name {params.membership_name:q} \
            --use-clade-as-default-haplotype \
            --output-table {output.haplotypes:q} 2>&1 | tee {log}
        """

rule count_recent_tips_by_clade:
    input:
        recency="tables/{lineage}/recency.json",
        clades="data/{lineage}/ha/nextclade_with_haplotypes.tsv",
    output:
        counts="tables/{lineage}/counts_of_recent_sequences_by_clade.md",
    conda: "../../workflow/envs/nextstrain.yaml"
    params:
        clade_column="emerging_haplotype",
    shell:
        """
        python3 scripts/count_recent_tips_by_clade.py \
            --recency {input.recency} \
            --clades {input.clades} \
            --clade-column {params.clade_column:q} \
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
        clade_column="subclade",
        mutations_column="founderMuts['subclade'].aaSubstitutions",
        derived_haplotype_column="derived_haplotype",
    shell:
        """
        python3 scripts/add_derived_haplotypes.py \
            --nextclade {input.nextclade} \
            --genes {params.genes:q} \
            --strip-genes \
            --clade-column {params.clade_column:q} \
            --mutations-column {params.mutations_column:q} \
            --attribute-name {params.derived_haplotype_column:q} \
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
        tsv-join -H -f {input.nextclade} -a derived_haplotype -k seqName -d strain {input.metadata} > {output.metadata}
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
        haplotype_column="derived_haplotype",
    shell:
        """
        python3 scripts/summarize_haplotypes.py \
            --metadata {input.metadata} \
            --frequencies {input.frequencies} \
            --titers {input.titers:q} \
            --titer-names {params.titer_names:q} \
            --haplotype-column {params.haplotype_column:q} \
            --output-table {output.table} \
            --output-markdown-table {output.markdown_table}
        """

rule create_narrative:
    input:
        template="narratives_template/nextstrain-cdc.md.jinja",
        h1n1pdm_clade_counts="tables/h1n1pdm/counts_of_recent_sequences_by_clade.md",
        h1n1pdm_haplotypes="tables/h1n1pdm/derived_haplotypes.md",
        h3n2_clade_counts="tables/h3n2/counts_of_recent_sequences_by_clade.md",
        h3n2_haplotypes="tables/h3n2/derived_haplotypes.md",
        vic_clade_counts="tables/vic/counts_of_recent_sequences_by_clade.md",
        vic_haplotypes="tables/vic/derived_haplotypes.md",
    output:
        narrative="narratives/nextstrain-cdc.md",
    conda: "../../workflow/envs/nextstrain.yaml"
    params:
        date=config.get("build_date", datetime.date.today().strftime("%Y-%m-%d")),
        earliest_date_arg=lambda wildcards: f"--earliest-date '{config['earliest_reporting_date']}'" if "earliest_reporting_date" in config else "",
    shell:
        r"""
        python3 scripts/create_narrative.py \
            --template {input.template} \
            --date {params.date:q} \
            {params.earliest_date_arg} \
            --markdown-includes \
                h1n1pdm_clade_counts={input.h1n1pdm_clade_counts} \
                h1n1pdm_haplotypes={input.h1n1pdm_haplotypes} \
                h3n2_clade_counts={input.h3n2_clade_counts} \
                h3n2_haplotypes={input.h3n2_haplotypes} \
                vic_clade_counts={input.vic_clade_counts} \
                vic_haplotypes={input.vic_haplotypes} \
            --output {output.narrative}
        """
