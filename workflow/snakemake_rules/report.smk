'''
Generate figure and statistics for the bi-annual influenza consultations
'''

localrules: get_nextclade_dataset
rule get_nextclade_dataset:
    output:
        nextclade = "report_statistics/nextclade_{lineage}/dataset/reference.fasta"
    conda: "../envs/nextstrain.yaml"
    params:
        nextclade_dir = "report_statistics/nextclade_{lineage}/dataset",
        dataset = "flu_{lineage}_ha"
    shell:
        """
        nextclade dataset get --name {params.dataset} --output-dir {params.nextclade_dir}
        """

rule nextclade_all:
    input:
        sequences = "data/{lineage}/ha.fasta",
        nextclade = "report_statistics/nextclade_{lineage}/dataset/reference.fasta"
    output:
        nextclade = "report_statistics/nextclade_{lineage}/nextclade.tsv"
    conda: "../envs/nextstrain.yaml"
    params:
        nextclade_dataset = "report_statistics/nextclade_{lineage}/dataset",
        nextclade_outdir = "report_statistics/nextclade_{lineage}/alignments",
    threads: 4
    shell:
        """
        nextclade run --verbosity=error --input-dataset {params.nextclade_dataset} -j {threads} --input-fasta {input.sequences} --output-tsv {output.nextclade} --output-dir {params.nextclade_outdir}
        """


rule clade_frequencies:
    input:
        nextclade = "report_statistics/nextclade_{lineage}/nextclade.tsv",
        metadata = "data/{lineage}/metadata.tsv"
    output:
        total_counts = "report_statistics/figures/{lineage}_total_counts.png"
        # clade_frequency_by_region = "report_statistics/figures/{lineage}_clade_by_region.png",
        # clade_frequency_by_clade = "report_statistics/figures/{lineage}_clade_by_clade.png"
    conda: "../envs/nextstrain.yaml"
    params:
        min_date = "2020-01-01"
    shell:
        """
        python3 scripts/graph_frequencies.py --metadata {input.metadata} --nextclade {input.nextclade} \
            --output-total-counts {output.total_counts}
        """
