'''
This file contains rules that run the canonical augur pipeline on a set of sequences
and generate a tree and a number of node_jsons

input:
 - builds/{build_name}/strains.txt
 - builds/{build_name}/metadata.tsv
 - builds/{build_name}/{segment}/sequences.fasta

output:
 - builds/{build_name}/{segment}/tree_raw.nwk
 - builds/{build_name}/{segment}/tree.nwk
 - builds/{build_name}/{segment}/branch-lengths.json
 - builds/{build_name}/{segment}/aa-muts.json
 - builds/{build_name}/{segment}/nt-muts.json
 - builds/{build_name}/{segment}/traits.json
'''


rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = "builds/{build_name}/{segment}/sequences.fasta",
        reference = "defaults/{lineage}/{segment}/reference.fasta",
        annotation = "defaults/{lineage}/{segment}/annotation.gff",
    output:
        alignment = "builds/{build_name}/{segment}/aligned.fasta"
    threads: 1
    resources:
        mem_mb=16000
    shell:
        """
        nextalign -r {input.reference} \
                  -m {input.annotation} \
                  -i {input.sequences} \
                  -o {output.alignment}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment,
        exclude_sites = "defaults/{lineage}/{segment}/exclude-sites.txt"
    output:
        tree = "builds/{build_name}/{segment}/tree_raw.nwk"
    threads: 8
    resources:
        mem_mb=16000
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads} \
            --exclude-sites {input.exclude_sites}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = "builds/{build_name}/metadata.tsv"
    output:
        tree = "builds/{build_name}/{segment}/tree.nwk",
        node_data = "builds/{build_name}/{segment}/branch-lengths.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = clock_rate,
        clock_std_dev = clock_std_dev
    conda: "environment.yaml"
    resources:
        mem_mb=16000
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --no-covariance \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "builds/{build_name}/{segment}/nt-muts.json"
    params:
        inference = "joint"
    conda: "environment.yaml"
    resources:
        mem_mb=4000
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "builds/{build_name}/{segment}/aa-muts.json",
    conda: "environment.yaml"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = "builds/{build_name}/metadata.tsv"
    output:
        node_data = "builds/{build_name}/{segment}/traits.json",
    params:
        columns = "region"
    conda: "environment.yaml"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """
