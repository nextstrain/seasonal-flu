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

def genes(segment):
    return {'ha':['HA'], 'na':['NA']}[segment]

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = "builds/{build_name}/{segment}/sequences.fasta",
        reference =  lambda w: f"{config['builds'][w.build_name]['reference']}",
        annotation = lambda w: f"{config['builds'][w.build_name]['annotation']}",
    output:
        alignment = "builds/{build_name}/{segment}/aligned.fasta"
    params:
        genes = lambda w: ','.join(genes(w.segment)),
        outdir =  "builds/{build_name}/{segment}/nextalign",
    threads: 1
    resources:
        mem_mb=16000
    shell:
        """
        nextalign -r {input.reference} \
                  -m {input.annotation} \
                  --genes {params.genes} \
                  -i {input.sequences} \
                  -o {output.alignment} \
                  --output-dir {params.outdir}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment,
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
            --nthreads {threads}
        """

def clock_rate(w):
    # these rates are from 12y runs on 2019-10-18
    rate = {
     ('h1n1pdm', 'ha'): 0.00329,
 	 ('h1n1pdm', 'na'): 0.00326,
 	 ('h1n1pdm', 'np'): 0.00221,
 	 ('h1n1pdm', 'pa'): 0.00217,
	 ('h1n1pdm', 'pb1'): 0.00205,
 	 ('h1n1pdm', 'pb2'): 0.00277,
 	 ('h3n2', 'ha'): 0.00382,
 	 ('h3n2', 'na'): 0.00267,
	 ('h3n2', 'np'): 0.00157,
 	 ('h3n2', 'pa'): 0.00178,
 	 ('h3n2', 'pb1'): 0.00139,
 	 ('h3n2', 'pb2'): 0.00218,
 	 ('vic', 'ha'): 0.00145,
 	 ('vic', 'na'): 0.00133,
 	 ('vic', 'np'): 0.00132,
 	 ('vic', 'pa'): 0.00178,
 	 ('vic', 'pb1'): 0.00114,
 	 ('vic', 'pb2'): 0.00106,
 	 ('yam', 'ha'): 0.00176,
 	 ('yam', 'na'): 0.00177,
 	 ('yam', 'np'): 0.00133,
 	 ('yam', 'pa'): 0.00112,
 	 ('yam', 'pb1'): 0.00092,
 	 ('yam', 'pb2'): 0.00113}
    return rate.get((w.lineage, w.segment), 0.001)

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
        clock_std_dev = 0.0001
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

# rule translate:
#     message: "Translating amino acid sequences"
#     input:
#         tree = rules.refine.output.tree,
#         node_data = rules.ancestral.output.node_data,
#         reference = files.reference
#     output:
#         node_data = "builds/{build_name}/{segment}/aa-muts.json",
#     conda: "environment.yaml"
#     shell:
#         """
#         augur translate \
#             --tree {input.tree} \
#             --ancestral-sequences {input.node_data} \
#             --reference-sequence {input.reference} \
#             --output {output.node_data} \
#         """

# rule traits:
#     message:
#         """
#         Inferring ancestral traits for {params.columns!s}
#         """
#     input:
#         tree = rules.refine.output.tree,
#         metadata = "builds/{build_name}/metadata.tsv"
#     output:
#         node_data = "builds/{build_name}/{segment}/traits.json",
#     params:
#         columns = "region"
#     conda: "environment.yaml"
#     shell:
#         """
#         augur traits \
#             --tree {input.tree} \
#             --metadata {input.metadata} \
#             --output {output.node_data} \
#             --columns {params.columns} \
#             --confidence
#         """
