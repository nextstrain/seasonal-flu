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

localrules: clades, sanitize_trees

build_dir = config.get("build_dir", "builds")

GENES = {
    'ha': ['SigPep', 'HA1', 'HA2'],
    'na': ['NA'],
}

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = build_dir + "/{build_name}/{segment}/sequences.fasta",
        reference =  lambda w: config['builds'][w.build_name]['reference'],
        annotation = lambda w: config['builds'][w.build_name]['annotation'],
    output:
        alignment = build_dir + "/{build_name}/{segment}/aligned.fasta"
    params:
        genes = lambda w: ','.join(GENES[w.segment]),
        outdir =  build_dir + "/{build_name}/{segment}/nextalign",
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
        tree = build_dir + "/{build_name}/{segment}/tree_raw.nwk"
    params:
        tree_builder_args = config["tree"]["tree-builder-args"]
    threads: 8
    resources:
        mem_mb=16000
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.tree_builder_args} \
            --output {output.tree} \
            --nthreads {threads}
        """

rule sanitize_trees:
    input:
        trees = lambda w: [f"{build_dir}/{w.build_name}/{segment}/tree_raw.nwk" for segment in config['segments']],
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree_common.nwk",  segment=config['segments'], build_dir=[build_dir])
    run:
        from Bio import Phylo

        trees = [Phylo.read(fname, 'newick') for fname in input.trees]
        common_leaves = set.intersection(*[set([x.name for x in tree.get_terminals()]) for tree in trees])
        for ti,tree in enumerate(trees):
            for leaf in set([x.name for x in tree.get_terminals()]).difference(common_leaves):
                tree.prune(leaf)

            tree.root_at_midpoint()
            tree.ladderize()
            Phylo.write(tree, output.trees[ti], 'newick')


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
    return rate.get((config['builds'][w.build_name]["lineage"], w.segment), 0.001)

def clock_std_dev(w):
    return clock_rate(w)/5

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
        tree = build_dir + "/{build_name}/{segment}/tree_common.nwk",
        alignment = rules.align.output,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
        node_data = build_dir + "/{build_name}/{segment}/branch-lengths.json"
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
        node_data = build_dir + "/{build_name}/{segment}/nt-muts.json"
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
        reference =  lambda w: f"{config['builds'][w.build_name]['reference']}",
        annotation = lambda w: f"{config['builds'][w.build_name]['annotation']}"
    output:
        node_data = build_dir + "/{build_name}/{segment}/aa_muts.json",
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    params:
        translations = lambda w: [f"{build_dir}/{w.build_name}/{w.segment}/nextalign/sequences.gene.{gene}.fasta" for gene in GENES[w.segment]],
        genes = lambda w: GENES[w.segment]
    conda: "environment.yaml"
    shell:
        """
        python3 scripts/translations_aamuts.py \
            --tree {input.tree} \
            --annotation {input.annotation} \
            --reference {input.reference} \
            --translations {params.translations:q} \
            --genes {params.genes} \
            --output {output.node_data} 2>&1 | tee {log} && touch {output.translations_done}
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        node_data = build_dir + "/{build_name}/{segment}/traits.json",
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

rule clades:
    message: "Annotating clades"
    input:
        tree = build_dir + "/{build_name}/ha/tree.nwk",
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = lambda w: config[build_dir + ""][w.build_name]["clades"] if w.segment=='ha' else f"{build_dir}/{w.build_name}/ha/clades.json"
    output:
        node_data = build_dir + "/{build_name}/{segment}/clades.json"
    run:
        if wildcards.segment == 'ha':
            shell("""
                augur clades \
                    --tree {input.tree} \
                    --mutations {input.nt_muts} {input.aa_muts} \
                    --clades {input.clades} \
                    --output {output.node_data}
            """)
        else:
            shell("""
                python3 scripts/import_tip_clades.py \
                    --tree {input.tree} \
                    --clades {input.clades} \
                    --output {output.node_data}
            """)

rule tip_frequencies:
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
        weights = "config/frequency_weights_by_region.json"
    params:
        narrow_bandwidth = 2 / 12.0,
        wide_bandwidth = 3 / 12.0,
        proportion_wide = 0.0,
        weight_attribute = "region",
        min_date = lambda w: config["builds"][w.build_name]["min-date"],
        max_date = lambda w: datetime.datetime.today().strftime("%Y-%m-%d"),
        pivot_interval = 2
    output:
        tip_freq = "auspice/{build_name}_{segment}_tip-frequencies.json"
    conda: "environment.yaml"
    shell:
        """
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --weights {input.weights} \
            --weights-attribute {params.weight_attribute} \
            --pivot-interval {params.pivot_interval} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output}
        """