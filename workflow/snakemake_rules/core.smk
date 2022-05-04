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

def genes(segment):
    return {'ha':['SigPep','HA1', 'HA2'], 'na':['NA']}[segment]

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = build_dir + "/{build_name}/{segment}/sequences.fasta",
        reference =  lambda w: f"{config['builds'][w.build_name]['reference']}",
        annotation = lambda w: f"{config['builds'][w.build_name]['annotation']}"
    output:
        alignment = build_dir + "/{build_name}/{segment}/aligned.fasta"
    params:
        genes = lambda w: ','.join(genes(w.segment)),
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
            tree.collapse_all(lambda c: c.branch_length < 1e-4)

            tree.root_at_midpoint()
            tree.ladderize()
            Phylo.write(tree, output.trees[ti], 'newick')

rule treeknit:
    input:
        trees = rules.sanitize_trees.output.trees,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        trees = expand("{build_dir}/{{build_name}}/TreeKnit/tree_{segment}.resolved.nwk",  segment=config['segments'], build_dir=[build_dir]),
        mccs = build_dir + "/{build_name}/TreeKnit/MCCs.dat"
    params:
        treetime_tmpdir = build_dir + "/{build_name}/TreeTime_tmp",
        tmp_trees = expand("{build_dir}/{{build_name}}/TreeKnit/tree_{segment}.nwk",  segment=config['segments'], build_dir=[build_dir]),
        treeknit_tmpdir = build_dir + "/{build_name}/TreeKnit"
    shell:
        """
        treetime clock --tree {input.trees[0]} --dates {input.metadata} --sequence-length 1000 --outdir {params.treetime_tmpdir}
        mv {params.treetime_tmpdir}/rerooted.newick {params.tmp_trees[0]}

        treetime clock --tree {input.trees[1]} --dates {input.metadata} --sequence-length 1000 --outdir {params.treetime_tmpdir}
        mv {params.treetime_tmpdir}/rerooted.newick {params.tmp_trees[1]}

        treeknit {params.tmp_trees} --outdir {params.treeknit_tmpdir}
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
    return rate.get((config['builds'][w.build_name]["lineage"], w.segment), 0.001)

def clock_std_dev(w):
    return clock_rate(w)/5

rule treetime_arg:
    message:
        """
        Refining tree
          - estimate timetree
        """
    input:
        trees = rules.treeknit.output.trees,
        mccs = rules.treeknit.output.mccs,
        alignments = expand("{build_dir}/{{build_name}}/{segment}/aligned.fasta",  segment=config['segments'], build_dir=[build_dir]),
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        directory(build_dir + "/{build_name}/treetime_arg")
    params:
        coalescent = "const",
        date_inference = "always"
    conda: "environment.yaml"
    resources:
        mem_mb=16000
    shell:
        """
        treetime arg --trees {input.trees} \
                      --mccs {input.mccs} \
                      --alignments {input.alignments} \
                      --confidence --clock-std-dev 0.001 \
                      --outdir  {output} \
		      --time-marginal {params.date_inference} \
                      --dates {input.metadata} --keep-root --keep-polytomies
        """


rule refine:
    message:
        """
        Refining tree
          - estimate timetree
        """
    input:
        treetime_arg = rules.treetime_arg.output,
        mccs = rules.treeknit.output.mccs,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree.nwk",  segment=config['segments'], build_dir=[build_dir]),
        node_data = expand("{build_dir}/{{build_name}}/{segment}/branch-lengths.json",  segment=config['segments'], build_dir=[build_dir]),
    params:
        coalescent = "const",
        date_inference = "marginal"
    conda: "environment.yaml"
    resources:
        mem_mb=16000
    shell:
        """
        python scripts/make-branch-length-json.py --timetree {input.treetime_arg}/timetree_1.nexus \
                --divtree {input.treetime_arg}/divergence_tree_1.nexus \
                --dates {input.treetime_arg}/dates_1.tsv --mccs {input.mccs} \
                --output-tree {output.trees[0]} --output-node-data {output.node_data[0]}

        python scripts/make-branch-length-json.py --timetree {input.treetime_arg}/timetree_2.nexus \
                --divtree {input.treetime_arg}/divergence_tree_2.nexus \
                --dates {input.treetime_arg}/dates_2.tsv  --mccs {input.mccs}  \
                --output-tree {output.trees[1]} --output-node-data {output.node_data[1]}

        """


rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = build_dir+"/{build_name}/{segment}/tree.nwk",
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
	    --inference marginal \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = build_dir+"/{build_name}/{segment}/tree.nwk",
        reference =  lambda w: f"{config['builds'][w.build_name]['reference']}",
        annotation = lambda w: f"{config['builds'][w.build_name]['annotation']}"
    output:
        node_data = build_dir + "/{build_name}/{segment}/aa_muts.json",
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    params:
        translations = lambda w: [f"{build_dir}/{w.build_name}/{w.segment}/nextalign/sequences.gene.{gene}.fasta" for gene in genes(w.segment)],
        genes = lambda w: genes(w.segment)
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
        tree = build_dir+"/{build_name}/{segment}/tree.nwk",
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
        tree = build_dir+"/{build_name}/{segment}/tree.nwk",
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