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

build_dir = config.get("build_dir", "builds")

GENES = {
    'ha': ['SigPep', 'HA1', 'HA2'],
    'na': ['NA'],
    'pb1': ['PB1'],
    'pb2': ['PB2']
}

rule mask:
    input:
        sequences = build_dir + "/{build_name}/{segment}/sequences.fasta",
    output:
        sequences = build_dir + "/{build_name}/{segment}/masked.fasta",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/mask_{build_name}_{segment}.txt"
    log:
        "logs/mask_{build_name}_{segment}.txt"
    shell:
        """
        augur mask \
            --sequences {input.sequences} \
            --mask-invalid \
            --output {output.sequences} 2>&1 | tee {log}
        """

checkpoint align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = build_dir + "/{build_name}/{segment}/masked.fasta",
        reference =  lambda w: config['builds'][w.build_name]['reference'],
        annotation = lambda w: config['builds'][w.build_name]['annotation'],
    output:
        alignment = build_dir + "/{build_name}/{segment}/aligned.fasta",
        translations = directory(build_dir + "/{build_name}/{segment}/nextalign"),
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/align_{build_name}_{segment}.txt"
    log:
        "logs/align_{build_name}_{segment}.txt"
    params:
        genes = lambda w: ','.join(GENES[w.segment]),
    threads: 8
    resources:
        mem_mb=16000
    shell:
        """
        nextalign run\
            -r {input.reference} \
            -m {input.annotation} \
            --genes {params.genes} \
            --jobs {threads} \
            {input.sequences} \
            --output-fasta {output.alignment} \
            --output-translations "{output.translations}/masked.gene.{{gene}}.fasta" 2>&1 | tee {log}
        """

def aggregate_translations(wildcards):
    """The alignment rule produces multiple outputs that we cannot easily name prior
    to running the rule. The names of the outputs depend on the segment being
    aligned and Snakemake's `expand` function does not provide a way to lookup
    the gene names per segment. Instead, we use Snakemake's checkpoint
    functionality to determine the names of the output files after alignment
    runs. Downstream rules refer to this function to specify the translations
    for a given segment.

    """
    checkpoint_output = checkpoints.align.get(**wildcards).output.translations
    return expand(build_dir + "/{build_name}/{segment}/nextalign/masked.gene.{gene}.fasta",
                  build_name=wildcards.build_name,
                  segment=wildcards.segment,
                  gene=GENES[wildcards.segment])

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment,
    output:
        tree = build_dir + "/{build_name}/{segment}/tree_raw.nwk"
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/tree_{build_name}_{segment}.txt"
    log:
        "logs/tree_{build_name}_{segment}.txt"
    params:
        tree_builder_args = config["tree"]["tree-builder-args"],
        override_default_args = lambda wildcards: "--override-default-args" if config["tree"].get("override_default_args", False) else "",
    threads: 8
    resources:
        mem_mb=16000
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.tree_builder_args} \
            {params.override_default_args} \
            --output {output.tree} \
            --nthreads {threads} 2>&1 | tee {log}
        """

rule sanitize_trees:
    input:
        trees = lambda w: [f"{build_dir}/{w.build_name}/{segment}/tree_raw.nwk" for segment in config['segments']],
        alignments = lambda w: [f"{build_dir}/{w.build_name}/{segment}/aligned.fasta" for segment in config['segments']],
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree_common.nwk",  segment=config['segments'], build_dir=[build_dir])
    params:
        clock_filter = 4
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/sanitize_trees_{build_name}.txt"
    log:
        "logs/sanitize_trees_{build_name}.txt"
    shell:
        """
        python3 scripts/sanitize_trees.py \
            --trees {input.trees:q} \
            --alignments {input.alignments:q} \
            --metadata {input.metadata} \
            --clock-filter {params.clock_filter} \
            --output {output.trees:q} 2>&1 | tee {log}
        """

rule treeknit:
    input:
        trees = rules.sanitize_trees.output.trees,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        trees = expand("{build_dir}/{{build_name}}/TreeKnit/tree_common_{segment}_resolved.nwk",
                       segment=config['segments'], build_dir=[build_dir]),
        mccs = build_dir + "/{build_name}/TreeKnit/MCCs.json"
    params:
        treetime_tmpdir = build_dir + "/{build_name}/TreeTime_tmp",
        tmp_trees = expand("{build_dir}/{{build_name}}/TreeKnit/tree_{segment}.nwk",
                            segment=config['segments'], build_dir=[build_dir]),
        treeknit_tmpdir = build_dir + "/{build_name}/TreeKnit",
        clock_filter=4
    shell:
        """
        treeknit {input.trees} --outdir {params.treeknit_tmpdir}
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
        alignments = expand("{build_dir}/{{build_name}}/{segment}/aligned.fasta",
                            segment=config['segments'], build_dir=[build_dir]),
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        time_trees = expand("{build_dir}/{{build_name}}/treetime_arg/treetime_{segment}.nexus", 
                        segment=config['segments'], build_dir=[build_dir]),
        divergence_trees = expand("{build_dir}/{{build_name}}/treetime_arg/divergence_tree_{segment}.nexus", 
                        segment=config['segments'], build_dir=[build_dir]),
        dates = expand("{build_dir}/{{build_name}}/treetime_arg/dates_{segment}.tsv", 
                        segment=config['segments'], build_dir=[build_dir]),
    params:
        coalescent = "const",
        date_inference = "always",
        output_dir = build_dir + "/{build_name}/treetime_arg",
    resources:
        mem_mb=16000
    shell:
        """
        treetime arg --trees {input.trees} \
                      --mccs {input.mccs} \
                      --alignments {input.alignments} \
                      --confidence --clock-std-dev 0.001 \
                      --clock-filter 0 \
                      --outdir  {params.output_dir} \
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
        treetime_arg_time_trees = rules.treetime_arg.output.time_trees,
        treetime_arg_divergence_trees = rules.treetime_arg.output.divergence_trees,
        treetime_arg_dates = rules.treetime_arg.output.dates,
        mccs = rules.treeknit.output.mccs,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree.nwk",
                        segment=config['segments'], build_dir=[build_dir]),
        node_data = expand("{build_dir}/{{build_name}}/{segment}/branch-lengths.json",
                            segment=config['segments'], build_dir=[build_dir]),
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/refine_{build_name}.txt"
    log:
        "logs/refine_{build_name}.txt"
    resources:
        mem_mb=16000
    shell:
        """
        python scripts/make-branch-length-json.py --timetree {input.treetime_arg_time_trees[0]} \
                --divtree {input.treetime_arg_divergence_trees[0]} \
                --dates {input.treetime_arg_dates[0]} --mccs {input.mccs} \
                --output-tree {output.trees[0]} --output-node-data {output.node_data[0]}

        python scripts/make-branch-length-json.py --timetree {input.treetime_arg_time_trees[1]} \
                --divtree {input.treetime_arg_divergence_trees[1]} \
                --dates {input.treetime_arg_dates[1]} --mccs {input.mccs}  \
                --output-tree {output.trees[1]} --output-node-data {output.node_data[1]}

        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = build_dir+"/{build_name}/{segment}/tree.nwk",
        alignment = rules.align.output.alignment,
    output:
        node_data = build_dir + "/{build_name}/{segment}/nt-muts.json"
    params:
        inference = "joint"
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/ancestral_{build_name}_{segment}.txt"
    log:
        "logs/ancestral_{build_name}_{segment}.txt"
    resources:
        mem_mb=4000
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} 2>&1 | tee {log}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        translations = aggregate_translations,
        tree = build_dir+"/{build_name}/{segment}/tree.nwk",
        reference =  lambda w: f"{config['builds'][w.build_name]['reference']}",
        annotation = lambda w: f"{config['builds'][w.build_name]['annotation']}",
    output:
        node_data = build_dir + "/{build_name}/{segment}/aa_muts.json",
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    params:
        genes = lambda w: GENES[w.segment]
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/translate_{build_name}_{segment}.txt"
    log:
        "logs/translate_{build_name}_{segment}.txt"
    shell:
        """
        python3 scripts/translations_aamuts.py \
            --tree {input.tree} \
            --annotation {input.annotation} \
            --reference {input.reference} \
            --translations {input.translations:q} \
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
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/traits_{build_name}_{segment}.txt"
    log:
        "logs/traits_{build_name}_{segment}.txt"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence 2>&1 | tee {log}
        """

# Determine clades with HA mutations.
rule clades:
    input:
        tree = build_dir + "/{build_name}/ha/tree.nwk",
        nt_muts = build_dir + "/{build_name}/ha/nt-muts.json",
        aa_muts = build_dir + "/{build_name}/ha/aa_muts.json",
        clades = lambda wildcards: config[build_dir + ""][wildcards.build_name]["clades"],
    output:
        node_data = build_dir + "/{build_name}/ha/clades.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/clades_{build_name}.txt"
    log:
        "logs/clades_{build_name}.txt"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nt_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.node_data} 2>&1 | tee {log}
        """

# Assign clade annotations to non-HA segments from HA.
rule import_clades:
    input:
        tree = build_dir + "/{build_name}/ha/tree.nwk",
        nt_muts = build_dir + "/{build_name}/{segment}/nt-muts.json",
        aa_muts = build_dir + "/{build_name}/{segment}/aa_muts.json",
        clades = build_dir + "/{build_name}/ha/clades.json",
    output:
        node_data = build_dir + "/{build_name}/{segment}/clades.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/import_clades_{build_name}_{segment}.txt"
    log:
        "logs/import_clades_{build_name}_{segment}.txt"
    shell:
        """
        python3 scripts/import_tip_clades.py \
            --tree {input.tree} \
            --clades {input.clades} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule tip_frequencies:
    input:
        tree = build_dir+"/{build_name}/{segment}/tree.nwk",
        metadata = build_dir + "/{build_name}/metadata.tsv",
    params:
        narrow_bandwidth = 2 / 12.0,
        wide_bandwidth = 3 / 12.0,
        proportion_wide = 0.0,
        min_date_arg = lambda w: f"--min-date {config['builds'][w.build_name]['min_date']}" if "min_date" in config["builds"].get(w.build_name, {}) else "",
        max_date = lambda w: config['builds'][w.build_name]['max_date'] if "max_date" in config["builds"].get(w.build_name, {}) else "0D",
        pivot_interval = 2
    output:
        tip_freq = "builds/{build_name}/{segment}/tip-frequencies.json"
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/tip_frequencies_{build_name}_{segment}.txt"
    log:
        "logs/tip_frequencies_{build_name}_{segment}.txt"
    shell:
        """
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --pivot-interval {params.pivot_interval} \
            {params.min_date_arg} \
            --max-date {params.max_date} \
            --output {output} 2>&1 | tee {log}
        """

rule annotate_recency_of_submissions:
    input:
        metadata = "builds/{build_name}/metadata.tsv",
    output:
        node_data = "builds/{build_name}/recency.json",
    params:
        submission_date_field=config["submission_date_field"],
        date_bins=config["recency"]["date_bins"],
        date_bin_labels=config["recency"]["date_bin_labels"],
        upper_bin_label=config["recency"]["upper_bin_label"],
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/recency_{build_name}.txt"
    log:
        "logs/recency_{build_name}.txt"
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

rule annotate_epiweeks:
    input:
        metadata="builds/{build_name}/metadata.tsv",
    output:
        node_data="builds/{build_name}/epiweeks.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/annotate_epiweeks_{build_name}.txt"
    log:
        "logs/annotate_epiweeks_{build_name}.txt"
    shell:
        """
        python3 scripts/calculate_epiweek.py \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """
