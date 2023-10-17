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
    'ha': ['SigPep', 'HA1', 'HA2', 'HA'],
    'na': ['NA'],
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
        mem_mb=16000,
        time="0:30:00",
    shell:
        """
        nextalign run\
            -r {input.reference} \
            -m {input.annotation} \
            --genes {params.genes} \
            --jobs {threads} \
            --include-reference \
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
        method = config["tree"].get("method", "iqtree"),
        tree_builder_args = lambda wildcards: f"--tree-builder-args {config['tree']['tree-builder-args']}" if config["tree"].get("tree-builder-args") else "",
        override_default_args = lambda wildcards: "--override-default-args" if config["tree"].get("override_default_args", False) else "",
    threads: 8
    resources:
        mem_mb=16000,
        time="2:00:00",
    shell:
        """
        augur tree \
            --method {params.method} \
            --alignment {input.alignment} \
            {params.tree_builder_args} \
            {params.override_default_args} \
            --output {output.tree} \
            --nthreads {threads} 2>&1 | tee {log}
        """

rule prune_reference:
    input:
        tree = build_dir + "/{build_name}/{segment}/tree_raw.nwk",
        reference = lambda wildcards: config['builds'][wildcards.build_name]['reference'],
    output:
        tree = build_dir + "/{build_name}/{segment}/tree_without_outgroup.nwk",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/prune_reference.py \
            --tree {input.tree} \
            --reference {input.reference} \
            --output {output.tree}
        """

rule prune_outliers:
    input:
        tree = build_dir + "/{build_name}/{segment}/tree_without_outgroup.nwk",
        aln = build_dir+"/{build_name}/{segment}/aligned.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        tree = build_dir + "/{build_name}/{segment}/tree_without_outgroup_clean.nwk",
        outliers = build_dir + "/{build_name}/{segment}/outliers.tsv"
    params:
        keep_strains_argument=lambda wildcards: "--keep-strains " + config["builds"][wildcards.build_name]["include"] if "include" in config["builds"][wildcards.build_name] else "",
    shell:
        """
        python3 scripts/flag_outliers.py \
            --tree {input.tree:q} \
            --aln {input.aln} \
            --dates {input.metadata} \
            --cutoff 4.0 \
            {params.keep_strains_argument} \
            --output-tree {output.tree:q} --output-outliers {output.outliers} 2>&1 | tee {log}
        """

rule sanitize_trees:
    input:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree_without_outgroup_clean.nwk",  segment=config['segments'], build_dir=[build_dir]),
        alignments = expand("{build_dir}/{{build_name}}/{segment}/aligned.fasta",  segment=config['segments'], build_dir=[build_dir]),
    output:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree_common.nwk",  segment=config['segments'], build_dir=[build_dir]),
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
            --output {output.trees:q} 2>&1 | tee {log}
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
        alignment = build_dir + "/{build_name}/{segment}/aligned.fasta",
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
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/refine_{build_name}_{segment}.txt"
    log:
        "logs/refine_{build_name}_{segment}.txt"
    resources:
        mem_mb=16000,
        time="2:00:00",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --keep-root \
            --stochastic-resolve \
            --timetree \
            --use-fft \
            --no-covariance \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} 2>&1 | tee {log}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
        alignment = build_dir + "/{build_name}/{segment}/aligned.fasta",
        translations = aggregate_translations,
        reference =  lambda w: f"{config['builds'][w.build_name]['reference']}",
        annotation = lambda w: f"{config['builds'][w.build_name]['annotation']}",
    output:
        node_data = build_dir + "/{build_name}/{segment}/muts.json",
        translations_done = build_dir + "/{build_name}/{segment}/translations.done",
    params:
        inference = "joint",
        genes = lambda w: GENES[w.segment],
        input_translations = lambda w: build_dir + f"/{w.build_name}/{w.segment}/nextalign/masked.gene.%GENE.fasta",
        output_translations = lambda w: build_dir + f"/{w.build_name}/{w.segment}/nextalign/masked.gene.%GENE_withInternalNodes.fasta",
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
            --root-sequence {input.reference} \
            --annotation {input.annotation} \
            --genes {params.genes} \
            --translations "{params.input_translations}" \
            --output-node-data {output.node_data} \
            --output-translations "{params.output_translations}" \
            --inference {params.inference} 2>&1 | tee {log} && touch {output.translations_done}
        """

# rule variant_escape_prediction:
#     input:
#         translations=build_dir + "/{build_name}/ha/translations.done",
#     output:
#         node_data = build_dir + "/{build_name}/{segment}/dmsa-phenotype/{experiment}_variant_escape_prediction.json",
#         pred_data = build_dir + "/{build_name}/{segment}/dmsa-phenotype/{experiment}_variant_escape_prediction.csv",
#     log:
#         "logs/{build_name}/{segment}/{experiment}_variant_escape_prediction.txt"
#     params:
#         dms_wt_seq_id = lambda w: config["escape_models"][f"{w.experiment}"]["dms_wt_seq_id"],
#         mut_effects_df = lambda w: config["escape_models"][f"{w.experiment}"]["mut_effects_df"],
#         mut_effect_col = lambda w: config["escape_models"][f"{w.experiment}"]["mut_effect_col"],
#         mutation_col = lambda w: config["escape_models"][f"{w.experiment}"]["mutation_col"],
#         # alignment = lambda w: build_dir + f"/{w.build_name}/{w.segment}/nextalign/masked.gene.%GENE_withInternalNodes.fasta",
#         alignment = lambda w: build_dir + f"/{w.build_name}/{w.segment}/nextalign/masked.gene.HA_withInternalNodes.fasta"
#         #allow_aa_subs_at_unmeasured_sites = lambda w: config["escape_models"][f"{w.experiment}"]["allow_aa_subs_at_unmeasured_sites"],
#         #allow_unmeasured_aa_subs_at_these_sites = lambda w: config["escape_models"][f"{w.experiment}"]["allow_unmeasured_aa_subs_at_these_sites"],
#         #min_pred_pheno = lambda w: config["escape_models"][f"{w.experiment}"]["min_pred_pheno"],
#         #max_pred_pheno = lambda w: config["escape_models"][f"{w.experiment}"]["max_pred_pheno"],
#     conda:
#         "../../profiles/dmsa-phenotype/dmsa-pred/dmsa_env.yaml"
#     resources:
#         mem_mb=2000
#     shell:
#         """
#         python profiles/dmsa-phenotype/dmsa-pred/dmsa_pred.py phenotype-prediction \
#             --model-type additive \
#             --alignment {params.alignment} \
#             --dms-wt-seq-id {params.dms_wt_seq_id} \
#             --mut-effects-df {params.mut_effects_df} \
#             --mut-effect-col {params.mut_effect_col} \
#             --mutation-col {params.mutation_col} \
#             --experiment-label {wildcards.experiment} \
#             --output-json {output.node_data} \
#             --output-df {output.pred_data} 2>&1 | tee {log}
#         """

#--allow-aa-subs-at-unmeasured-sites {params.allow_aa_subs_at_unmeasured_sites} \
#--allow-unmeasured-aa-subs-at-these-sites {params.allow_unmeasured_aa_subs_at_these_sites} \
#--min-pred-pheno {params.min_pred_pheno} \
#--max-pred-pheno {params.max_pred_pheno} \

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
        muts = build_dir + "/{build_name}/ha/muts.json",
        clades = lambda wildcards: config["builds"][wildcards.build_name]["clades"],
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
            --mutations {input.muts} \
            --clades {input.clades} \
            --output {output.node_data} 2>&1 | tee {log}
        """

# Determine subclades for na and ha.
rule subclades:
    input:
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
        muts = build_dir + "/{build_name}/{segment}/muts.json",
        clades = lambda wildcards: config["builds"][wildcards.build_name].get("subclades"),
    output:
        node_data = build_dir + "/{build_name}/{segment}/subclades.json",
    params:
        membership_name = "subclade",
        label_name = "Subclade",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/subclades_{build_name}_{segment}.txt"
    log:
        "logs/subclades_{build_name}_{segment}.txt"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.muts} \
            --clades {input.clades} \
            --membership-name {params.membership_name} \
            --label-name {params.label_name} \
            --output {output.node_data} 2>&1 | tee {log}
        """

# Assign clade annotations to non-HA segments from HA.
rule import_clades:
    input:
        tree = build_dir + "/{build_name}/ha/tree.nwk",
        muts = build_dir + "/{build_name}/{segment}/muts.json",
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

rule annotate_haplotypes:
    input:
        tree=build_dir + "/{build_name}/ha/tree.nwk",
        translations=build_dir + "/{build_name}/ha/translations.done",
        clades=build_dir + "/{build_name}/ha/clades.json",
    output:
        haplotypes=build_dir + "/{build_name}/ha/haplotypes.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/annotate_haplotypes_{build_name}_ha.txt"
    log:
        "logs/annotate_haplotypes_{build_name}_ha.txt"
    params:
        min_tips=config.get("haplotypes", {}).get("min_tips", 5),
        alignment=build_dir + "/{build_name}/ha/nextalign/masked.gene.HA1_withInternalNodes.fasta",
    shell:
        """
        python3 scripts/annotate_haplotypes.py \
            --tree {input.tree} \
            --alignment {params.alignment} \
            --clades {input.clades} \
            --min-tips {params.min_tips} \
            --output-node-data {output.haplotypes} 2>&1 | tee {log}
        """

rule tip_frequencies:
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
    params:
        narrow_bandwidth = 1 / 12.0,
        wide_bandwidth = 3 / 12.0,
        proportion_wide = 0.0,
        min_date_arg = lambda w: f"--min-date {config['builds'][w.build_name]['min_date']}" if "min_date" in config["builds"].get(w.build_name, {}) else "",
        max_date = lambda w: config['builds'][w.build_name]['max_date'] if "max_date" in config["builds"].get(w.build_name, {}) else "0D",
        pivot_interval = 1
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
