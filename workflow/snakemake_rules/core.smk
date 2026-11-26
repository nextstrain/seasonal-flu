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

import glob

build_dir = config.get("build_dir", "builds")

# coding regions translated and reconstructed per segment. For the internal segments those
# that the Nextclade datasets of all lineages have: PB1-F2 and PA-X are not in every one,
# and the second protein of MP is M2 in influenza A but BM2 in the B/Vic dataset
GENES = {
    'ha': ['SigPep', 'HA1', 'HA2'],
    'na': ['NA'],
    'pb2': ['PB2'],
    'pb1': ['PB1'],
    'pa': ['PA'],
    'np': ['NP'],
    'mp': ['M1'],
    'ns': ['NS1', 'NEP'],
}

def get_nextclade_dataset(w):
    """Name of the Nextclade dataset of a build and segment, or None."""
    build = config['builds'][w.build_name]
    if build.get('nextclade_dataset'):
        return build['nextclade_dataset'].format(segment=w.segment)
    return nextclade_dataset_by_lineage_and_segment.get(build['lineage'], {}).get(w.segment)


def _dataset_or_build_file(w, dataset_file, key):
    dataset = get_nextclade_dataset(w)
    if dataset:
        return f"nextclade_datasets/{dataset}/{dataset_file}"
    build = config['builds'][w.build_name]
    if key not in build:
        raise ValueError(f"Build {w.build_name} has no Nextclade dataset for {w.segment} and "
                         f"no {key}: give one of nextclade_dataset or {key}")
    return build[key].format(segment=w.segment)


def get_reference(w):
    """Reference sequence of a build and segment, from its Nextclade dataset."""
    return _dataset_or_build_file(w, "reference.fasta", "reference")


def get_annotation(w):
    """Genome annotation of a build and segment, from its Nextclade dataset."""
    return _dataset_or_build_file(w, "genome_annotation.gff3", "annotation")


rule get_nextclade_dataset:
    message:
        "Downloading the Nextclade dataset {wildcards.dataset}"
    output:
        reference = "nextclade_datasets/{dataset}/reference.fasta",
        annotation = "nextclade_datasets/{dataset}/genome_annotation.gff3",
    wildcard_constraints:
        dataset = r"nextstrain/flu/.+",
    shell:
        """
        nextclade3 dataset get --name {wildcards.dataset} --output-dir nextclade_datasets/{wildcards.dataset}
        """

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


def uses_tidy_tree(build_name, segment):
    """
    Whether to build this segment with the lineage guided tidy tree.

    Tidy tree needs a hierarchical lineage nomenclature with a guide tree and founder
    sequences, which exists for the HA and NA segments of the lineages listed in
    nomenclature_repo_by_lineage_and_segment. Everything else, in particular the
    internal segments and the lineages without a nomenclature such as yam, is built
    with augur tree.
    """
    if not config.get("tidy_tree", {}).get("enabled", True):
        return False
    lineage = config['builds'][build_name]['lineage']
    return nomenclature_url(lineage, segment, "founders.fasta") is not None


def get_founder_sequences(w):
    lineage = config['builds'][w.build_name]['lineage']
    segment = w.segment
    if not uses_tidy_tree(w.build_name, segment):
        return []
    return f"nomenclature/{lineage}/{segment}/founders.fasta"


def get_raw_tree(w):
    """The tree to prune outliers from, built either by tidy tree or by augur tree."""
    if uses_tidy_tree(w.build_name, w.segment):
        return build_dir + f"/{w.build_name}/{w.segment}/tree_tidy.nwk"
    return build_dir + f"/{w.build_name}/{w.segment}/tree_without_outgroup.nwk"

checkpoint align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = build_dir + "/{build_name}/{segment}/masked.fasta",
        founder_sequences= get_founder_sequences,
        reference = get_reference,
        annotation = get_annotation,
    output:
        alignment = build_dir + "/{build_name}/{segment}/aligned.fasta",
        translations = directory(build_dir + "/{build_name}/{segment}/translations"),
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/align_{build_name}_{segment}.txt"
    log:
        "logs/align_{build_name}_{segment}.txt"
    params:
        genes = lambda w: ','.join(GENES[w.segment]),
        # tidy tree needs the founders in the alignment; the other builds keep the
        # reference instead, which prune_reference removes again after the tree is built
        include_reference = lambda w: "" if uses_tidy_tree(w.build_name, w.segment) else "--include-reference",
    threads: 1
    resources:
        mem_mb=2000,
        time="0:30:00",
    shell:
        """
        nextclade3 run\
            {input.sequences} {input.founder_sequences} \
            -r {input.reference} \
            -m {input.annotation} \
            --gap-alignment-side right \
            --cds-selection {params.genes} \
            --jobs {threads} \
            {params.include_reference} \
            --output-fasta {output.alignment} \
            --output-translations "{output.translations}/{{cds}}.fasta" 2>&1 | tee {log}
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
    return expand(build_dir + "/{build_name}/{segment}/translations/{gene}.fasta",
                  build_name=wildcards.build_name,
                  segment=wildcards.segment,
                  gene=GENES[wildcards.segment])

rule download_nomenclature:
    message:
        "Fetching {wildcards.file} of the {wildcards.lineage} {wildcards.segment} nomenclature -> {output}"
    input:
        # lineages without a nomenclature repository, such as yam, keep their frozen
        # definitions in config
        local = lambda w: [] if nomenclature_url(w.lineage, w.segment, w.file) else f"config/{w.lineage}/{w.segment}/{w.file}",
    output:
        nomenclature_file = "nomenclature/{lineage}/{segment}/{file}"
    wildcard_constraints:
        file = r"clades\.tsv|subclades\.tsv|subclade\.nwk"
    params:
        url = lambda w: nomenclature_url(w.lineage, w.segment, w.file),
    shell:
        """
        if [ -n "{input.local}" ]; then
            cp {input.local} {output.nomenclature_file}
        else
            curl -fsSL -o {output.nomenclature_file} "{params.url}"
        fi
        """

rule download_founders:
    message:
        "Fetching lineage founders of the {wildcards.lineage} {wildcards.segment} nomenclature -> {output}"
    output:
        founders = "nomenclature/{lineage}/{segment}/founders.fasta"
    params:
        url = lambda w: nomenclature_url(w.lineage, w.segment, "founders.fasta"),
    shell:
        # The founders are aligned to the reference of the nomenclature. Remove the
        # gaps so that align places them against the reference of the build.
        """
        curl -fsSL {params.url} | sed '/^>/!s/-//g' > {output.founders}
        """


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
        method = config.get("tree", {}).get("method", "iqtree"),
        tree_builder_args = lambda wildcards: f"--tree-builder-args={config['tree']['tree-builder-args']}" if config.get("tree", {}).get("tree-builder-args") else "",
        override_default_args = lambda wildcards: "--override-default-args" if config.get("tree", {}).get("override_default_args", False) else "",
    threads: 8
    resources:
        mem_mb=4000,
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
        reference = get_reference,
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

def get_root_clade(w):
    """
    Root lineage of the tidy tree. By default tidy tree chooses the most specific
    lineage containing all sequences of the build; a lineage can be set instead with
    tidy_tree.root_lineage.{lineage}.{segment}.{resolution}, e.g. 2y, in which case
    sequences outside it are left out of the tree.
    """
    lineage = config['builds'][w.build_name]['lineage']
    resolution = f"{config['builds'][w.build_name].get('resolution', '6')}y"
    roots = config.get("tidy_tree", {}).get("root_lineage", {})
    return roots.get(lineage, {}).get(w.segment, {}).get(resolution, "auto")

def get_guide_tree_name(w):
    lineage = config['builds'][w.build_name]['lineage']
    segment = w.segment
    return f"nomenclature/{lineage}/{segment}/subclade.nwk"

# the tidy tree scripts, so that editing them reruns the rules that use them
TIDY_TREE_LIBRARY = sorted(glob.glob("scripts/tidytree/*.py"))

rule tidytree:
    message: "Building lineage guided tree"
    input:
        scripts = ["scripts/tidy_tree.py", *TIDY_TREE_LIBRARY],
        alignment = rules.align.output.alignment,
        guide_tree = get_guide_tree_name,
        # the lineage of each sequence is the nextclade call that ingest adds to the metadata
        metadata = build_dir + "/{build_name}/metadata.tsv",
    output:
        tree = build_dir + "/{build_name}/{segment}/tree_tidy.nwk",
        lineages = build_dir + "/{build_name}/{segment}/tidy_lineages.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/tidytree_{build_name}_{segment}.txt"
    log:
        "logs/tidytree_{build_name}_{segment}.txt"
    params:
        tree_builder_args = lambda wildcards: f"--tree-builder-args={config['tree']['tree-builder-args']}" if config.get("tree", {}).get("tree-builder-args") else "",
        override_default_args = lambda wildcards: "--override-default-args" if config.get("tree", {}).get("override_default_args", False) else "",
        root_lineage = lambda wildcards: "--root-lineage " + get_root_clade(wildcards),
        max_dropped_fraction = config.get("tidy_tree", {}).get("max_dropped_fraction", 0.05),
        unknown_lineage_policy = config.get("tidy_tree", {}).get("unknown_lineage_policy", "parent"),
    # independent subtrees are built at the same time
    threads: 4
    resources:
        mem_mb=2000,
        time="2:00:00",
    shell:
        """
        python3 scripts/tidy_tree.py \
            --alignment {input.alignment} \
            --guide-tree {input.guide_tree} \
            --lineage-assignments {input.metadata} \
            --seq-id-column strain \
            --lineage-column subclade_nextclade_{wildcards.segment} \
            --threads {threads} \
            --ignore-missing-founders \
            --unknown-lineage-policy {params.unknown_lineage_policy} \
            --max-dropped-fraction {params.max_dropped_fraction} \
            {params.root_lineage} \
            --metadata {input.metadata} \
            --output-lineages {output.lineages} \
            --output-tree {output.tree} 2>&1 | tee {log}
        """

rule prune_outliers:
    input:
        tree = get_raw_tree,
        aln = build_dir+"/{build_name}/{segment}/aligned.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        tree = build_dir + "/{build_name}/{segment}/tree_without_outliers_clean.nwk",
        outliers = build_dir + "/{build_name}/{segment}/outliers.tsv"
    params:
        keep_strains_argument=lambda wildcards: "--keep-strains " + config["builds"][wildcards.build_name]["include"] if "include" in config["builds"][wildcards.build_name] else "",
        cutoff=config.get("prune_outliers_z_score_cutoff", 4.0),
        keep_short_branches=lambda w: "--keep-short-branches" if uses_tidy_tree(w.build_name, w.segment) else "",
    resources:
        mem_mb=2000,
    shell:
        """
        python3 scripts/flag_outliers.py \
            --tree {input.tree:q} \
            --aln {input.aln} \
            --dates {input.metadata} \
            --cutoff {params.cutoff} \
            {params.keep_short_branches} \
            {params.keep_strains_argument} \
            --output-tree {output.tree:q} --output-outliers {output.outliers} 2>&1 | tee {log}
        """

rule sanitize_trees:
    input:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree_without_outliers_clean.nwk",  segment=config['segments'], build_dir=[build_dir]),
        alignments = expand("{build_dir}/{{build_name}}/{segment}/aligned.fasta",  segment=config['segments'], build_dir=[build_dir]),
    output:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree_common.nwk",  segment=config['segments'], build_dir=[build_dir]),
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/sanitize_trees_{build_name}.txt"
    log:
        "logs/sanitize_trees_{build_name}.txt"
    params:
        # tidy tree collapses short branches itself, with the lineages in view
        keep_short_branches = lambda w: [
            f"{build_dir}/{w.build_name}/{segment}/tree_without_outliers_clean.nwk"
            for segment in config['segments'] if uses_tidy_tree(w.build_name, segment)
        ],
    resources:
        mem_mb=2000,
    shell:
        """
        python3 scripts/sanitize_trees.py \
            --trees {input.trees:q} \
            --alignments {input.alignments:q} \
            --keep-short-branches {params.keep_short_branches:q} \
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

# TreeTime v1, the rewrite of TreeTime as a standalone executable, which the Nextstrain
# Docker image installs as `treetime1`; `treetime` is TreeTime v0, which augur imports
TREETIME = config.get("treetime", {}).get("executable", "treetime1")

def clock_filter_iqd():
    """Interquartile ranges beyond which the clock filter removes tips; 0 for no filter."""
    return config.get("refine", {}).get("clock_filter_iqd", 5)

rule refine:
    message:
        """
        Inferring the time tree with TreeTime
          - {params.coalescent} coalescent timescale
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation (0: no filter)
        """
    input:
        tree = build_dir + "/{build_name}/{segment}/tree_common.nwk",
        alignment = build_dir + "/{build_name}/{segment}/aligned.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        # the time tree with branch lengths in years, and the same tree in substitutions
        # per site (years times the clock rate), as augur refine writes it. TreeTime's
        # ancestral normalizes its inferred model to a rate of 1 per unit branch length and
        # cannot read the model of timetree, so it needs the rescaled tree; drop the
        # rescaling once ancestral fits the rate or reads the timetree model
        time_tree = build_dir + "/{build_name}/{segment}/tree_time.nwk",
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
        node_data = build_dir + "/{build_name}/{segment}/branch-lengths.json"
    params:
        coalescent = "opt",
        clock_filter_iqd = lambda wildcards: clock_filter_iqd(),
        # 0 switches the filter off, which needs the flag: TreeTime filters at 4 IQDs by default
        clock_filter = lambda wildcards: (f"--clock-filter {clock_filter_iqd()} --prune-filtered-tips"
                                          if clock_filter_iqd() else "--clock-filter 0"),
        clock_rate = clock_rate,
        clock_std_dev = clock_std_dev,
        # tidy tree has collapsed the short branches already. The polytomies are where tips
        # branch off at a founder; resolving them would move such tips into the subtree of
        # the lineage. Reoptimized branch lengths would let TreeTime collapse more branches.
        # Other trees get their polytomies resolved by the node dates.
        topology = lambda w: ("--no-collapse-short-branches"
                              if uses_tidy_tree(w.build_name, w.segment) else ""),
        divergence_units = lambda wildcards: config.get("refine", {}).get("divergence_units", "mutations-per-site"),
        treetime = TREETIME,
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/refine_{build_name}_{segment}.txt"
    log:
        "logs/refine_{build_name}_{segment}.txt"
    threads: 2
    resources:
        mem_mb=2000,
        time="2:00:00",
    shell:
        # the clock filter removes the flagged and the undated tips, as augur refine does
        """
        {params.treetime} timetree \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns strain \
            --date-column date \
            --keep-root \
            {params.topology} \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --confidence \
            {params.clock_filter} \
            --prune-undated-tips \
            --divergence-units {params.divergence_units} \
            --topology-order descendant-count \
            --seed 1 \
            --jobs {threads} \
            --no-progress \
            --output-tree-nwk {output.time_tree} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        python3 scripts/scale_time_tree.py \
            --tree {output.time_tree} \
            --node-data {output.node_data} \
            --output {output.tree} 2>&1 | tee -a {log}
        """

# TreeTime reads the coding regions from CDS rows, the annotations have gene rows
rule annotation_with_cds:
    input:
        annotation = get_annotation,
    output:
        annotation = build_dir + "/{build_name}/{segment}/annotation_cds.gff",
    shell:
        """
        python3 scripts/annotation_with_cds.py \
            --annotation {input.annotation} \
            --output {output.annotation}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations with TreeTime"
    input:
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
        alignment = build_dir + "/{build_name}/{segment}/aligned.fasta",
        translations = aggregate_translations,
        annotation = build_dir + "/{build_name}/{segment}/annotation_cds.gff",
    output:
        node_data = build_dir + "/{build_name}/{segment}/muts.json",
        translations_done = build_dir + "/{build_name}/{segment}/translations.done",
    params:
        genes = lambda w: GENES[w.segment],
        input_translations = lambda w: build_dir + f"/{w.build_name}/{w.segment}/translations/%GENE.fasta",
        output_translations = lambda w: build_dir + f"/{w.build_name}/{w.segment}/translations/%GENE_withInternalNodes.fasta",
        treetime = TREETIME,
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/ancestral_{build_name}_{segment}.txt"
    log:
        "logs/ancestral_{build_name}_{segment}.txt"
    threads: 2
    resources:
        mem_mb=2000
    shell:
        # the translations of the tips come from the alignment by nextclade. Ambiguous tip
        # states are resolved, which leaves changes to ambiguity codes out of the mutations,
        # as augur ancestral does. The mutations of the root are relative to the inferred
        # root sequence, which is the reference of the node data.
        """
        {params.treetime} ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --cds-selection {params.genes} \
            --translations "{params.input_translations}" \
            --impute-missing-data \
            --jobs {threads} \
            --no-progress \
            --output-node-data {output.node_data} \
            --output-translations "{params.output_translations}" 2>&1 | tee {log} && touch {output.translations_done}
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
        columns = config.get("traits", {}).get("columns", "region"),
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

# Determine subclades for na and ha, from their signature mutations.
rule subclades_from_mutations:
    input:
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
        muts = build_dir + "/{build_name}/{segment}/muts.json",
        clades = lambda wildcards: config["builds"][wildcards.build_name].get("subclades"),
    output:
        node_data = build_dir + "/{build_name}/{segment}/subclades_mutations.json",
    params:
        membership_name = "subclade",
        label_name = "Subclade",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/subclades_from_mutations_{build_name}_{segment}.txt"
    log:
        "logs/subclades_from_mutations_{build_name}_{segment}.txt"
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

# Determine subclades of a tidy tree from the lineages its tips were placed in. With
# tidy_tree.compare_subclades, the labels from signature mutations are computed alongside
# and the disagreements reported in the log.
rule subclades_from_tree:
    input:
        scripts = ["scripts/assign_subclades_from_tree.py", *TIDY_TREE_LIBRARY],
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
        lineages = build_dir + "/{build_name}/{segment}/tidy_lineages.tsv",
        guide_tree = get_guide_tree_name,
        compare = lambda w: (build_dir + f"/{w.build_name}/{w.segment}/subclades_mutations.json"
                             if config.get("tidy_tree", {}).get("compare_subclades", False) else []),
    output:
        node_data = build_dir + "/{build_name}/{segment}/subclades_tree.json",
    params:
        membership_name = "subclade",
        label_name = "Subclade",
        compare = lambda w, input: f"--compare {input.compare}" if input.compare else "",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/subclades_from_tree_{build_name}_{segment}.txt"
    log:
        "logs/subclades_from_tree_{build_name}_{segment}.txt"
    shell:
        """
        python3 scripts/assign_subclades_from_tree.py \
            --tree {input.tree} \
            --lineages {input.lineages} \
            --guide-tree {input.guide_tree} \
            --membership-name {params.membership_name} \
            --label-name {params.label_name} \
            {params.compare} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

rule subclades:
    input:
        node_data = lambda w: build_dir + (f"/{w.build_name}/{w.segment}/subclades_tree.json"
                                           if uses_tidy_tree(w.build_name, w.segment)
                                           else f"/{w.build_name}/{w.segment}/subclades_mutations.json"),
    output:
        node_data = build_dir + "/{build_name}/{segment}/subclades.json",
    shell:
        """
        cp {input.node_data} {output.node_data}
        """

rule annotate_derived_haplotypes:
    input:
        tree=build_dir + "/{build_name}/ha/tree.nwk",
        translations=build_dir + "/{build_name}/ha/translations.done",
        clades=build_dir + "/{build_name}/ha/subclades.json",
    output:
        haplotypes=build_dir + "/{build_name}/ha/derived_haplotypes.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/annotate_derived_haplotypes_{build_name}_ha.txt"
    log:
        "logs/annotate_derived_haplotypes_{build_name}_ha.txt"
    params:
        min_tips=config.get("haplotypes", {}).get("min_tips", 5),
        alignment=build_dir + "/{build_name}/ha/translations/HA1_withInternalNodes.fasta",
        clade_label_attribute="Subclade",
        clade_node_attribute="subclade",
    shell:
        """
        python3 scripts/annotate_derived_haplotypes.py \
            --tree {input.tree} \
            --alignment {params.alignment} \
            --clades {input.clades} \
            --clade-label-attribute {params.clade_label_attribute} \
            --clade-node-attribute {params.clade_node_attribute} \
            --min-tips {params.min_tips} \
            --output-node-data {output.haplotypes} 2>&1 | tee {log}
        """

rule tip_frequencies:
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
    params:
        narrow_bandwidth=lambda wildcards: config["builds"][wildcards.build_name].get("frequencies", {}).get("narrow_bandwidth", 1 / 20.0),
        wide_bandwidth=lambda wildcards: config["builds"][wildcards.build_name].get("frequencies", {}).get("wide_bandwidth", 3 / 12.0),
        proportion_wide=lambda wildcards: config["builds"][wildcards.build_name].get("frequencies", {}).get("proportion_wide", 0.0),
        min_date_arg = lambda w: f"--min-date {config['builds'][w.build_name]['min_date']}" if "min_date" in config["builds"].get(w.build_name, {}) else "",
        max_date = lambda w: config['builds'][w.build_name]['max_date'] if "max_date" in config["builds"].get(w.build_name, {}) else "0D",
        pivot_interval=lambda wildcards: config["builds"][wildcards.build_name].get("frequencies", {}).get("pivot_interval", 1),
        pivot_interval_units=lambda wildcards: config["builds"][wildcards.build_name].get("frequencies", {}).get("pivot_interval_units", "weeks"),
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
            --pivot-interval-units {params.pivot_interval_units} \
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
        submission_date_field=config.get("submission_date_field"),
        date_bins=config.get("recency", {}).get("date_bins"),
        date_bin_labels=config.get("recency", {}).get("date_bin_labels"),
        upper_bin_label=config.get("recency", {}).get("upper_bin_label"),
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

rule get_strains_in_alignment:
    input:
        alignment = "builds/{build_name}/{segment}/aligned.fasta",
    output:
        alignment_strains = "builds/{build_name}/{segment}/aligned_strains.txt",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        seqkit fx2tab -n -i {input.alignment} | sort -k 1,1 > {output.alignment_strains}
        """

rule get_shared_strains_in_alignments:
    input:
        alignment_strains = expand("builds/{{build_name}}/{segment}/aligned_strains.txt", segment=config["segments"]),
    output:
        shared_strains = "builds/{build_name}/shared_strains_in_alignment.txt",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/intersect_items.py \
            --items {input.alignment_strains:q} \
            --output {output.shared_strains}
        """

rule select_shared_strains_from_alignment_and_sort:
    input:
        shared_strains = "builds/{build_name}/shared_strains_in_alignment.txt",
        alignment = "builds/{build_name}/{segment}/aligned.fasta",
    output:
        alignment = "builds/{build_name}/{segment}/aligned.sorted.fasta",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        seqkit grep -f {input.shared_strains} {input.alignment} \
            | seqkit sort -n > {output.alignment}
        """

rule calculate_pairwise_distances:
    input:
        alignment = "builds/{build_name}/{segment}/aligned.sorted.fasta",
    output:
        distances = "builds/{build_name}/{segment}/distances.csv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/calculate_pairwise_distances_{build_name}_{segment}.txt"
    shell:
        """
        pathogen-distance \
            --alignment {input.alignment} \
            --output {output.distances}
        """

rule embed_with_tsne:
    input:
        alignments = expand("builds/{{build_name}}/{segment}/aligned.sorted.fasta", segment=config["segments"]),
        distances = expand("builds/{{build_name}}/{segment}/distances.csv", segment=config["segments"]),
    output:
        embedding = "builds/{build_name}/embed_tsne.csv",
    conda: "../envs/nextstrain.yaml"
    params:
        perplexity=config.get("embedding", {}).get("perplexity", 200),
    benchmark:
        "benchmarks/embed_with_tsne_{build_name}.txt"
    shell:
        """
        pathogen-embed \
            --alignment {input.alignments} \
            --distance-matrix {input.distances} \
            --output-dataframe {output.embedding} \
            t-sne \
                --perplexity {params.perplexity}
        """

rule cluster_tsne_embedding:
    input:
        embedding = "builds/{build_name}/embed_tsne.csv",
    output:
        clusters = "builds/{build_name}/cluster_embed_tsne.csv",
    conda: "../envs/nextstrain.yaml"
    params:
        label_attribute="tsne_cluster",
        distance_threshold=1.0,
    benchmark:
        "benchmarks/cluster_tsne_embedding_{build_name}.txt"
    shell:
        """
        pathogen-cluster \
            --embedding {input.embedding} \
            --label-attribute {params.label_attribute:q} \
            --distance-threshold {params.distance_threshold} \
            --output-dataframe {output.clusters}
        """

rule convert_embedding_clusters_to_node_data:
    input:
        clusters = "builds/{build_name}/cluster_embed_tsne.csv",
    output:
        node_data = "builds/{build_name}/cluster_embed_tsne.json",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/table_to_node_data.py \
            --table {input.clusters} \
            --output {output.node_data}
        """
