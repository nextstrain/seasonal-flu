FREQUENCY_REGIONS = [
    'Africa',
    'Europe',
    'North America',
    'China',
    'South Asia',
    'Japan Korea',
    'Oceania',
    'South America',
    'Southeast Asia',
    'West Asia',
]

rule all_who:
    input:
        [
            "auspice-who/" + build.get("auspice_name", f"{build_name}_{{segment}}").format(segment=segment) + "_" + suffix + ".json"
            for build_name, build in config["builds"].items()
            for segment in config["segments"]
            for suffix in ['tree', 'meta', 'frequencies', 'titers', 'titer-tree-model', 'titer-sub-model', 'entropy', 'sequences']
        ],

def _get_file_by_auspice_name(wildcards):
    for build_name, build_params in config["builds"].items():
        for segment in config["segments"]:
            if build_params.get("auspice_name", f"{build_name}_{{segment}}").format(segment=segment) == wildcards.auspice_name:
                return f"auspice/{build_name}_{segment}_{wildcards.suffix}.json"

    return ""

rule rename_auspice_file:
    input:
        _get_file_by_auspice_name,
    output:
        "auspice-who/{auspice_name}_{suffix}.json",
    shell:
        """
        ln {input} {output}
        """

rule tree_frequencies:
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    params:
        min_date = lambda wildcards: config["builds"][wildcards.build_name].get("min_date"),
        pivot_interval = 1,
        regions = ['global'] + FREQUENCY_REGIONS,
        min_clade = 20,
    output:
        frequencies = build_dir + "/{build_name}/{segment}/tree_frequencies.json",
    conda: "environment.yaml"
    shell:
        """
        augur frequencies \
            --method diffusion \
            --include-internal-nodes \
            --tree {input.tree} \
            --regions {params.regions:q} \
            --metadata {input.metadata} \
            --pivot-interval {params.pivot_interval} \
            --minimal-clade-size {params.min_clade} \
            --min-date {params.min_date} \
            --output {output}
        """

rule filter_translations_by_region:
    input:
        translations=build_dir + "/{build_name}/{segment}/translations.done",
        metadata = build_dir + "/{build_name}/metadata.tsv",
        exclude = lambda wildcards: config["builds"][wildcards.build_name]["exclude"],
    output:
        translations = build_dir + "/{build_name}/{segment}/translations_by_region/{region}/{gene}.fasta",
    params:
        translations = build_dir + "/{build_name}/{segment}/nextalign/masked.gene.{gene}.fasta",
        min_date = lambda wildcards: config["builds"][wildcards.build_name].get("min_date"),
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        augur filter \
            --sequences {params.translations} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --query "region == '{wildcards.region}'" \
            --min-date {params.min_date} \
            --output-sequences {output.translations:q}
        """

def region_translations(wildcards):
    return [
        f"{build_dir}/{wildcards.build_name}/{wildcards.segment}/translations_by_region/{wildcards.region}/{gene}.fasta"
        for gene in GENES[wildcards.segment]
    ]

rule complete_mutation_frequencies_by_region:
    input:
        metadata = build_dir + "/{build_name}/metadata.tsv",
        alignment = region_translations,
    params:
        genes = lambda w: ','.join(GENES[w.segment]),
        min_date = lambda wildcards: config["builds"][wildcards.build_name].get("min_date"),
        min_freq = 0.003,
        pivot_interval = 1,
        stiffness = 20,
        inertia = 0.2,
    output:
        mut_freq = build_dir + "/{build_name}/{segment}/mutation_frequencies/{region}.json"
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/mutation_frequencies_{build_name}_{segment}_{region}.txt"
    log:
        "logs/mutation_frequencies_{build_name}_{segment}_{region}.txt"
    resources:
        mem_mb=4000,
    shell:
        """
        augur frequencies \
            --method diffusion \
            --alignments {input.alignment} \
            --metadata {input.metadata} \
            --gene-names {params.genes:q} \
            --pivot-interval {params.pivot_interval} \
            --stiffness {params.stiffness} \
            --inertia {params.inertia} \
            --ignore-char X \
            --min-date {params.min_date} \
            --minimal-frequency {params.min_freq} \
            --output {output.mut_freq:q} &> {log:q}
        """

rule global_mutation_frequencies:
    input:
        frequencies = expand(build_dir + "/{{build_name}}/{{segment}}/mutation_frequencies/{region}.json", region = FREQUENCY_REGIONS),
        tree_freq = rules.tree_frequencies.output,
    params:
        regions = FREQUENCY_REGIONS
    output:
        auspice="auspice/{build_name}_{segment}_frequencies.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/global_frequencies.py --region-frequencies {input.frequencies:q} \
                                              --tree-frequencies {input.tree_freq} \
                                              --regions {params.regions:q} \
                                              --output-auspice {output.auspice}
        """

rule scores:
    input:
        metadata = build_dir + "/{build_name}/metadata.tsv",
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
    output:
        node_data = build_dir + "/{build_name}/{segment}/scores.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/scores.py  --metadata {input.metadata} \
                                  --tree {input.tree} \
                                  --output {output}
        """

rule export_titers:
    input:
        sub = "builds/{build_name}/{segment}/titers-sub-model/titers.json",
        tree = "builds/{build_name}/{segment}/titers-tree-model/titers.json",
    output:
        raw = "auspice/{build_name}_{segment}_titers.json",
        tree = "auspice/{build_name}_{segment}_titer-tree-model.json",
        sub = "auspice/{build_name}_{segment}_titer-sub-model.json",
    run:
        import json
        with open(input.sub) as fh:
            sub = json.load(fh)

        with open(output.sub, 'wt') as sub_file:
            json.dump({'avidity': sub['avidity'],
                       'potency': sub['potency'],
                       'substitution': sub['substitution']},
                      sub_file, indent=1)

        with open(output.raw, 'wt') as raw_file:
            json.dump(sub['titers'], raw_file, indent=1)

        with open(input.tree) as fh:
            tree = json.load(fh)

        with open(output.tree, 'wt') as tree_file:
            json.dump({'avidity': tree['avidity'],
                       'potency': tree['potency'],
                       'dTiter': {k:v['dTiter'] for k,v in tree['nodes'].items()}},
                      tree_file, indent=1)

rule export_entropy:
    input:
        aln = rules.align.output.alignment,
        gene_map = lambda w: config['builds'][w.build_name]['annotation'],
    params:
        genes = lambda w: ','.join(GENES[w.segment]),
    output:
        "auspice/{build_name}_{segment}_entropy.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        '''
        python3 scripts/entropy.py --alignment {input.aln} \
                --genes {params.genes} \
                --gene_map {input.gene_map} \
                --output {output}
        '''

rule export_sequence_json:
    input:
        aln = rules.ancestral.output.node_data,
        tree = rules.refine.output.tree,
        aa_seqs = aggregate_translations,
    params:
        genes = lambda w: GENES[w.segment]
    output:
        "auspice/{build_name}_{segment}_sequences.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        '''
        python3 scripts/sequence_export.py --alignment {input.aln} \
                --genes {params.genes} \
                --tree {input.tree} \
                --translations {input.aa_seqs} \
                --output {output}
        '''

def _get_node_data_for_report_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.annotate_epiweeks.output.node_data,
        rules.annotate_recency_of_submissions.output.node_data,
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.clades.output.node_data,
        rules.traits.output.node_data,
        rules.scores.output.node_data,
    ]

    # Only request a distance file for builds that have mask configurations
    # defined.
    if _get_build_distance_map_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)

    if config["builds"][wildcards.build_name].get('enable_titer_models', False) and wildcards.segment == 'ha':
        for collection in config["builds"][wildcards.build_name]["titer_collections"]:
            inputs.append(rules.titers_sub.output.titers_model.format(titer_collection=collection["name"], **wildcards))
            inputs.append(rules.titers_tree.output.titers_model.format(titer_collection=collection["name"], **wildcards))

    if config["builds"][wildcards.build_name].get('enable_glycosylation', False) and wildcards.segment in ['ha', 'na']:
        inputs.append(rules.glyc.output.glyc)

    if config["builds"][wildcards.build_name].get('enable_lbi', False) and wildcards.segment in ['ha', 'na']:
        inputs.append(rules.lbi.output.lbi)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule export_who:
    input:
        tree = build_dir + "/{build_name}/{segment}/tree.nwk",
        metadata = build_dir + "/{build_name}/metadata.tsv",
        auspice_config = lambda w: config['builds'][w.build_name]['auspice_config'],
        node_data = _get_node_data_for_report_export,
        colors = "config/colors.tsv",
    output:
        tree = "auspice/{build_name}_{segment}_tree.json",
        meta = "auspice/{build_name}_{segment}_meta.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        augur export v1 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --output-tree {output.tree} \
            --output-meta {output.meta} \
            --minify-json
        """
