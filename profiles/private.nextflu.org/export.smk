rule all_who:
    input:
        [
            "auspice-who/" + build.get("auspice_name", f"{build_name}_{{segment}}").format(segment=segment) + "_" + suffix + ".json"
            for build_name, build in config["builds"].items()
            for segment in config["segments"]
            for suffix in ['tree', 'meta', 'titers', 'titer-tree-model', 'titer-sub-model', 'entropy', 'sequences']
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

# def region_translations(w):
#     genes = gene_names(w)
#     return ["results/full-aaseq-%s_%s_%s_%s_%s.fasta"%(g, w.region, w.lineage, w.segment, w.resolution)
#             for g in genes]

# for seg, genes in genes_to_translate.items():
#     rule:
#         input:
#             metadata = rules.parse.output.metadata,
#             sequences = rules.parse.output.sequences,
#             exclude = files.outliers,
#             reference = files.reference
#         params:
#             genes=genes,
#             region="{region}"
#         output:
#             alignments = expand("results/full-aaseq-{gene}_{{region}}_{{lineage}}_{{segment}}_{{resolution}}.fasta",
#                                 gene=genes)
#         conda: "../../workflow/envs/nextstrain.yaml"
#         shell:
#             """
#             python3 scripts/full_region_alignments.py  --sequences {input.sequences}\
#                                                  --metadata {input.metadata} \
#                                                  --exclude {input.exclude} \
#                                                  --genes {params.genes} \
#                                                  --region {params.region:q} \
#                                                  --resolution {wildcards.resolution} \
#                                                  --reference {input.reference} \
#                                                  --output {output.alignments:q}
#             """

# rule complete_mutation_frequencies_by_region:
#     input:
#         metadata = rules.parse.output.metadata,
#         alignment = region_translations
#     params:
#         genes = gene_names,
#         min_date = min_date,
#         max_date = max_date,
#         min_freq = 0.003,
#         pivot_interval = pivot_interval,
#         stiffness = 20,
#         inertia = 0.2
#     output:
#         mut_freq = "results/mutation_frequencies_{region}_{lineage}_{segment}_{resolution}.json"
#     conda: "../../workflow/envs/nextstrain.yaml"
#     benchmark:
#         "benchmarks/mutation_frequencies_{region}_{lineage}_{segment}_{resolution}.txt"
#     log:
#         "logs/mutation_frequencies_{region}_{lineage}_{segment}_{resolution}.txt"
#     resources:
#         mem_mb=4000,
#     run:
#         import os
#         alignments = [alignment
#             for alignment in input.alignment
#             if os.path.getsize(alignment) > 0]

#         genes = [filename.split('results/full-aaseq-',1)[1].split('_', 1)[0]
#             for filename in alignments]

#         # Make sure our filename splitting worked as expected and we got expected gene names
#         assert all(gene in params.genes for gene in genes), \
#             "Gene parsed from file path did not match any expected gene names."

#         if alignments:
#             shell("""
#                 augur frequencies --method diffusion \
#                                   --alignments {alignments:q} \
#                                   --metadata {input.metadata} \
#                                   --gene-names {genes:q} \
#                                   --pivot-interval {params.pivot_interval} \
#                                   --stiffness {params.stiffness} \
#                                   --inertia {params.inertia} \
#                                   --ignore-char X \
#                                   --min-date {params.min_date} \
#                                   --max-date {params.max_date} \
#                                   --minimal-frequency {params.min_freq} \
#                                   --output {output.mut_freq:q} &> {log:q}
#             """)
#         else:
#             # Create an empty JSON file if there are no alignments
#             shell("""
#                 echo {{}} > {output.mut_freq:q}
#             """)

# rule global_mutation_frequencies:
#     input:
#         frequencies = expand(build_dir + "/{{build_name}}/{{segment}}/mutation_frequencies/{region}.json", region = frequency_regions),
#         tree_freq = rules.tree_frequencies.output,
#     params:
#         regions = frequency_regions
#     output:
#         auspice="auspice/{build_name}_{segment}_frequencies.json",
#     conda: "../../workflow/envs/nextstrain.yaml"
#     shell:
#         '''
#         python3 scripts/global_frequencies.py --region-frequencies {input.frequencies:q} \
#                                               --tree-frequencies {input.tree_freq} \
#                                               --regions {params.regions:q} \
#                                               --output-auspice {output.auspice} \
#                                               --output-augur {output.augur}
#         '''

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
