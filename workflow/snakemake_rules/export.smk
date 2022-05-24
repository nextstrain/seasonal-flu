build_dir = config.get("build_dir", "builds")

def _get_node_data_by_wildcards(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.clades.output.node_data,
        rules.traits.output.node_data
    ]

    if config.get('titer-models',False):
        inputs.append(rules.titers_sub.output.titers_model)
        inputs.append(rules.titers_tree.output.titers_model)
    if config.get('glycosylation', False) and wildcards.segment in ['ha', 'na']:
        inputs.append(rules.glyc.output.glyc)
    if config.get('lbi', False) and wildcards.segment in ['ha', 'na']:
        inputs.append(rules.lbi.output.lbi)
    if config['builds'][wildcards.build_name].get('vaccines', False):
        inputs.append(config['builds'][wildcards.build_name].get('vaccines'))

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs

rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
        node_data = _get_node_data_by_wildcards,
        auspice_config = lambda w: config['builds'][w.build_name]['auspice_config'],
        lat_longs = config['lat-longs']
    output:
        auspice_json = "auspice/{build_name}_{segment}.json",
        root_sequence_json = "auspice/{build_name}_{segment}_root-sequence.json",
    log:
        "logs/export_{build_name}_{segment}.txt"
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000;
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --include-root-sequence \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
