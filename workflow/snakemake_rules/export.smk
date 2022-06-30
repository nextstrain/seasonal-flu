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
        rules.traits.output.node_data,
        rules.annotate_epiweeks.output.node_data,
        rules.annotate_recency_of_submissions.output.node_data,
    ]

    # Only request a distance file for builds that have distance map
    # configurations defined.
    if _get_build_distance_map_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)

    if config["builds"][wildcards.build_name].get('enable_titer_models', False) and wildcards.segment == 'ha':
        inputs.append(rules.titers_sub.output.titers_model)
        inputs.append(rules.titers_tree.output.titers_model)

    if config["builds"][wildcards.build_name].get('enable_glycosylation', False) and wildcards.segment in ['ha', 'na']:
        inputs.append(rules.glyc.output.glyc)

    if config["builds"][wildcards.build_name].get('enable_lbi', False) and wildcards.segment in ['ha', 'na']:
        inputs.append(rules.lbi.output.lbi)

    if config["builds"][wildcards.build_name].get('enable_forecasts', False) and config["builds"][wildcards.build_name].get("lineage") in ["h3n2"] and wildcards.segment in ["ha"]:
        wildcards_dict["model"] = config["fitness_model"]["best_model"]
        inputs.append(rules.titer_tree_cross_immunities.output.cross_immunities)
        inputs.append(rules.forecast_tips.output.node_data)
        inputs.append(rules.merge_weighted_distances_to_future.output.node_data)

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
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/export_{build_name}_{segment}.txt"
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

def _get_final_tip_frequencies(wildcards):
    if config["builds"][wildcards.build_name].get("enable_forecasts", False) and wildcards.segment == "ha":
        model = config["fitness_model"]["best_model"]
        return f"auspice/{wildcards.build_name}_{wildcards.segment}_{model}_forecast-tip-frequencies.json"
    else:
        return "builds/{build_name}/{segment}/tip-frequencies.json"

rule final_tip_frequencies:
    input:
        frequencies=_get_final_tip_frequencies,
    output:
        frequencies="auspice/{build_name}_{segment}_tip-frequencies.json",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        cp -f {input.frequencies} {output.frequencies}
        """
