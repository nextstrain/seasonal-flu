configfile: "config/config.json"

segments = ['ha', 'na']
lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam']
resolutions = ['6m', '2y','3y', '6y', '12y']

passages = ['cell']
centers = ['cdc']
assays = ['hi']

wildcard_constraints:
    lineage = "[A-Za-z0-9]{3,7}",
    segment = "[A-Za-z0-9]{2,3}",
    resolution = "[A-Za-z0-9]{2,3}",
    passage = "[a-z]{3,4}",
    center = "[a-z]{3,5}",
    assay = "[a-z]{2,3}"

localrules: download_all, simplify_auspice_names, targets, clean, clobber
include: "Snakefile_base"

rule all_live:
    input:
        auspice_main = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json",
                              lineage=lineages, segment=segments, resolution=resolutions),
        auspice_tip_frequencies = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json",
                              lineage=lineages, segment=segments, resolution=resolutions)

# separate rule for interaction with fauna
rule download_all:
    input:
        titers = expand("data/{center}_{lineage}_{passage}_{assay}_titers.tsv",
                         center=centers, lineage=lineages, passage=passages, assay=assays),
        sequences = expand("data/{lineage}_{segment}.fasta", lineage=lineages, segment=segments)


def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    inputs = [
        rules.annotate_epiweeks.output.node_data,
        rules.annotate_recency_of_submissions.output.node_data,
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model,
        rules.titer_tree_cross_immunities.output.cross_immunities,
        rules.clades.output.clades,
        rules.traits.output.node_data,
        rules.lbi.output.lbi,
        files.vaccine_json
    ]

    if wildcards.lineage == "h3n2" and wildcards.segment == "ha" and wildcards.resolution == "2y":
        wildcards_dict["model"] = config["fitness_model"]["best_model"]
        inputs.append(rules.forecast_tips.output.node_data)
        inputs.append(rules.merge_weighted_distances_to_future.output.node_data)

    # Only request a distance file for builds that have distance map
    # configurations defined.
    if _get_build_distance_map_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs

def _get_node_data_for_predictors(wildcards):
    """Return a list of node data files for fitness predictors
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model,
        rules.lbi.output.lbi,
        rules.convert_translations_to_json.output.translations,
        rules.titer_tree_cross_immunities.output.cross_immunities
    ]

    # Only request a distance file for builds that have distance map
    # configurations defined.
    if _get_build_distance_map_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule convert_node_data_to_table:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        node_data = _get_node_data_for_predictors
    output:
        table = "results/node-data_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.tsv"
    params:
        annotations = _get_annotations_for_node_data,
        excluded_fields_arg = _get_excluded_fields_arg
    conda: "environment.yaml"
    shell:
        """
        python3 flu-forecasting/scripts/node_data_to_table.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --jsons {input.node_data} \
            --output {output} \
            --annotations {params.annotations} \
            {params.excluded_fields_arg} \
        """


rule convert_frequencies_to_table:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.tip_frequencies.output.tip_freq
    output:
        table = "results/frequencies_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.tsv"
    conda: "environment.yaml"
    shell:
        """
        python3 scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --output {output}
        """


rule merge_node_data_and_frequencies:
    input:
        node_data = rules.convert_node_data_to_table.output.table,
        frequencies = rules.convert_frequencies_to_table.output.table
    output:
        table = "results/tip-attributes_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.tsv"
    run:
        node_data = pd.read_table(input.node_data)
        frequencies = pd.read_table(input.frequencies)
        df = node_data.merge(
            frequencies,
            how="inner",
            on=["strain", "is_terminal"]
        )
        df.to_csv(output.table, sep="\t", index=False, header=True)


rule target_distances:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table
    output:
        distances = "results/target-distances_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.tsv"
    params:
        delta_months = _get_delta_months_to_forecast
    conda: "environment.yaml"
    shell:
        """
        python3 flu-forecasting/scripts/calculate_target_distances.py \
            --tip-attributes {input.attributes} \
            --delta-months {params.delta_months} \
            --output {output}
        """


rule forecast_tips:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table,
        distances = rules.target_distances.output.distances,
        frequencies = rules.tip_frequencies.output.tip_freq,
        model = "models/{model}.json"
    output:
        node_data = "results/forecast_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_{model}.json",
        frequencies = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_{model}_forecast-tip-frequencies.json",
        table = "results/forecast_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_{model}.tsv"
    params:
        delta_months = _get_delta_months_to_forecast
    conda: "environment.yaml"
    shell:
        """
        python3 flu-forecasting/src/forecast_model.py \
            --tip-attributes {input.attributes} \
            --distances {input.distances} \
            --frequencies {input.frequencies} \
            --model {input.model} \
            --delta-months {params.delta_months} \
            --output-node-data {output.node_data} \
            --output-frequencies {output.frequencies} \
            --output-table {output.table}
        """


rule calculate_weighted_distance_to_future:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table,
        forecasts = rules.forecast_tips.output.table
    output:
        node_data = "results/weighted_distances_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_{model}.json",
    conda: "environment.yaml"
    shell:
        """
        python3 flu-forecasting/scripts/calculate_weighted_distances.py \
            --tip-attributes {input.attributes} \
            --forecasts {input.forecasts} \
            --distance-attribute-name weighted_distance_to_future_by_{wildcards.model} \
            --output {output.node_data}
        """


rule merge_weighted_distances_to_future:
    input:
        distances = expand("results/weighted_distances_{{center}}_{{lineage}}_{{segment}}_{{resolution}}_{{passage}}_{{assay}}_{model}.json", model=config["fitness_model"]["models"])
    output:
        node_data = "results/weighted_distances_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json"
    run:
        # Start with a single base JSON in node data format.
        with open(input.distances[0], "r") as fh:
            base_json = json.load(fh)

        # Update the base JSON with each subsequent model's distances to the future.
        for json_file in input.distances[1:]:
            with open(json_file, "r") as fh:
                model_json = json.load(fh)
                for strain, distances in model_json["nodes"].items():
                    base_json["nodes"][strain].update(distances)

        # Save merged data.
        with open(output.node_data, "w") as oh:
            json.dump(base_json, oh)


rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = files.auspice_config,
        node_data = _get_node_data_for_export,
        description = files.description,
        colors = files.colors,
        lat_longs = files.lat_longs,
    output:
        auspice_json = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
        root_sequence = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_root-sequence.json"
    conda: "environment.yaml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice_json} \
            --description {input.description} \
            --include-root-sequence \
            --minify-json
        """

def get_assay_by_lineage(wildcards):
    if wildcards.lineage == "h3n2":
        return "fra"
    else:
        return "hi"

def get_wildcards_dict(wildcards):
    wildcards_dict = dict(wildcards)
    assay = get_assay_by_lineage(wildcards)
    wildcards_dict["assay"] = assay

    return wildcards_dict

def get_main_auspice_json(wildcards):
    wildcards_dict = get_wildcards_dict(wildcards)
    return "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_{assay}.json".format(
        **wildcards_dict
    )

def get_root_sequence_json(wildcards):
    wildcards_dict = get_wildcards_dict(wildcards)
    return "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_{assay}_root-sequence.json".format(
        **wildcards_dict
    )

def get_tip_frequencies(wildcards):
    wildcards_dict = get_wildcards_dict(wildcards)

    if wildcards.lineage == "h3n2" and wildcards.segment == "ha" and wildcards.resolution == "2y":
        return "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_{assay}_{model}_forecast-tip-frequencies.json".format(
            model=config["fitness_model"]["best_model"],
            **wildcards_dict
        )
    else:
        return "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_{assay}_tip-frequencies.json".format(
            **wildcards_dict
        )

rule simplify_auspice_names:
    input:
        main = get_main_auspice_json,
        frequencies = get_tip_frequencies,
        root_sequence = get_root_sequence_json,
    output:
        main = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json",
        root_sequence = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_root-sequence.json"
    shell:
        '''
        mv {input.main} {output.main} &
        cp -f {input.frequencies} {output.frequencies} &
        mv {input.root_sequence} {output.root_sequence} &
        '''

rule targets:
    input:
        main = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json",
        root_sequence = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_root-sequence.json"
    output:
        target = "targets/flu_seasonal_{lineage}_{segment}_{resolution}"
    shell:
        '''
        touch {output.target}
        '''

rule forecasts:
    input: expand(rules.forecast_tips.output.frequencies, center="who", lineage="h3n2", segment="ha", resolution="2y", passage="cell", assay="hi", model=config["fitness_model"]["models"])

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "targets ",
        "auspice ",
        "auspice-who ",
        "logs"
    shell:
        "rm -rfv {params}"

rule clobber:
    message: "Removing directories: {params}"
    params:
        "results ",
        "targets ",
        "auspice ",
        "auspice-who ",
        "logs ",
        "data"
    shell:
        "rm -rfv {params}"
