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
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model,
        rules.clades.output.clades,
        rules.traits.output.node_data,
        rules.lbi.output.lbi,
        files.vaccine_json
    ]

    if wildcards.lineage == "h3n2" and wildcards.segment == "ha" and wildcards.resolution == "2y":
        inputs.append(rules.forecast_tips.output.node_data)

    # Only request a distance file for builds that have distance map
    # configurations defined.
    if _get_build_distance_map_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

def _get_node_data_for_predictors(wildcards):
    """Return a list of node data files for fitness predictors
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model,
        rules.lbi.output.lbi,
        rules.delta_frequency.output.delta_frequency,
        rules.convert_translations_to_json.output.translations,
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
        model = "models/lbi-ne_star.json"
    output:
        node_data = "results/forecast_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
        frequencies = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_forecast-tip-frequencies.json"
    params:
        delta_months = _get_delta_months_to_forecast
    shell:
        """
        python3 flu-forecasting/src/forecast_model.py \
            --tip-attributes {input.attributes} \
            --distances {input.distances} \
            --frequencies {input.frequencies} \
            --model {input.model} \
            --delta-months {params.delta_months} \
            --output-node-data {output.node_data} \
            --output-frequencies {output.frequencies}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = files.auspice_config,
        node_data = _get_node_data_for_export
    output:
        auspice_json = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
            --minify-json
        """

def get_tip_frequencies(wildcards):
    if wildcards.lineage == "h3n2" and wildcards.segment == "ha" and wildcards.resolution == "2y":
        return "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_hi_forecast-tip-frequencies.json"
    else:
        return "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_hi_tip-frequencies.json"

rule simplify_auspice_names:
    input:
        main = "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_hi.json",
        frequencies = get_tip_frequencies
    output:
        main = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json"
    shell:
        '''
        mv {input.main} {output.main} &
        mv {input.frequencies} {output.frequencies} &
        '''

rule targets:
    input:
        main = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json"
    output:
        target = "targets/flu_seasonal_{lineage}_{segment}_{resolution}"
    shell:
        '''
        touch {output.target}
        '''

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
