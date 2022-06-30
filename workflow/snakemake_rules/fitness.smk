
build_dir = config.get("build_dir", "builds")
localrules: glyc, lbi

glyc_gene = {'ha':'HA1', 'na':'NA'}

rule glyc:
    input:
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    output:
        glyc = build_dir + "/{build_name}/{segment}/glyc.json"
    params:
        alignment = lambda w: f"{build_dir}/{w.build_name}/{w.segment}/nextalign/masked.gene.{glyc_gene.get(w.segment)}_withInternalNodes.fasta",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/glyc_{build_name}_{segment}.txt"
    log:
        "logs/glyc_{build_name}_{segment}.txt"
    shell:
        """
        python3 scripts/glyc.py \
            --tree {input.tree} \
            --alignment {params.alignment} \
            --output {output.glyc} 2>&1 | tee {log}
        """

rule lbi:
    message: "Calculating LBI"
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = 0.5,
        window = 0.5,
        names = "lbi"
    output:
        lbi =  build_dir + "/{build_name}/{segment}/lbi.json"
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/lbi_{build_name}_{segment}.txt"
    log:
        "logs/lbi_{build_name}_{segment}.txt"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --output {output} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window} 2>&1 | tee {log}
        """

# Rules for forecasting models.

rule convert_translations_to_json:
    input:
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done",
        translations = lambda w: [f"{build_dir}/{w.build_name}/{w.segment}/nextalign/masked.gene.{gene}_withInternalNodes.fasta" for gene in GENES[w.segment]],
    output:
        translations = "builds/{build_name}/{segment}/aa-seq.json",
    params:
        gene_names = lambda w: GENES[w.segment],
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/convert_translations_to_json_{build_name}_{segment}.txt"
    log:
        "logs/convert_translations_to_json_{build_name}_{segment}.txt"
    shell:
        """
        python3 flu-forecasting/scripts/convert_translations_to_json.py \
            --tree {input.tree} \
            --alignment {input.translations} \
            --gene-names {params.gene_names} \
            --output {output.translations} 2>&1 | tee {log}
        """

rule pairwise_titer_tree_distances:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.tip_frequencies.output.tip_freq,
        model = rules.titers_tree.output.titers_model,
        date_annotations = rules.refine.output.node_data,
    output:
        distances = "builds/{build_name}/{segment}/pairwise-titer-tree-distances.json",
    params:
        attribute_names = "cTiter_pairwise",
        months_back_for_current_samples = config.get("fitness_model", {}).get("months_back_for_current_samples"),
        years_back_to_compare = config.get("fitness_model", {}).get("max_years_for_distances"),
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/pairwise_titer_tree_distances_{build_name}_{segment}.txt"
    log:
        "logs/pairwise_titer_tree_distances_{build_name}_{segment}.txt"
    resources:
        mem_mb=4000
    shell:
        """
        python3 flu-forecasting/scripts/pairwise_titer_tree_distances.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --model {input.model} \
            --attribute-name {params.attribute_names} \
            --date-annotations {input.date_annotations} \
            --months-back-for-current-samples {params.months_back_for_current_samples} \
            --years-back-to-compare {params.years_back_to_compare} \
            --output {output} 2>&1 | tee {log}
        """

rule titer_tree_cross_immunities:
    input:
        frequencies = rules.tip_frequencies.output.tip_freq,
        distances = rules.pairwise_titer_tree_distances.output.distances,
        date_annotations = rules.refine.output.node_data,
    output:
        cross_immunities = "builds/{build_name}/{segment}/titer-tree-cross-immunity.json",
    params:
        distance_attributes = "cTiter_pairwise",
        immunity_attributes = "cTiter_x",
        decay_factors = "14.0",
        years_to_wane = config.get("fitness_model", {}).get("max_years_for_distances"),
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/titer_tree_cross_immunities_{build_name}_{segment}.txt"
    log:
        "logs/titer_tree_cross_immunities_{build_name}_{segment}.txt"
    resources:
        mem_mb=8000
    shell:
        """
        python3 flu-forecasting/src/cross_immunity.py \
            --frequencies {input.frequencies} \
            --distances {input.distances} \
            --date-annotations {input.date_annotations} \
            --distance-attributes {params.distance_attributes} \
            --immunity-attributes {params.immunity_attributes} \
            --decay-factors {params.decay_factors} \
            --output {output} 2>&1 | tee {log}
        """

#
# Configure distance maps (for amino acid and other distances).
#

def _get_build_distance_map_config(wildcards):
    distance_config = distance_map_config[
        (distance_map_config["lineage"] == config["builds"][wildcards.build_name].get("lineage")) &
        (distance_map_config["segment"] == wildcards.segment)
    ]

    if distance_config.shape[0] > 0:
        return distance_config
    else:
        return None

def _get_distance_comparisons_by_lineage_and_segment(wildcards):
    config = _get_build_distance_map_config(wildcards)
    return " ".join(config.loc[:, "compare_to"].values)

def _get_distance_attributes_by_lineage_and_segment(wildcards):
    config = _get_build_distance_map_config(wildcards)
    return " ".join(config.loc[:, "attribute"].values)

def _get_distance_maps_by_lineage_and_segment(wildcards):
    distance_config = _get_build_distance_map_config(wildcards)
    lineage = config["builds"][wildcards.build_name].get("lineage")
    return [
        "config/distance_maps/{lineage}/{wildcards.segment}/{distance_map}.json".format(wildcards=wildcards, lineage=lineage, distance_map=distance_map)
        for distance_map in distance_config.loc[:, "distance_map"].values
    ]

rule distances:
    input:
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done",
        alignments = lambda w: [f"{build_dir}/{w.build_name}/{w.segment}/nextalign/masked.gene.{gene}_withInternalNodes.fasta" for gene in GENES[w.segment]],
        distance_maps = _get_distance_maps_by_lineage_and_segment,
    output:
        distances = "builds/{build_name}/{segment}/distances.json",
    params:
        genes = lambda w: GENES[w.segment],
        comparisons = _get_distance_comparisons_by_lineage_and_segment,
        attribute_names = _get_distance_attributes_by_lineage_and_segment,
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/distances_{build_name}_{segment}.txt"
    log:
        "logs/distances_{build_name}_{segment}.txt"
    resources:
        mem_mb=8000,
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --output {output.distances} 2>&1 | tee {log}
        """

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

#
# Define node data table functions.
#

def float_to_datestring(time):
    """Convert a floating point date from TreeTime `numeric_date` to a date string
    """
    # Extract the year and remainder from the floating point date.
    year = int(time)
    remainder = time - year

    # Calculate the day of the year (out of 365 + 0.25 for leap years).
    tm_yday = int(remainder * 365.25)
    if tm_yday == 0:
        tm_yday = 1

    # Construct a date object from the year and day of the year.
    date = datetime.datetime.strptime("%s-%s" % (year, tm_yday), "%Y-%j")

    # Build the date string with zero-padded months and days.
    date_string = "%s-%.2i-%.2i" % (date.year, date.month, date.day)

    return date_string

def _get_end_timepoint(wildcards):
    return float_to_datestring(numeric_date(datetime.date.today()))

def _get_annotations_for_node_data(wildcards):
    annotations = ["%s=%s" % (key, value) for key, value in wildcards.items()]
    annotations.append("timepoint=%s" % _get_end_timepoint(wildcards))
    return " ".join(annotations)

def _get_excluded_fields_arg(wildcards):
    if config.get("fitness_model", {}).get("excluded_node_data_fields"):
        return "--excluded-fields %s" % " ".join(config["fitness_model"]["excluded_node_data_fields"])
    else:
        return ""

rule convert_node_data_to_table:
    input:
        tree = rules.refine.output.tree,
        metadata = "builds/{build_name}/metadata.tsv",
        node_data = _get_node_data_for_predictors,
    output:
        table = "builds/{build_name}/{segment}/node-data-table.tsv",
    params:
        annotations = _get_annotations_for_node_data,
        excluded_fields_arg = _get_excluded_fields_arg,
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/convert_node_data_to_table_{build_name}_{segment}.txt"
    log:
        "logs/convert_node_data_to_table_{build_name}_{segment}.txt"
    shell:
        """
        python3 flu-forecasting/scripts/node_data_to_table.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --jsons {input.node_data} \
            --output {output} \
            --annotations {params.annotations} \
            {params.excluded_fields_arg} 2>&1 | tee {log}
        """

rule convert_frequencies_to_table:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.tip_frequencies.output.tip_freq,
    output:
        table = "builds/{build_name}/{segment}/frequencies.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/convert_frequencies_to_table_{build_name}_{segment}.txt"
    log:
        "logs/convert_frequencies_to_table_{build_name}_{segment}.txt"
    shell:
        """
        python3 scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --output {output} 2>&1 | tee {log}
        """

rule merge_node_data_and_frequencies:
    input:
        node_data = rules.convert_node_data_to_table.output.table,
        frequencies = rules.convert_frequencies_to_table.output.table,
    output:
        table = "builds/{build_name}/{segment}/tip-attributes.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/merge_node_data_and_frequencies_{build_name}_{segment}.txt",
    log:
        "logs/merge_node_data_and_frequencies_{build_name}_{segment}.txt"
    params:
        how="inner",
        on=["strain", "is_terminal"],
    shell:
        """
        python3 scripts/join_tables.py \
            --left {input.node_data} \
            --right {input.frequencies} \
            --how {params.how} \
            --on {params.on} \
            --output {output.table} 2>&1 | tee {log}
        """

def _get_delta_months_to_forecast(wildcards):
    return " ".join([str(month) for month in config["fitness_model"]["delta_months"]])

rule target_distances:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table,
    output:
        distances = "builds/{build_name}/{segment}/target-distances.tsv",
    params:
        delta_months = _get_delta_months_to_forecast,
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/target_distances_{build_name}_{segment}.txt"
    log:
        "logs/target_distances_{build_name}_{segment}.txt"
    shell:
        """
        python3 flu-forecasting/scripts/calculate_target_distances.py \
            --tip-attributes {input.attributes} \
            --delta-months {params.delta_months} \
            --output {output} 2>&1 | tee {log}
        """

rule forecast_tips:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table,
        distances = rules.target_distances.output.distances,
        frequencies = rules.tip_frequencies.output.tip_freq,
        model = "models/{model}.json",
    output:
        node_data = "builds/{build_name}/{segment}/forecast_{model}.json",
        frequencies = "auspice/{build_name}_{segment}_{model}_forecast-tip-frequencies.json",
        table = "builds/{build_name}/{segment}/forecast_{model}.tsv",
    params:
        delta_months = _get_delta_months_to_forecast,
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/forecast_tips_{build_name}_{segment}_{model}.txt"
    log:
        "logs/forecast_tips_{build_name}_{segment}_{model}.txt"
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
            --output-table {output.table} 2>&1 | tee {log}
        """

rule calculate_weighted_distance_to_future:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table,
        forecasts = rules.forecast_tips.output.table,
    output:
        node_data = "builds/{build_name}/{segment}/weighted_distances_{model}.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/calculate_weighted_distance_to_future_{build_name}_{segment}_{model}.txt"
    log:
        "logs/calculate_weighted_distance_to_future_{build_name}_{segment}_{model}.txt"
    shell:
        """
        python3 flu-forecasting/scripts/calculate_weighted_distances.py \
            --tip-attributes {input.attributes} \
            --forecasts {input.forecasts} \
            --distance-attribute-name weighted_distance_to_future_by_{wildcards.model} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule merge_weighted_distances_to_future:
    input:
        distances = expand("builds/{{build_name}}/{{segment}}/weighted_distances_{model}.json", model=config.get("fitness_model", {}).get("models", [])),
    output:
        node_data = "builds/{build_name}/{segment}/weighted_distances.json",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/merge_weighted_distances_to_future_{build_name}_{segment}.txt"
    log:
        "logs/merge_weighted_distances_to_future_{build_name}_{segment}.txt"
    shell:
        """
        python3 flu-forecasting/scripts/merge_weighted_distances_to_future.py \
            --distances {input.distances:q} \
            --output {output.node_data} 2>&1 | tee {log}
        """
