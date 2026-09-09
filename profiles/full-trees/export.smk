ruleorder: export_full_trees > export
ruleorder: export_full_trees > export_private

human_antigenic_advance_attribute = "kikawa_2026_cTiterSub"
ferret_antigenic_advance_attribute = "cell_hi_cTiterSub"

rule calculate_antigenic_advance_from_human_titers:
    input:
        substitutions=lambda wildcards: f"profiles/full-trees/{config['builds'][wildcards.build_name]['lineage']}/{wildcards.titer_collection}.parquet",
        tree="builds/{build_name}/ha/tree.nwk",
        translations_done="builds/{build_name}/{segment}/translations.done",
    params:
        genes="HA1",
        translations="builds/{build_name}/{segment}/translations/HA1_withInternalNodes.fasta",
    output:
        antigenic_advance="builds/{build_name}/{segment}/titers-sub-human/{titer_collection}.json",
    log:
        "logs/calculate_antigenic_advance_from_human_titers_{build_name}_{segment}_{titer_collection}.txt",
    shell:
        r"""
        python scripts/calculate_antigenic_advance.py \
            --substitution-weights {input.substitutions} \
            --tree {input.tree} \
            --alignment {params.translations:q} \
            --gene-names {params.genes:q} \
            --attribute-prefix "{wildcards.titer_collection}_" \
            --output-node-data {output.antigenic_advance} 2>&1 | tee {log}
        """

rule download_mlr_json:
    output:
        mlr="builds/{build_name}/ha/mlr/model.json",
    params:
        model_url=lambda wildcards: f"https://data.nextstrain.org/files/workflows/forecasts-flu/gisaid/emerging_haplotype/{config['builds'][wildcards.build_name]['lineage']}/region/mlr/MLR_results.json",
    shell:
        r"""
        curl --compressed -o {output.mlr:q} {params.model_url:q}
        """

rule parse_frequencies_and_ga_from_mlr_json:
    input:
        mlr="builds/{build_name}/ha/mlr/model.json",
    output:
        forecasts="builds/{build_name}/ha/mlr/freq_forecast.tsv",
        fitnesses="builds/{build_name}/ha/mlr/fitnesses.tsv",
    shell:
        r"""
        python scripts/parse-json.py \
            --input {input.mlr} \
            --outfreqforecast {output.forecasts} \
            --outga {output.fitnesses}
        """

rule calculate_human_antigenic_distance_to_the_future:
    input:
        titer_model="builds/{build_name}/ha/titers-sub-human/kikawa_2026.json",
        titers=lambda wildcards: f"data/{config['builds'][wildcards.build_name]['lineage']}/who_ferret_cell_hi_titers.tsv",
        forecasts="builds/{build_name}/ha/mlr/freq_forecast.tsv",
        tip_attributes="builds/{build_name}/metadata.tsv",
        ha1_sequences_dir="builds/{build_name}/ha/translations",
    output:
        node_data="builds/{build_name}/human_antigenic_distance_to_future.json",
        table="builds/{build_name}/human_antigenic_distance_to_future.tsv",
    params:
        min_date="2026-02-01",
        min_reference_year=2022,
        attribute_name="antigenic_distance_to_future_human",
        ha1_sequences="builds/{build_name}/ha/translations/HA1.fasta",
    shell:
        r"""
        python scripts/calculate_titer_distance_from_candidates.py \
            --titer-model {input.titer_model} \
            --titers {input.titers} \
            --tip-attributes {input.tip_attributes} \
            --forecasts {input.forecasts} \
            --ha1-sequences {params.ha1_sequences} \
            --min-date {params.min_date} \
            --min-reference-year {params.min_reference_year} \
            --attribute-name {params.attribute_name} \
            --output-node-data {output.node_data} \
            --output-table {output.table}
        """

rule calculate_ferret_antigenic_distance_to_the_future:
    input:
        titer_model="builds/{build_name}/ha/titers-sub-model/cell_hi.json",
        titers=lambda wildcards: f"data/{config['builds'][wildcards.build_name]['lineage']}/who_ferret_cell_hi_titers.tsv",
        forecasts="builds/{build_name}/ha/mlr/freq_forecast.tsv",
        tip_attributes="builds/{build_name}/metadata.tsv",
        ha1_sequences_dir="builds/{build_name}/ha/translations",
    output:
        node_data="builds/{build_name}/ferret_antigenic_distance_to_future.json",
        table="builds/{build_name}/ferret_antigenic_distance_to_future.tsv",
    params:
        min_date="2026-02-01",
        min_reference_year=2022,
        attribute_name="antigenic_distance_to_future_ferret",
        ha1_sequences="builds/{build_name}/ha/translations/HA1.fasta",
    shell:
        r"""
        python scripts/calculate_titer_distance_from_candidates.py \
            --titer-model {input.titer_model} \
            --titers {input.titers} \
            --tip-attributes {input.tip_attributes} \
            --forecasts {input.forecasts} \
            --ha1-sequences {params.ha1_sequences} \
            --min-date {params.min_date} \
            --min-reference-year {params.min_reference_year} \
            --attribute-name {params.attribute_name} \
            --output-node-data {output.node_data} \
            --output-table {output.table}
        """

def get_antigenic_advance_from_human_titers(wildcards):
    if "vic" not in wildcards.build_name and wildcards.segment == "ha":
        return [
            f"builds/{wildcards.build_name}/{wildcards.segment}/titers-sub-human/{titer_collection}.json"
            for titer_collection in ["kikawa_2026"]
        ]
    else:
        return []

def get_antigenic_distances_to_future(wildcards):
    if "vic" not in wildcards.build_name and wildcards.segment == "ha":
        return [
            f"{build_dir}/{wildcards.build_name}/human_antigenic_distance_to_future.json",
            f"{build_dir}/{wildcards.build_name}/ferret_antigenic_distance_to_future.json",
        ]
    else:
        return []

rule export_full_trees:
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
        node_data = _get_node_data_by_wildcards,
        antigenic_advance_from_human_titers=get_antigenic_advance_from_human_titers,
        antigenic_distance_to_future=get_antigenic_distances_to_future,
        auspice_config = lambda w: config['builds'][w.build_name]['auspice_config'],
        description = lambda w: config['builds'][w.build_name].get("description", "config/description.md"),
        lat_longs = config.get('lat-longs', "config/lat_longs.tsv"),
    output:
        auspice_json = "auspice/{build_name}_{segment}.json"
    benchmark:
        "benchmarks/export_{build_name}_{segment}.txt"
    log:
        "logs/export_{build_name}_{segment}.txt"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} {input.antigenic_advance_from_human_titers} {input.antigenic_distance_to_future} \
            --include-root-sequence-inline \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule get_antigenic_advance_attributes_per_tip:
    input:
        auspice_json="auspice/{build_name}_{segment}.json",
    output:
        tip_attributes="builds/{build_name}/{segment}/tip_attributes.tsv",
    params:
        attributes=[
            "emerging_haplotype_ha",
            human_antigenic_advance_attribute,
            ferret_antigenic_advance_attribute,
        ],
    shell:
        r"""
        python scripts/auspice_tree_to_table.py \
            --tree {input.auspice_json} \
            --output-metadata {output.tip_attributes} \
            --attributes {params.attributes}
        """

rule plot_fitness_by_antigenic_advance:
    input:
        tip_attributes="builds/{build_name}/{segment}/tip_attributes.tsv",
        fitnesses="builds/{build_name}/{segment}/mlr/fitnesses.tsv",
        auspice_config = lambda w: config['builds'][w.build_name]['auspice_config'],
    output:
        figure="figures/{build_name}_{segment}_fitness_by_antigenic_advance.png",
    params:
        variant_attribute="emerging_haplotype_ha",
        human_antigenic_advance_attribute=human_antigenic_advance_attribute,
        ferret_antigenic_advance_attribute=ferret_antigenic_advance_attribute,
    shell:
        r"""
        python scripts/plot_fitness_by_antigenic_advance.py \
            --tip-attributes {input.tip_attributes} \
            --fitnesses {input.fitnesses} \
            --auspice-config {input.auspice_config} \
            --variant-attribute {params.variant_attribute} \
            --human-antigenic-advance-attribute {params.human_antigenic_advance_attribute} \
            --ferret-antigenic-advance-attribute {params.ferret_antigenic_advance_attribute} \
            --output-figure {output.figure}
        """
