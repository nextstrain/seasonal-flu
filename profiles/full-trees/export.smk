ruleorder: export_full_trees > export

rule calculate_human_antigenic_distance_to_the_future:
    input:
        titer_model="builds/{build_name}/ha/titers-sub-model/kikawa_2025b.json",
        titers=lambda wildcards: f"data/{config['builds'][wildcards.build_name]['lineage']}/who_ferret_cell_hi_titers.tsv",
        forecasts=lambda wildcards: f"../forecasts-flu/results/gisaid/emerging_haplotype/{config['builds'][wildcards.build_name]['lineage']}/region/mlr/freq_forecast.tsv",
        tip_attributes="builds/{build_name}/ha/emerging_haplotypes.tsv",
        ha1_sequences_dir="builds/{build_name}/ha/translations",
    output:
        node_data="builds/{build_name}/human_antigenic_distance_to_future.json",
        table="builds/{build_name}/human_antigenic_distance_to_future.tsv",
    params:
        min_date="2025-05-01",
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
        forecasts=lambda wildcards: f"../forecasts-flu/results/gisaid/emerging_haplotype/{config['builds'][wildcards.build_name]['lineage']}/region/mlr/freq_forecast.tsv",
        tip_attributes="builds/{build_name}/ha/emerging_haplotypes.tsv",
        ha1_sequences_dir="builds/{build_name}/ha/translations",
    output:
        node_data="builds/{build_name}/ferret_antigenic_distance_to_future.json",
        table="builds/{build_name}/ferret_antigenic_distance_to_future.tsv",
    params:
        min_date="2025-05-01",
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

def get_antigenic_distances_to_future(wildcards):
    if "vic" in wildcards.build_name:
        return []
    else:
        return [
            f"{build_dir}/{wildcards.build_name}/human_antigenic_distance_to_future.json",
            f"{build_dir}/{wildcards.build_name}/ferret_antigenic_distance_to_future.json",
        ]

rule export_full_trees:
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
        node_data = _get_node_data_by_wildcards,
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
            --node-data {input.node_data} {input.antigenic_distance_to_future} \
            --include-root-sequence-inline \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
