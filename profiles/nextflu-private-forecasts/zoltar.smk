rule calculate_clade_frequency_forecasts:
    input:
        tree="auspice/{build_name}_ha.json",
        forecast_frequencies="auspice/{build_name}_ha_{model}_forecast-tip-frequencies.json",
    output:
        forecast="forecasts/by_build_and_model/{build_name}/{model}.tsv",
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/calculate_clade_frequency_forecasts_{build_name}_{model}.txt"
    log:
        "logs/calculate_clade_frequency_forecasts_{build_name}_{model}.txt"
    shell:
        """
        python3 scripts/forecast_frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.forecast_frequencies} \
            --annotations sample={wildcards.build_name} model={wildcards.model} \
            --output {output.forecast} 2>&1 | tee {log}
        """

rule aggregate_forecasts_by_model:
    input:
        forecasts=expand("forecasts/by_build_and_model/{build_name}/{{model}}.tsv", build_name=list(config["builds"].keys()))
    output:
        forecasts="forecasts/by_model/{model}.tsv",
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/aggregate_forecasts_by_model_{model}.txt"
    log:
        "logs/aggregate_forecasts_by_model_{model}.txt"
    shell:
        """
        tsv-append -H {input.forecasts} > {output.forecasts} 2> {log}
        """

rule aggregate_forecasts:
    input:
        forecasts=expand("forecasts/by_model/{model}.tsv", model=config["fitness_model"]["models"]),
    output:
        forecasts="forecasts/all.tsv",
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/aggregate_forecasts.txt"
    log:
        "logs/aggregate_forecasts.txt"
    shell:
        """
        tsv-append -H {input.forecasts} > {output.forecasts} 2> {log}
        """

rule prepare_zoltar_predictions:
    input:
        forecasts="forecasts/by_model/{model}.tsv",
    output:
        forecasts="forecasts/zoltar/{model}.json",
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/prepare_zoltar_predictions_{model}.txt"
    log:
        "logs/prepare_zoltar_predictions_{model}.txt"
    shell:
        """
        python3 flu-forecasting/scripts/prepare_zoltar_predictions.py \
            --forecasts {input.forecasts:q} \
            --output {output.forecasts} 2>&1 | tee {log}
        """

rule all_forecasts:
    input:
        "forecasts/all.tsv",
        expand("forecasts/zoltar/{model}.json", model=config["fitness_model"]["models"]),
