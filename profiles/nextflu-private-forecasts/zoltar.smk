rule calculate_clade_frequency_forecasts:
    input:
        forecast="builds/{build_name}/ha/forecast_{model}.tsv",
        clades="builds/{build_name}/ha/clades.json",
    output:
        forecast="builds/{build_name}/clade_forecast_{model}.tsv",
    conda: "../../workflow/envs/nextstrain.yaml"
    benchmark:
        "benchmarks/calculate_clade_frequency_forecasts_{build_name}_{model}.txt"
    log:
        "logs/calculate_clade_frequency_forecasts_{build_name}_{model}.txt"
    shell:
        """
        python3 flu-forecasting/scripts/calculate_clade_frequency_forecasts.py \
            --forecasts {input.forecast} \
            --clades {input.clades} \
            --output {output.forecast} 2>&1 | tee {log}
        """

rule prepare_zoltar_predictions:
    input:
        forecasts=expand("builds/{build_name}/clade_forecast_{{model}}.tsv", build_name=list(config["builds"].keys()))
    output:
        forecasts="forecasts/predictions_{model}.json",
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
        expand("builds/{build_name}/clade_forecast_{model}.tsv", build_name=list(config["builds"].keys()), model=config["fitness_model"]["models"])
