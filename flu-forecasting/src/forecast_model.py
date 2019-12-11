"""Forecast given tip data into future using the given previously trained model.
"""
import argparse
import json
import numpy as np
import pandas as pd

from fit_model import DistanceExponentialGrowthModel
from weighted_distances import get_distances_by_sample_names


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints with standardized predictors")
    parser.add_argument("--distances", help="tab-delimited file of distances between pairs of samples")
    parser.add_argument("--frequencies", help="JSON representing historical frequencies to project from")
    parser.add_argument("--model", required=True, help="JSON representing the model fit with training and cross-validation results, beta coefficients for predictors, and summary statistics")
    parser.add_argument("--delta-months", required=True, type=int, nargs="+", help="number of months to project clade frequencies into the future")
    parser.add_argument("--output-node-data", required=True, help="node data JSON of forecasts for the given tips")
    parser.add_argument("--output-frequencies", help="frequencies JSON extended with forecasts for the given tips")
    parser.add_argument("--output-table", help="table of forecasts for the given tips")

    args = parser.parse_args()

    # Load standardized tip attributes subsetting to tip name, clade, frequency,
    # and requested predictors.
    tips = pd.read_csv(
        args.tip_attributes,
        sep="\t",
        parse_dates=["timepoint"]
    )

    # Load distances.
    distances = pd.read_csv(args.distances, sep="\t")
    distances_by_sample_names = get_distances_by_sample_names(distances)

    # Load model details
    with open(args.model, "r") as fh:
        model_json = json.load(fh)

    predictors = model_json["predictors"]
    cost_function = model_json["cost_function"]
    l1_lambda = model_json["l1_lambda"]
    coefficients = np.array(model_json["coefficients_mean"])
    mean_stds = np.array(model_json["mean_stds_mean"])

    delta_month = args.delta_months[-1]
    delta_time = delta_month / 12.0
    delta_offset = pd.DateOffset(months=delta_month)

    model = DistanceExponentialGrowthModel(
        predictors=predictors,
        delta_time=delta_time,
        cost_function=cost_function,
        l1_lambda=l1_lambda,
        distances=distances_by_sample_names
    )
    model.coef_ = coefficients
    model.mean_stds_ = mean_stds

    # collect fitness and projection
    forecasts_df = model.predict(tips)
    forecasts_df["weighted_distance_to_future_by_%s" % "-".join(predictors)] = forecasts_df["y"]

    # collect dicts from dataframe
    strain_to_fitness = {}
    strain_to_weighted_distance_to_future = {}
    for index, row in forecasts_df.iterrows():
        strain_to_fitness[row['strain']] = row['fitness']
        strain_to_weighted_distance_to_future[row['strain']] = row['y']

    # populate node data
    node_data = {}
    strains = list(tips['strain'])
    for strain in strains:
        node_data[strain] = {
            "fitness": strain_to_fitness[strain],
            "weighted_distance_to_future": strain_to_weighted_distance_to_future[strain]
        }

    # output to file
    with open(args.output_node_data, "w") as jsonfile:
        json.dump({"nodes": node_data}, jsonfile, indent=1)

    # load historic frequencies
    with open(args.frequencies, "r") as fh:
        frequencies = json.load(fh)
    frequencies.pop("generated_by")

    pivots = frequencies.pop("pivots")
    projection_pivot = pivots[-1]

    forecasts = []
    for delta_month in args.delta_months:
        delta_time = delta_month / 12.0
        delta_offset = pd.DateOffset(months=delta_month)

        model = DistanceExponentialGrowthModel(
            predictors=predictors,
            delta_time=delta_time,
            cost_function=cost_function,
            l1_lambda=l1_lambda,
            distances=distances_by_sample_names
        )
        model.coef_ = coefficients
        model.mean_stds_ = mean_stds

        # collect fitness and projection
        forecasts_df = model.predict(tips)
        forecasts_df["future_timepoint"] = forecasts_df["timepoint"] + delta_offset

        # collect dicts from dataframe
        strain_to_projected_frequency = {}
        for index, row in forecasts_df.iterrows():
            strain_to_projected_frequency[row['strain']] = row['projected_frequency']

        # extend frequencies
        for strain in frequencies.keys():
            trajectory = frequencies[strain]['frequencies']
            if strain in strain_to_projected_frequency:
                trajectory.append(strain_to_projected_frequency[strain])
            else:
                trajectory.append(0.0)

        # extend pivots
        pivots.append(projection_pivot + delta_time)

        # Collect forecast data frames, if requested.
        if args.output_table:
            forecasts.append(forecasts_df)

    # reconnect pivots and label projection pivot
    frequencies['pivots'] = pivots
    frequencies['projection_pivot'] = projection_pivot

    # output to file
    with open(args.output_frequencies, "w") as jsonfile:
        json.dump(frequencies, jsonfile, indent=1)

    # Save forecasts table, if requested.
    if args.output_table:
        all_forecasts = pd.concat(forecasts, ignore_index=True)
        all_forecasts.to_csv(args.output_table, sep="\t", index=False, header=True, na_rep="N/A")
