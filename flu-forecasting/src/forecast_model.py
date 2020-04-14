"""Forecast given tip data into future using the given previously trained model.
"""
import argparse
import csv
import json
import numpy as np
import pandas as pd
import sys

from fit_model import DistanceExponentialGrowthModel
from weighted_distances import get_distances_by_sample_names


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints with standardized predictors")
    parser.add_argument("--distances", help="tab-delimited file of distances between pairs of samples")
    parser.add_argument("--frequencies", help="JSON representing historical frequencies to project from")
    parser.add_argument("--model", required=True, help="JSON representing the model fit with training and cross-validation results, beta coefficients for predictors, and summary statistics")
    parser.add_argument("--delta-months", required=True, type=int, nargs="+", help="number of months to project clade frequencies into the future")
    parser.add_argument("--output-node-data", help="node data JSON of forecasts for the given tips")
    parser.add_argument("--output-frequencies", help="frequencies JSON extended with forecasts for the given tips")
    parser.add_argument("--output-table", help="table of forecasts for the given tips")

    args = parser.parse_args()

    # Confirm that at least one output file has been specified.
    outputs = [
        args.output_node_data,
        args.output_frequencies,
        args.output_table
    ]
    outputs_missing =[output is None for output in outputs]
    if all(outputs_missing):
        print("ERROR: No output files were specified", file=sys.stderr)
        sys.exit(1)

    # Load standardized tip attributes subsetting to tip name, clade, frequency,
    # and requested predictors.
    tips = pd.read_csv(
        args.tip_attributes,
        sep="\t",
        parse_dates=["timepoint"]
    )
    tips = tips[tips["frequency"] > 0].copy()

    # Load distances.
    with open(args.distances, "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        distances_by_sample_names = get_distances_by_sample_names(reader)

    # Load model details
    with open(args.model, "r") as fh:
        model_json = json.load(fh)

    predictors = model_json["predictors"]

    # Confirm that the requested predictors have columns in the given tip
    # attributes table.
    missing_predictors = [
        predictor
        for predictor in predictors
        if predictor not in tips.columns
    ]
    if len(missing_predictors) > 0:
        print(
            "ERROR: the following predictors are missing columns in the tip attributes table:", ", ".join(missing_predictors),
            file=sys.stderr
        )
        sys.exit(1)

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
    final_forecasts = model.predict(tips)
    final_forecasts["weighted_distance_to_future_by_%s" % "-".join(predictors)] = final_forecasts["y"]
    final_forecasts["future_timepoint"] = final_forecasts["timepoint"] + delta_offset

    # collect dicts from dataframe
    strain_to_fitness = {}
    strain_to_future_timepoint = {}
    strain_to_projected_frequency = {}
    strain_to_weighted_distance_to_future = {}
    for index, row in final_forecasts.iterrows():
        strain_to_fitness[row['strain']] = row['fitness']
        strain_to_future_timepoint[row['strain']] = row["future_timepoint"].strftime("%Y-%m-%d")
        strain_to_projected_frequency[row['strain']] = row['projected_frequency']
        strain_to_weighted_distance_to_future[row['strain']] = row['y']

    # output to file
    if args.output_node_data:
        # populate node data
        node_data = {}
        strains = list(tips['strain'])
        for strain in strains:
            node_data[strain] = {
                "fitness": strain_to_fitness[strain],
                "future_timepoint": strain_to_future_timepoint[strain],
                "projected_frequency": strain_to_projected_frequency[strain]
            }

        with open(args.output_node_data, "w") as jsonfile:
            json.dump({"nodes": node_data}, jsonfile, indent=1)

    # load historic frequencies
    if args.frequencies:
        with open(args.frequencies, "r") as fh:
            frequencies = json.load(fh)

        pivots = frequencies.pop("pivots")
        generated_by = frequencies.pop("generated_by", None)
        projection_pivot = pivots[-1]
    else:
        frequencies = None

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

        if frequencies is not None:
            # extend frequencies
            for strain in frequencies.keys():
                trajectory = frequencies[strain]['frequencies']
                if strain in strain_to_projected_frequency:
                    trajectory.append(strain_to_projected_frequency[strain])
                else:
                    trajectory.append(0.0)

            # extend pivots
            pivots.append(projection_pivot + delta_time)

    # reconnect pivots and label projection pivot
    if frequencies is not None:
        frequencies['pivots'] = pivots
        frequencies['projection_pivot'] = projection_pivot

    # output to file
    if args.output_frequencies:
        with open(args.output_frequencies, "w") as jsonfile:
            json.dump(frequencies, jsonfile, indent=1)

    # Save forecasts table, if requested.
    if args.output_table:
        # Annotate the model used to make this forecast for self-documentation
        # and debugging.
        final_forecasts["model"] = "-".join(predictors)

        # Annotate amino acid sequences along with forecasts to enable weighted
        # distance calculations to the estimated future population.
        tip_sequences = tips.loc[:, ["strain", "aa_sequence"]].copy()
        final_forecasts = final_forecasts.merge(
            tip_sequences,
            on="strain"
        )

        final_forecasts.to_csv(args.output_table, sep="\t", index=False, header=True, na_rep="N/A")
