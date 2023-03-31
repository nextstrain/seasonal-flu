#!/usr/bin/env python3
import argparse
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--forecasts", required=True, help="TSV of forecasts per strain")
    parser.add_argument("--output", required=True, help="TSV of forecasts per clade")

    args = parser.parse_args()

    # Load forecasts.
    forecasts = pd.read_csv(
        args.forecasts,
        sep="\t",
    )

    # Prepare empty predictions.
    predictions = {
        "meta": {},
        "predictions": [],
    }

    # Get clade frequencies at last observed timepoint.
    last_observed_pivot = forecasts.query("observed == True")["pivot"].max()
    current_clade_frequencies = forecasts.query(f"pivot == {last_observed_pivot}").groupby(["model", "sample", "clade"])["frequency"].sum().reset_index()

    # Get number of samples.
    total_samples = forecasts["sample"].drop_duplicates().size

    # Select all active clades. These are detectable but not fixed in all
    # samples.
    active_clades = current_clade_frequencies.query(
        "(frequency >= 0.01) & (frequency <= 0.99)"
    ).groupby("clade")["sample"].count().reset_index().query(
        f"sample == {total_samples}"
    )["clade"].values

    # Get predicted clade frequencies for active clades at the last timepoint.
    last_pivot = forecasts["pivot"].max()
    is_last_pivot = (forecasts["pivot"] == last_pivot)
    is_active_clade = forecasts["clade"].isin(active_clades)
    future_clade_frequencies = forecasts.loc[
        (is_last_pivot) & (is_active_clade)
    ].groupby([
        "model",
        "sample",
        "clade",
    ])["frequency"].sum().reset_index()

    # Build predictions per clade.
    for clade, clade_forecasts in future_clade_frequencies.groupby("clade"):
        predictions["predictions"].append({
            "unit": f"global:{clade}",
            "target": "frequency in one year",
            "class": "sample",
            "prediction": {
                "sample": clade_forecasts["frequency"].values.tolist(),
            }
        })

    # Save forecasts.
    with open(args.output, "w", encoding="utf8") as oh:
        json.dump(predictions, oh)
