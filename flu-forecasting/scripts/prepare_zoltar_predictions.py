#!/usr/bin/env python3
import argparse
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--forecasts", nargs="+", required=True, help="TSV of forecasts per strain")
    parser.add_argument("--output", required=True, help="TSV of forecasts per clade")

    args = parser.parse_args()

    # Load forecasts.
    forecasts = pd.concat(
        [
            pd.read_csv(
                forecast_file,
                sep="\t",
            )
            for forecast_file in args.forecasts
        ],
        ignore_index=True)

    # Prepare empty predictions.
    predictions = {
        "meta": {},
        "predictions": [],
    }

    # Build predictions per clade.
    for clade, clade_forecasts in forecasts.groupby("clade"):
        predictions["predictions"].append({
            "unit": f"global:{clade}",
            "target": "frequency in one year",
            "class": "sample",
            "prediction": {
                "sample": clade_forecasts["projected_frequency"].values.tolist(),
            }
        })

    # Save forecasts.
    with open(args.output, "w", encoding="utf8") as oh:
        json.dump(predictions, oh)
