import argparse
import json
import csv

def write_outfile(output_site, grouped_site):
    with open(output_site, "w") as outfile:
            first_record = list(grouped_site.values())[0]
            fieldnames = list(first_record.keys())
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()

            for record in grouped_site.values():
                writer.writerow(record)


def parse_json(input_file, output_ga, output_rf, output_freq, output_raw, output_freq_forecast, model_version):
    """
    Extracts empirical frequency and inferred frequency and fitness values from <model>_results.json.
    Parses "freq" and "ga" (growth advantage) from MLR model results.
    Parses "freq" and "delta" (relative fitness) and "ga" (growth advantage) from Latent model results.
    Parses "raw_freq" from MLR or Latent model as raw_freq.
    """

    # Read in JSON file
    with open(input_file, "r") as file:
        data = json.load(file)

        # Groups by site
        grouped_ga = {}
        grouped_rf = {}
        grouped_freq = {}
        grouped_raw_freq = {}

        if output_freq:
            # Parse freq from MLR or Latent model
            for record in data["data"]:
                if record["site"] == "freq" and record["ps"] in ["mean", "median", "HDI_95_upper", "HDI_95_lower"]:
                    key = (record["location"], record["variant"], record["date"])
                    if key not in grouped_freq:
                        grouped_freq[key] = {"location": record["location"], "date": record["date"], "variant": record["variant"]}
                    grouped_freq[key][record["ps"]] = record["value"]

            # Write output <model>/freq.tsv
            print("Parsing freq from model results.")
            write_outfile(output_freq, grouped_freq)

        # Parse forecast freq from MLR model
        if output_freq_forecast:
            print("Parsing forecast freq from model results.")
            grouped_freq_forecast = {}
            for record in data["data"]:
                if record["site"] == "freq_forecast" and record["ps"] in ["mean", "median", "HDI_95_upper", "HDI_95_lower"]:
                    key = (record["location"], record["variant"], record["date"])
                    if key not in grouped_freq_forecast:
                        grouped_freq_forecast[key] = {"location": record["location"], "date": record["date"], "variant": record["variant"]}
                    grouped_freq_forecast[key][record["ps"]] = record["value"]

            # Write output <model>/freq.tsv
            write_outfile(output_freq_forecast, grouped_freq_forecast)

        # Parse raw_freq from model if requested
        if output_raw:
            print("Parsing raw_freq from model results.")
            for record in data["data"]:
                if record["site"] == "raw_freq":
                    key = (record["location"], record["variant"], record["date"])
                    if key not in grouped_raw_freq:
                        grouped_raw_freq[key] = {"location": record["location"], "date": record["date"], "variant": record["variant"]}
                    grouped_raw_freq[key]["raw_freq"] = record["value"]

            # Write output raw_freq.tsv
            write_outfile(output_raw, grouped_raw_freq)

         # Parse ga (MLR model)
        if model_version == "MLR":
            print("Parsing ga (growth advantage) from MLR model results.")
            for record in data["data"]:
                if record["site"] == "ga" and record["ps"] in ["mean", "median", "HDI_95_upper", "HDI_95_lower"]:
                    key = (record["location"], record["variant"])
                    if key not in grouped_ga:
                        grouped_ga[key] = {"location": record["location"], "variant": record["variant"]}
                    grouped_ga[key][record["ps"]] = record["value"]

            # Write output mlr/ga.tsv
            write_outfile(output_ga, grouped_ga)

        # Parse rf and ga (Latent model)
        elif model_version == "Latent":
            print("Parsing delta (relative fitness) from Latent model results.")
            print("Parsing ga (growth advantage) from Latent model results.")
            for record in data["data"]:
                if record["site"] == "delta" and record["ps"] in ["mean", "median", "HDI_95_upper", "HDI_95_lower"]:
                    key = (record["location"], record["variant"], record["date"])
                    if key not in grouped_rf:
                        grouped_rf[key] = {"location": record["location"], "date": record["date"], "variant": record["variant"]}
                    grouped_rf[key][record["ps"]] = record["value"]
                if record["site"] == "ga" and record["ps"] in ["mean", "median", "HDI_95_upper", "HDI_95_lower"]:
                    key = (record["location"], record["variant"], record["date"])
                    if key not in grouped_ga:
                        grouped_ga[key] = {"location": record["location"], "date": record["date"], "variant": record["variant"]}
                    grouped_ga[key][record["ps"]] = record["value"]

            # Write outout latent/rf.tsv
            write_outfile(output_rf, grouped_rf)
            # Write output latent/ga.tsv
            write_outfile(output_ga, grouped_ga)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and parse MLR-model JSON data")
    parser.add_argument("--input", required=True, help="Path to the MLR output JSON file (<model>_results.json)")
    parser.add_argument("--outga", help="Path to filtered and parsed GA (growth advantage) TSV file (mlr/ga.tsv)")
    parser.add_argument("--outrf", help="Path to filtered and parsed RF (relative fitness) TSV file (latent/rf.tsv)")
    parser.add_argument("--outfreq", help="Path to filtered and parsed freq TSV file (<model>/freq.tsv)")
    parser.add_argument("--outraw", help="Path to empirical freq TSV file (raw_freq.tsv)")
    parser.add_argument("--outfreqforecast", help="Path to forecast frequencies TSV file")
    parser.add_argument("--model", default="MLR", help="Model version ['Latent', 'MLR']")
    args = parser.parse_args()
    parse_json(args.input, args.outga, args.outrf, args.outfreq, args.outraw, args.outfreqforecast, args.model)
