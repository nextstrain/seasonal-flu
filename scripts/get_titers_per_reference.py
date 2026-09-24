"""Split titers into separate files per reference virus.
"""
import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--titers", required=True, help="TSV of titers to split by reference")
    parser.add_argument("--output-references", required=True, help="text file listing the references with titer outputs")
    parser.add_argument("--output-titers-directory", required=True, help="directory where split titers TSV are placed per reference")

    args = parser.parse_args()

    titers = pd.read_csv(
        args.titers,
        sep="\t",
    )

    # Find references with autologous and heterologous measurements.
    distinct_pairs = titers.loc[:, ["virus_strain", "serum_strain"]].drop_duplicates()
    print(f"Found {distinct_pairs.shape[0]} distinct pairs")

    has_autologous_measurement = (distinct_pairs["virus_strain"] == distinct_pairs["serum_strain"])
    autologous_references = set(distinct_pairs.loc[has_autologous_measurement, "serum_strain"].drop_duplicates().values)
    print(f"Found {len(autologous_references)} autologous references")

    has_heterologous_measurement = (distinct_pairs["virus_strain"] != distinct_pairs["serum_strain"])
    heterologous_references = set(distinct_pairs.loc[has_heterologous_measurement, "serum_strain"].drop_duplicates().values)
    print(f"Found {len(heterologous_references)} heterologous references")

    selected_references = autologous_references & heterologous_references
    print(f"Found {len(selected_references)} references")

    selected_titers = titers[titers["serum_strain"].isin(selected_references)].copy()
    selected_titers["reference_path"] = selected_titers["serum_strain"].apply(
        lambda strain: strain.replace("/", "_")
    )
    selected_reference_paths = selected_titers["reference_path"].drop_duplicates().values

    for reference, reference_titers in selected_titers.groupby("reference_path"):
        reference_titers.to_csv(
            f"{args.output_titers_directory}/{reference}.tsv",
            sep="\t",
            index=False,
        )

    with open(args.output_references, "w", encoding="utf-8") as oh:
        for reference in selected_reference_paths:
            print(reference, file=oh)
