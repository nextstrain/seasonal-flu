import argparse
import matplotlib.pyplot as plt
import json
import numpy as np
import pandas as pd
import seaborn as sns


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tip-attributes", required=True, help="TSV of tip attributes from Auspice JSON including an antigenic advance attribute")
    parser.add_argument("--fitnesses", required=True, help="TSV of MLR-based fitnesses with pivot record included")
    parser.add_argument("--auspice-config", required=True, help="Auspice config JSON with colorings per variant")
    parser.add_argument("--variant-attribute", required=True, help="column in the tip attributes table for the variants to summarize advance and fitness by")
    parser.add_argument("--human-antigenic-advance-attribute", required=True, help="column in the tip attributes table with the antigenic advance from human serology")
    parser.add_argument("--ferret-antigenic-advance-attribute", required=True, help="column in the tip attributes table with the antigenic advance from ferret serology")
    parser.add_argument("--output-figure", required=True, help="output figure")

    args = parser.parse_args()

    variant_attribute = args.variant_attribute
    human_advance_attribute = args.human_antigenic_advance_attribute
    ferret_advance_attribute = args.ferret_antigenic_advance_attribute

    ga_path = args.fitnesses
    antigenic_advance_path = args.tip_attributes
    auspice_config_path = args.auspice_config

    output_path = args.output_figure

    with open(auspice_config_path, "r", encoding="utf-8") as fh:
        colorings = json.load(fh).get("colorings")

    haplotype_colorings = [dict(coloring["scale"]) for coloring in colorings if coloring["key"] == variant_attribute][0]
    haplotype_colorings_domain = list(haplotype_colorings.keys())
    haplotype_colorings_range = list(haplotype_colorings.values())

    ga = pd.read_csv(
        ga_path,
        sep="\t",
    ).query(
        "location == 'hierarchical'"
    ).rename(
        columns={
            "variant": "emerging_haplotype",
            "median": "ga_median",
            "HDI_95_lower": "ga_HDI_95_lower",
            "HDI_95_upper": "ga_HDI_95_upper",
        }
    )

    ga.loc[ga["emerging_haplotype"] == "other", "emerging_haplotype"] = "unassigned"
    ga = ga[ga["emerging_haplotype"] != "unassigned"].copy()

    antigenic_advance = pd.read_csv(
        antigenic_advance_path,
        sep="\t",
    )

    antigenic_advance_per_haplotype = antigenic_advance.groupby(variant_attribute).aggregate(
        human_advance_median=(human_advance_attribute, "median"),
        human_advance_std=(human_advance_attribute, "std"),
        ferret_advance_median=(ferret_advance_attribute, "median"),
        ferret_advance_std=(ferret_advance_attribute, "std"),
    ).reset_index()

    antigenic_advance_per_haplotype["human_advance_lower"] = (
        antigenic_advance_per_haplotype["human_advance_median"] - antigenic_advance_per_haplotype["human_advance_std"]
    )
    antigenic_advance_per_haplotype["human_advance_upper"] = (
        antigenic_advance_per_haplotype["human_advance_median"] + antigenic_advance_per_haplotype["human_advance_std"]
    )

    antigenic_advance_per_haplotype["ferret_advance_lower"] = (
        antigenic_advance_per_haplotype["ferret_advance_median"] - antigenic_advance_per_haplotype["ferret_advance_std"]
    )
    antigenic_advance_per_haplotype["ferret_advance_upper"] = (
        antigenic_advance_per_haplotype["ferret_advance_median"] + antigenic_advance_per_haplotype["ferret_advance_std"]
    )

    ga_by_antigenic_advance = ga.merge(
        antigenic_advance_per_haplotype,
        left_on="emerging_haplotype",
        right_on=variant_attribute,
    ).round(3)

    def plot_panel(ax, df, x_median_col, x_lower_col, x_upper_col, x_title, color_map):
        for haplotype, group in df.groupby("emerging_haplotype"):
            color = color_map.get(haplotype, "gray")

            x_err = np.array([
                group[x_median_col] - group[x_lower_col],
                group[x_upper_col] - group[x_median_col],
            ])
            y_err = np.array([
                group["ga_median"] - group["ga_HDI_95_lower"],
                group["ga_HDI_95_upper"] - group["ga_median"],
            ])

            ax.errorbar(
                group[x_median_col],
                group["ga_median"],
                xerr=x_err,
                yerr=y_err,
                fmt="o",
                color=color,
                ecolor=color,
                elinewidth=1,
                capsize=0,
            )

            # Text labels, offset slightly above the point
            for _, row in group.iterrows():
                ax.text(
                    row[x_median_col] + 0.01,
                    row["ga_median"] + 0.01,
                    row["emerging_haplotype"],
                    fontsize=10,
                    color=color,
                    horizontalalignment="left",
                    verticalalignment="bottom",
                )

        ax.set_xlabel(x_title, fontsize=12)
        ax.tick_params(labelsize=12)
        ax.grid(False)

    # %%
    color_map = dict(zip(haplotype_colorings_domain, haplotype_colorings_range))

    fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=300, constrained_layout=True, sharey=True)

    plot_panel(
        axes[0],
        ga_by_antigenic_advance,
        "human_advance_median",
        "human_advance_lower",
        "human_advance_upper",
        "antigenic advance by human sera",
        color_map,
    )
    axes[0].set_ylabel("growth advantage", fontsize=12)

    plot_panel(
        axes[1],
        ga_by_antigenic_advance,
        "ferret_advance_median",
        "ferret_advance_lower",
        "ferret_advance_upper",
        "antigenic advance by ferret sera",
        color_map,
    )

    sns.despine()
    plt.savefig(output_path)
