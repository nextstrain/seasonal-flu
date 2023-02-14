"""
python3 plot_forecasts.py \
    --tree auspice/flu.json \
    --models cTiter_x-ne_star ne_star-lbi \
    --frequencies cTiter_x-ne_star_tip-frequencies.json ne_star-lbi_tip-frequencies.json \
    --clades 'A1b/131K' 'A1b/135K'
"""
import argparse
from datetime import datetime
import matplotlib as mpl
mpl.use('agg')
import matplotlib.dates as mdates
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


fontsize = 14
matplotlib_params = {
    'axes.labelsize': fontsize,
    'font.size': fontsize,
    'legend.fontsize': 12,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'text.usetex': False,
    'figure.figsize': [8, 6],
    'savefig.dpi': 300,
    'figure.dpi': 300,
    'text.usetex': False
}
plt.rcParams.update(matplotlib_params)

old_epoch = '0000-12-31T00:00:00'
mdates.set_epoch(old_epoch)

# From axes.prop_cycle described at https://matplotlib.org/users/customizing.html
line_colors = [
    '#ff7f0e',
    '#2ca02c',
    '#9467bd',
    '#e377c2',
    '#8c564b',
    '#d62728',
    '#7f7f7f',
    '#bcbd22',
    '#17becf'
]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--forecasts", nargs="+", help="tables of tip frequencies annotated by tip name, clade, model, and date of model")
    parser.add_argument("--clades", nargs="+", help="names of clades to plot")
    parser.add_argument("--models", nargs="+", help="names of models as hyphenated lists of fitness metrics (e.g., ne_star-lbi)")
    parser.add_argument("--model-names", nargs="+", help="human-readable model names for each provided model (e.g., mutational load and LBI)")
    parser.add_argument("--groupby", default="date", help="column from the input table to group forecasts for plotting as separate lines (e.g., by date, sample, etc.)")
    parser.add_argument("--drop-samples", nargs="+", help="samples to drop from plots (e.g., due to incorrect phylogenetic nesting of clades, etc.)")
    parser.add_argument("--height-per-row", type=float, default=1.25, help="height per row in inches")
    parser.add_argument("--output", help="forecast plot")

    args = parser.parse_args()

    # Define primary model as the first model.
    primary_model = args.models[0]

    # Map model names to human-readable names.
    model_names_by_model = dict(zip(args.models, args.model_names))
    color_by_model = dict(zip(args.models, line_colors[:len(args.models)]))

    # Load frequencies.
    frequencies = pd.concat([pd.read_csv(forecast, sep="\t") for forecast in args.forecasts])

    # Filter samples.
    if args.drop_samples is not None and len(args.drop_samples) > 0:
        frequencies = frequencies[~frequencies[args.groupby].isin(args.drop_samples)]

    # Group frequencies by clade, model, and date.
    df = frequencies.groupby(["clade", "pivot", "model", args.groupby, "observed"])["frequency"].sum().reset_index()

    # Filter clades.
    df = df[df["clade"].isin(args.clades)].copy()

    # Get pivots.
    numeric_pivots = df["pivot"].drop_duplicates().sort_values().values
    timepoint = df[df["observed"]]["pivot"].max()
    timepoint_index = list(numeric_pivots).index(timepoint) + 1

    years = YearLocator()
    yearsFmt = DateFormatter('%Y')
    months = MonthLocator(range(1, 13), bymonthday=1, interval=3)
    monthsFmt = DateFormatter("%b")

    offset = datetime(2000,1,1).toordinal()
    pivots = [offset + (x - 2000) * 365.25 for x in numeric_pivots[:timepoint_index]]
    projected_pivots = numeric_pivots[timepoint_index:]
    ordinal_projected_pivots = [offset + (x - 2000) * 365.25 for x in projected_pivots]

    today = datetime.today().toordinal()

    # Plot frequencies.
    fig, axes = plt.subplots(int(np.ceil(len(args.clades) / 2.0)), 2, figsize=(16, args.height_per_row * len(args.clades)), squeeze=False,
                             sharex=True, sharey=True)

    delta_frequency_records = []
    for ci, clade_name in enumerate(args.clades):
        ax = axes.flatten()[ci]

        clade_df = df[df["clade"] == clade_name].copy()

        # Plot per model per groupby column per type
        for (model, group), clade_model_df in clade_df.groupby(["model", args.groupby]):
            if model not in args.models:
                continue

            observed_clade_model_df = clade_model_df[clade_model_df["observed"]]
            predicted_clade_model_df = clade_model_df[~clade_model_df["observed"]]
            delta_frequency = predicted_clade_model_df["frequency"].values[-1] - observed_clade_model_df["frequency"].values[-1]

            # Clades that only change by +/- 5% frequency neither grow nor decline.
            if np.abs(delta_frequency) < 0.05:
                status = "persist"
            elif delta_frequency < 0:
                status = "decline"
            else:
                status = "grow"

            delta_frequency_records.append({
                "clade": clade_name,
                args.groupby: group,
                "model": model,
                "delta_frequency": delta_frequency,
                "status": status
            })

            # Plot observed clade frequencies up to the present.
            markersize = 6
            linewidth = 1.0
            alpha = 0.2
            zorder = -1

            ax.plot(pivots, observed_clade_model_df["frequency"], "o-", c="#1f77b4", alpha=alpha, linewidth=linewidth, zorder=zorder)
            color = color_by_model[model]

            # Plot projected clade frequencies up to one year in the future.
            ax.plot(
                ordinal_projected_pivots,
                predicted_clade_model_df["frequency"],
                "o--",
                linewidth=linewidth,
                markersize=markersize,
                color=color,
                zorder=zorder,
                alpha=alpha
            )

        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(['{:3.0f}%'.format(x*100) for x in [0, 0.2, 0.4, 0.6, 0.8, 1.0]])
        ax.tick_params(axis='x', which='major', labelsize=fontsize, pad=20)
        ax.tick_params(axis='x', which='minor', pad=7)
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.xaxis.set_minor_locator(months)
        ax.xaxis.set_minor_formatter(monthsFmt)

        # Despine axes.
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_xlabel("Date")
        ax.set_ylabel("Frequency")

    delta_frequency_df = pd.DataFrame(delta_frequency_records)
    delta_frequency_df.to_csv("delta_frequency.tsv", sep="\t", index=False)

    proportion_by_status = (delta_frequency_df.groupby([
        "clade",
        "status"
    ])["sample"].count().reset_index().pivot(index="clade", columns=["status"], values=["sample"]).fillna(0) / delta_frequency_df["sample"].drop_duplicates().shape[0]) * 100

    for ci, clade_name in enumerate(args.clades):
        ax = axes.flatten()[ci]
        clade_proportions = proportion_by_status.query(f"clade == '{clade_name}'")
        # clade_text = f"{clade_name} ({clade_proportions[('sample', 'decline')][0]}% decline, {clade_proportions[('sample', 'persist')][0]}% persist, {clade_proportions[('sample', 'grow')][0]}% grow)"
        clade_text = f"{clade_name}"
        ax.text(
            0.025,
            0.95,
            clade_text,
            {"fontsize": 14},
            transform=ax.transAxes,
            horizontalalignment="left",
            verticalalignment="center"
        )

    if len(args.clades) < axes.flatten().shape[0]:
        axes.flatten()[-1].axis("off")

    patches = [mpatches.Patch(color="#1f77b4", label="Observed")]
    for model in args.models:
        patches.append(
            mpatches.Patch(
                color=color_by_model[model],
                label="Predicted by " + model_names_by_model[model]
            )
        )

    axes.flatten()[-1].legend(
        handles=patches,
        loc="upper right",
        frameon=False,
        fancybox=False
    )

    plt.tight_layout(w_pad=1.5)
    plt.savefig(args.output, bbox_inches="tight")
    plt.close()
