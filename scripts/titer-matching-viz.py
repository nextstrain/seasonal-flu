#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
import json
from typing import TypedDict
import sys
import argparse
import sys

YearStats = dict[str, tuple[int, int, int]]  # {year: (matched, maybe_matched, total)}

class Stats(TypedDict):
    lineage: str
    center: str
    passage: str
    assay: str
    strains: YearStats
    measurements: YearStats

def read_json(fname:str)->Stats:
    with open(fname) as fh:
        return json.load(fh)



def order_years(stats: list[Stats]) -> tuple[list[str], list[str]]:
    """
    Produces the x-values for visualisation so that we have consistency across plots
    """
    observed: set[str] = set()
    xvals: list[str] = []
    xlabels: list[str] = []
    
    for s in stats:
        observed.update(s['strains'].keys())
        observed.update(s['measurements'].keys())
    
    observed.remove('total')

    if 'unknown' in observed:
        observed.remove('unknown')
        xvals.append('unknown')
        xlabels.append('unknown')
        # add empty x-value/label to space the unknown a little more to the left
        xvals.append('<dummy>')
        xlabels.append('')
    try:
        min_cat = [v for v in observed if v.startswith('<=')][0]
        observed.remove(min_cat)
        xvals.append(min_cat)
        xlabels.append(min_cat)
    except IndexError:
        pass

    ordered = sorted(observed)
    years = [str(year) for year in range(int(ordered[0]), int(ordered[-1])+1)]

    xvals.extend(years)
    xlabels.extend([y if idx!=0 and idx%5==0 else '' for idx,y in enumerate(years) ])
    xlabels[-1] = xvals[-1] # ensure last (most recent) value is shown
    return (xvals, xlabels)


def plot_matches(stats: list[Stats], output_file: str):
    """
    Create small-multiples visualization of fauna matching percentages over time.

    Creates one subplot per subtype, with each contributor as a separate line.
    Two rows: fauna matches (top) and curated matches (bottom)
    X-axis: years
    Y-axis: percentage of strains matched (a/b*100)
    """

    lineages = set([s['lineage'] for s in stats])
    assert len(lineages)==1
    lineage = list(lineages)[0]
    centers = list(set([s['center'] for s in stats]))
    stats_per_center = [[s for s in stats if s['center']==center] for center in centers]
    max_stats_per_center = max([len(l) for l in stats_per_center])
    xvals, xlabels = order_years(stats)
    x_indices = list(range(len(xvals)))

    c_strain = '#c51b8a'
    c_measurement = '#2c7fb8'

    fig, axes = plt.subplots(len(centers), max_stats_per_center, figsize=(max_stats_per_center*5, len(centers) * 4), squeeze=False)

    circle_scalar = 1/10
    circla_alpha = 0.5

    for row_idx, center in enumerate(centers):
        for col_idx, data in enumerate(stats_per_center[row_idx]):
            ax = axes[row_idx, col_idx]

            percentages = {
                'strains': {x:data['strains'][x][0]/data['strains'][x][2]*100 for x in xvals if x in data['strains']},
                'strains_maybe': {x:data['strains'][x][1]/data['strains'][x][2]*100 for x in xvals if x in data['strains']},
                'measurements': {x:data['measurements'][x][0]/data['measurements'][x][2]*100 for x in xvals if x in data['measurements']},
                'measurements_maybe': {x:data['measurements'][x][1]/data['measurements'][x][2]*100 for x in xvals if x in data['measurements']},
            }

            # Extract y-values and sizes in the order of years
            strains_y = [percentages['strains'].get(x) for x in xvals]
            strains_maybe_y = [percentages['strains_maybe'].get(x) for x in xvals] # maybe: potential matches
            strains_sizes = [data['strains'][x][0]*circle_scalar if x in data['strains'] else 0 for x in xvals]
            strains_maybe_sizes = [data['strains'][x][1]*circle_scalar if x in data['strains'] else 0 for x in xvals]
            measurements_y = [percentages['measurements'].get(x) for x in xvals]
            measurements_maybe_y = [percentages['measurements_maybe'].get(x) for x in xvals]  # maybe: potential matches
            measurements_sizes = [data['measurements'][x][0]*circle_scalar if x in data['measurements'] else 0 for x in xvals]
            measurements_maybe_sizes = [data['measurements'][x][1]*circle_scalar if x in data['measurements'] else 0 for x in xvals]

            # Maybe (potential) matches of strains
            ax.plot(x_indices, strains_maybe_y,
                    linewidth=0.5, linestyle='--', color=c_strain)
            ax.scatter(x_indices, strains_maybe_y,
                    s=strains_maybe_sizes, color=c_strain, alpha=circla_alpha, zorder=5, clip_on=False)
        
            # True matches of strains
            ax.plot(x_indices, strains_y,
                    linewidth=0.5, color=c_strain)
            ax.scatter(x_indices, strains_y,
                    s=strains_sizes, color=c_strain, alpha=circla_alpha, zorder=5, clip_on=False)

            # Maybe (potential) matches of measurements
            ax.plot(x_indices, measurements_maybe_y,
                    linewidth=0.5, linestyle='--', color=c_measurement)
            ax.scatter(x_indices, measurements_maybe_y,
                    s=measurements_maybe_sizes, color=c_measurement, alpha=circla_alpha, zorder=5, clip_on=False)

            # True matches of measurements
            ax.plot(x_indices, measurements_y,
                    linewidth=0.5, color=c_measurement)
            ax.scatter(x_indices, measurements_y,
                    s=measurements_sizes, color=c_measurement, alpha=circla_alpha, zorder=5, clip_on=False)

            # Add horizontal dashed lines for overall percentages
            for idx in [0,1]:
                total_strains = data['strains']['total']
                total_measurements = data['measurements']['total']
                pct_strains = total_strains[idx] / total_strains[2] * 100
                pct_measurements = total_measurements[idx] / total_measurements[2] * 100

                if idx==0:
                    ax.axhline(pct_strains, linestyle='-', color=c_strain, alpha=0.7)
                    ax.axhline(pct_measurements, linestyle='-', color=c_measurement, alpha=0.7)
                    # Add text labels for total percentages in bottom left
                    ax.text(0.02, 0.42, f'Unique strains matched:', transform=ax.transAxes, color=c_strain, fontsize=12)
                    ax.text(0.02, 0.34, f'{total_strains[0]:,} / {total_strains[2]:,} ({pct_strains:.1f}%)', transform=ax.transAxes, color=c_strain, fontsize=12)
                    ax.text(0.02, 0.18, f'Unique measurements matched:', transform=ax.transAxes, color=c_measurement, fontsize=12)
                    ax.text(0.02, 0.10, f'{total_measurements[0]:,} / {total_measurements[2]:,} ({pct_measurements:.1f}%)', transform=ax.transAxes, color=c_measurement, fontsize=12)
                if idx==1:
                    ax.text(0.02, 0.26, f'maybes = {total_strains[1]:,} ({pct_strains:.1f}%)', transform=ax.transAxes, color=c_strain, fontsize=12)
                    ax.text(0.02, 0.02, f'maybes = {total_measurements[1]:,} ({pct_measurements:.1f}%)', transform=ax.transAxes, color=c_measurement, fontsize=12)

            ax.set_xticks(x_indices)
            ax.set_xticklabels(xlabels, rotation=45, ha='right')
            ax.set_xlim(-0.5, len(x_indices) - 0.5)

            if col_idx == 0:
                ax.set_ylabel(center, fontsize=16)

            ax.set_title(f"{lineage} | {center} | {data['passage']} | {data['assay']}")
            ax.set_ylim(0, 100)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

    # Add figure title
    fig.suptitle(f"Lineage: {lineage}", fontsize=16, fontweight='bold')

    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    # Save as high-resolution PNG
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to {output_file}", file=sys.stderr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", nargs="+", required=True, metavar='JSON', help="stats JSONs")
    parser.add_argument("--output", required=True, metavar='PNG', help="output viz")
    args = parser.parse_args()

    try:
        stats = [read_json(s) for s in sorted(args.stats)]
        plot_matches(stats, args.output)
    except Exception as e:
        print(f"Visualisation script failed with error:", file=sys.stderr)
        print(e, file=sys.stderr)
        print(f"Exiting with code 0 so that automated pipelines don't die", file=sys.stderr)
        with open(args.output, 'w') as fh:
            print("Script failed!", file=fh)

