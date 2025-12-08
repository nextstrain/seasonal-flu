#!/usr/bin/env python3

"""
TKTK
"""


import argparse
import sys
import zstandard as zstd
import gzip
import lzma
from collections import defaultdict
import io
import csv
import glob
import os
from augur.io import read_metadata
import matplotlib.pyplot as plt
import matplotlib

SUBTYPES = ['h3n2', 'h1n1pdm', 'vic', 'yam']

# SUBTYPES = ['h3n2']
UNKNOWN_YEAR = 'unknown'

type Subtype = str
type Contrib = str
type Strains = set[str]
type Year = str # '2014' etc, with 'unknown' the null
type TiterInfo = defaultdict[Subtype, defaultdict[Contrib, defaultdict[Year, Strains]]]
type MatchInfo = tuple[int, int] # match count, total number strains
type TiterMatches = defaultdict[Subtype, defaultdict[Contrib, defaultdict[Year, MatchInfo]]]


def open_file(file_path: str):
    """
    Open a file, decompressing if needed.
    Supports .zst, .gz, and .xz files.
    Returns a text-mode file handle.
    """
    if file_path.endswith('.zst'):
        dctx = zstd.ZstdDecompressor()
        fh = open(file_path, 'rb')
        reader = dctx.stream_reader(fh)
        # Wrap in TextIOWrapper to get text mode
        return io.TextIOWrapper(reader, encoding='utf-8')
    elif file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt', encoding='utf-8')

    elif file_path.endswith('.xz'):
        return lzma.open(file_path, 'rt', encoding='utf-8')
    else:
        return open(file_path, 'r')


def read_titer_file(fname: str) -> list[tuple[str,str]]:
    virus_serum_tuples: list[tuple[str,str]] = []
    with open_file(fname) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            virus_serum_tuples.append((row['virus_strain'], row['serum_strain']))

    return virus_serum_tuples

def extract_year(strain:str) -> str:
    try:
        year = strain.split('/')[-1]
        year = year.replace('-egg', '')
        if int(year) < 1960 or int(year)>2026:
            # don't bother logging, I've checked and they're not solvable.
            # BUT! If we cross-reference with fauna they might be?
            return UNKNOWN_YEAR
        return year
    except:
        print("Unknown strain year", strain, file=sys.stderr)
        return UNKNOWN_YEAR

def read_titers(dir: str)-> TiterInfo:
    """
    Read all the titers files and collapse them to a simple set of unique strains
    for each subtype, for each collaborating center.
    Don't differentiate between virus strain and serum strain.
    """
    data: TiterInfo = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    for subtype in SUBTYPES:
        for file in glob.glob(os.path.join(dir, subtype, '*.tsv.gz')):
            contributor = file.split(subtype)[1].replace('/', '').split('_')[0]
            pairs = read_titer_file(file)
            for pair in pairs:
                for el in pair:
                    data[subtype][contributor][extract_year(el)].add(el)
    print("-------------- titers ---------------", file=sys.stderr)
    for subtype in SUBTYPES:
        for contrib,strains_per_year in data[subtype].items():
            num_strains = sum([len(strains) for strains in strains_per_year.values()])
            years = sorted(list(strains_per_year.keys()))
            print(f"{subtype:<10}{contrib:<8} n(unique strains)={num_strains:,}", file=sys.stderr)
    return data


def parse_metadata(name: str, dir: str, filename: str, epi_isl_key:str = "gisaid_epi_isl") -> defaultdict[Subtype, list[tuple[str,str]]]:
    data: defaultdict[Subtype, list[tuple[str,str]]] = defaultdict(list)
    for subtype in SUBTYPES:
        df = read_metadata(os.path.join(dir, subtype, filename))
        data[subtype] = list(zip(df.index, df[epi_isl_key]))
    print(f"------ parsing {name} metadata  ------", file=sys.stderr)
    for subtype in SUBTYPES:
        print(f"{subtype:<10} n(strains/epi-isls)={len(data[subtype]):,}", file=sys.stderr)
    return data

def strain_set(tuple_data: defaultdict[Subtype, list[tuple[str,str]]]) -> defaultdict[Subtype, Strains]:
    data: defaultdict[Subtype, Strains] = defaultdict(set)
    for subtype in SUBTYPES:
        data[subtype] = set([pair[0]for pair in tuple_data[subtype]])
    return data



def compute_simple_matches(name: str,
        strain_data: defaultdict[Subtype, Strains],
        titers: TiterInfo) -> TiterMatches:
    print(f"------ {name} metadata vs titers (simple) matching --------", file=sys.stderr)
    data: TiterMatches = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for subtype in SUBTYPES:
        metadata_strains = strain_data[subtype]
        for contrib in titers[subtype]:
            for year, titer_strains in titers[subtype][contrib].items():
                num_matches = len(titer_strains & metadata_strains)
                num_titer_strains = len(titer_strains)
                data[subtype][contrib][year] = (num_matches, num_titer_strains)
                print(f"{subtype:<10}{contrib:<8}{year:<9} {num_matches:,}/{num_titer_strains:,} ({int(num_matches*100/num_titer_strains)}%)", file=sys.stderr)
    return data

punctuation = ['/','_', '(', ')', ',', '-', ' ', ';', '.']
def simplify_strain(s: str) -> str:
    """Normalize string for more fuzzy-like matching"""
    s = s.lower()
    for char in punctuation:
        s = s.replace(char, '')
    s = s.replace("a/a/", "a/")
    return s

def compute_complex_matches(
        curated_metadata: defaultdict[Subtype, list[tuple[str,str]]],
        fauna_metadata: defaultdict[Subtype, list[tuple[str,str]]],
        titers: TiterInfo) -> TiterMatches:
    print(f"------ curated metadata vs titers (complex) matching --------", file=sys.stderr)
    data: TiterMatches = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for subtype in SUBTYPES:
        curated_metadata_strains   = set([pair[0]for pair in curated_metadata[subtype]])
        fauna_strain_to_epi_isl    = dict(fauna_metadata[subtype])
        epi_isl_to_curated_strain  = dict([[x[1],x[0]] for x in curated_metadata[subtype]])
        normalized_curated_strains = {simplify_strain(s) for s in curated_metadata_strains}

        for contrib in titers[subtype]:
            for year, titer_strains in titers[subtype][contrib].items():
                num_titer_strains = len(titer_strains)
                num_matches:int = 0
                for titer_strain in titer_strains:
                    # first do the same as the simple matching
                    if titer_strain in curated_metadata_strains:
                        num_matches+=1
                        continue

                    # If it's a fauna match & we can use EPI ISL linkage then do so
                    if epi_isl:=fauna_strain_to_epi_isl.get(titer_strain, None):
                        if curated_strain:=epi_isl_to_curated_strain.get(epi_isl, None):
                            num_matches+=1
                            continue
                    
                    # Normalize
                    if simplify_strain(titer_strain) in normalized_curated_strains:
                        num_matches+=1
                        continue

                data[subtype][contrib][year] = (num_matches, num_titer_strains)
                print(f"{subtype:<10}{contrib:<8}{year:<9} {num_matches:,}/{num_titer_strains:,} ({int(num_matches*100/num_titer_strains)}%)", file=sys.stderr)
    return data


def plot_matches(fauna_matches: TiterMatches,
                 curated_matches: TiterMatches,
                 curated_matches_complex: TiterMatches,
                 output_file: str):
    """
    Create small-multiples visualization of fauna matching percentages over time.

    Creates one subplot per subtype, with each contributor as a separate line.
    Two rows: fauna matches (top) and curated matches (bottom)
    X-axis: years
    Y-axis: percentage of strains matched (a/b*100)
    """
    # Set up the figure with 3 rows and columns for each subtype
    n_subtypes = len(SUBTYPES)
    fig, axes = plt.subplots(3, n_subtypes, figsize=(5*n_subtypes, 12), squeeze=False)

    # Collect all unique contributors across all subtypes and both datasets
    _all_contributors: set[str] = set()
    for subtype in SUBTYPES:
        _all_contributors.update(fauna_matches[subtype].keys())
        _all_contributors.update(curated_matches[subtype].keys())
    all_contributors: list[str] = sorted(_all_contributors)

    n_contributors = len(all_contributors)
    cmap = plt.cm.tab10
    contributor_colors = {contrib: cmap(i % cmap.N)
                         for i, contrib in enumerate(all_contributors)}

    # Plot both fauna and curated matches
    datasets = [
        (0, fauna_matches, "Fauna metadata - titers match %"),
        (1, curated_matches, "Curated metadata - titers match %"),
        (2, curated_matches_complex, "Curated metadata - titers match %\nwith remapping layer"),
    ]

    for row_idx, match_data, y_label in datasets:
        for col_idx, subtype in enumerate(SUBTYPES):
            ax = axes[row_idx, col_idx]

            # Get data for this subtype
            subtype_data = match_data[subtype]

            # Collect all unique years across all contributors for this subtype
            all_years_set: set[str] = set()
            for year_data in subtype_data.values():
                all_years_set.update(year_data.keys())

            # Sort years, ensuring 'unknown' comes at the end
            numeric_years = sorted([y for y in all_years_set if y != UNKNOWN_YEAR])
            all_years = numeric_years + ([UNKNOWN_YEAR] if UNKNOWN_YEAR in all_years_set else [])

            # Track overall percentages for y-axis ticks
            overall_percentages = []
            overall_colors = []
            overall_contributors = []

            # Plot each contributor as a separate line
            for contributor, year_data in subtype_data.items():
                # Calculate percentages and totals only for years where this contributor has data
                contributor_years = []
                percentages = []
                totals = []
                for year in all_years:
                    if year in year_data:
                        matches, total = year_data[year]
                        contributor_years.append(year)
                        if total > 0:
                            percentages.append((matches / total) * 100)
                        else:
                            percentages.append(0)
                        totals.append(total)

                if not contributor_years:
                    continue

                # Calculate overall percentage across all years (including unknown)
                total_matches = sum(year_data[year][0] for year in year_data.keys())
                total_strains = sum(year_data[year][1] for year in year_data.keys())
                overall_pct = (total_matches / total_strains * 100) if total_strains > 0 else 0

                # Get the predefined color for this contributor
                line_color = contributor_colors[contributor]

                # Plot the line without markers
                line = ax.plot(contributor_years, percentages, label=contributor,
                              linewidth=0.5, color=line_color)

                # Store overall percentage, color, and contributor name for y-axis ticks
                overall_percentages.append(overall_pct)
                overall_colors.append(line_color)
                overall_contributors.append(contributor)

                # Add scatter points with sizes proportional to total strains
                # Scale the marker size: base size + proportional to total
                marker_sizes = [t / 10 for t in totals]  # Adjust divisor to scale appropriately
                ax.scatter(contributor_years, percentages, s=marker_sizes, color=line_color,
                          alpha=0.7, zorder=5)

            # Customize the subplot
            # Only show title on the top row
            if row_idx == 0:
                ax.set_title(f'{subtype.upper()}', fontsize=14, fontweight='bold')
            # Only show y-axis label on the first (leftmost) subplot of each row
            if col_idx == 0:
                ax.set_ylabel(y_label, fontsize=12)
            ax.set_ylim(0, 100)

            # Get the default y-ticks that matplotlib would use (for grid lines)
            existing_yticks = list(ax.get_yticks())

            # Draw grid at the original tick positions (both x and y)
            ax.grid(True, alpha=0.3)

            # Now add overall percentage ticks as minor ticks (won't create grid lines)
            # First, keep the major ticks as they were
            ax.set_yticks(existing_yticks)
            # Add overall percentages as minor ticks (tick marks only)
            ax.set_yticks(overall_percentages, minor=True)
            # Don't show labels on the minor ticks themselves
            ax.yaxis.set_tick_params(which='minor', labelleft=False)

            # Add custom text labels for overall percentages positioned to the right of the y-axis
            for contributor, overall_pct, color in zip(overall_contributors, overall_percentages, overall_colors):
                # Position text just to the right of the y-axis (inside the plot area)
                ax.text(0.01, overall_pct, f'{contributor} = {overall_pct:.1f}%',
                       transform=ax.get_yaxis_transform(),
                       ha='left', va='center',
                       color=color, fontweight='bold',
                       fontsize=9)

            # Set the x-axis limits and tick positions
            # This ensures the categorical axis follows our year ordering
            ax.set_xlim(-0.5, len(all_years) - 0.5)

            # Set x-axis to only show every 5th year
            # Build tick list ensuring 'unknown' is at the end if it exists
            tick_positions = []
            tick_labels = []
            for idx, year in enumerate(all_years):
                if year == UNKNOWN_YEAR or (year != UNKNOWN_YEAR and int(year) % 5 == 0):
                    tick_positions.append(idx)
                    tick_labels.append(year)

            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels)

            # Rotate x-axis labels for better readability
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save as high-resolution PNG
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to {output_file}", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--titers", required=True, metavar='DIR',
                        help="Folder which is expected to have subfolders 'h3n2' etc each with titers TSVs in them")
    parser.add_argument("--fauna", required=True, metavar='DIR',
                        help="Folder which is expected to have subfolders 'h3n2' etc each with a fauna `metadata.tsv.xz`")
    parser.add_argument("--curated", required=True, metavar='DIR',
                        help="Folder which is expected to have subfolders 'h3n2' etc each with a fauna `metadata.tsv.xz`")
    parser.add_argument("--output", default="fauna_matches.png",
                        help="Output file for the visualization (default: fauna_matches.png)")
    # parser.add_argument("--query-strains", required=True, help="Text file with query strain names (e.g. include.txt)")
    # parser.add_argument("--strain-map", required=False, nargs='*', help="Hardcoded map of query strain (fauna) to new curated strain")
    # parser.add_argument("--no-fuzz", required=False,action='store_true', help="Skip fuzzing for efficiency reasons")

    args = parser.parse_args()

    titers = read_titers(args.titers)

    fauna = parse_metadata("fauna", args.fauna, "metadata.tsv.xz")
    curated = parse_metadata("curated", args.curated, "metadata.tsv")
    fauna_strains = strain_set(fauna)
    curated_strains = strain_set(curated)

    # simple (identical) matching
    fauna_matches = compute_simple_matches("fauna", fauna_strains, titers)
    curated_matches = compute_simple_matches("curated", curated_strains, titers)

    # complex matching 
    curated_matches_complex = compute_complex_matches(curated, fauna, titers)

    plot_matches(fauna_matches, curated_matches, curated_matches_complex, args.output)


