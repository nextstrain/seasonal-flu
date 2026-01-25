#!/usr/bin/env python3

"""
Given a titers TSV we cross reference the strain names (virus strain & serum strain)
against the strain names in our curated metadata, as well as referencing a hardcoded
map of fauna strain names to EPI_ISL. The output titers TSV has strain names corrected
where possible and thus has improved strain matching for downstream analyses.
"""

import argparse
import sys
import csv
import json
from collections import defaultdict
from augur.io import read_metadata, open_file
from typing import Any
from datetime import date
import re


def read_titers(fname: str) -> list[dict[str, str]]:
    with open_file(fname) as f:
        reader = csv.DictReader(f, delimiter='\t')
        titers = [row for row in reader]
    print(f"Parsed titers TSV ({fname}). n={len(titers):,}", file=sys.stderr)
    return titers

def write_titers(fname: str, titers: list[dict[str, str]]) -> None:
    with open(fname, 'w', encoding='utf8') as f:
        writer = csv.DictWriter(f, fieldnames=titers[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(titers)

def read_fauna_map(fname: str) -> dict[str,str]:
    with open_file(fname) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fauna_map = {row['strain']:row['gisaid_epi_isl'] for row in reader}
    print(f"Parsed fauna strain map ({fname}). n={len(fauna_map):,}", file=sys.stderr)
    return fauna_map

UNKNOWN_YEAR = 'unknown'
def parse_year(date:str) -> int|str:
    try:
        return int(date[0:4])
    except ValueError:
        return UNKNOWN_YEAR

type Metadata = dict[str, tuple[str, int|str]]
def parse_metadata(fname: str) -> Metadata:
    df = read_metadata(fname)
    m:dict[str, tuple[str, int|str]] = dict([el[0], (el[1], parse_year(el[2]))] for el in zip(df.index, df['gisaid_epi_isl'], df['date'])) # type: ignore
    print(f"Parsed curated metadata ({fname}). n={len(m):,}", file=sys.stderr)
    return m

punctuation = ['/','_', '(', ')', ',', '-', ' ', ';', '.']
def simplify_strain(s: str) -> str:
    """
    Normalize string for more fuzzy-like matching. If a simplified titer strain
    matches a simplified metadata strain this is considered a "maybe" match.

    Simplification means: lowercase, punctuation removed, leading zeros removed from numbers,
    and some other ad-hoc fixes.
    """
    s = s.lower()
    # Remove leading zeros after slashes (e.g., /01 -> /1) before punctuation is stripped
    s = re.sub(r'/0+(\d)', r'/\1', s)
    for char in punctuation:
        s = s.replace(char, '')
    s = s.replace("a/a/", "a/")
    return s

def simplify_metadata(metadata:Metadata) -> dict[str, str]:
    """
    Generate a simplified-strain-name → true-strain-name map which which we will compare
    simplified titer names with. If there are collisions we drop them all and
    """
    mapping: defaultdict[str, list] = defaultdict(list)
    for true_strain in metadata.keys():
        mapping[simplify_strain(true_strain)].append(true_strain)
    simplified: dict[str,str] = {}
    dropped: list[list[str]] = []
    for simple_strain, true_strains in mapping.items():
        if len(true_strains)==1:
            simplified[simple_strain] = true_strains[0]
        else:
            dropped.append(true_strains)
    if dropped:
        count = sum([len(x) for x in dropped])
        print(f"[WARNING] In total we dropped {count:,} curated strains for simplified matching due to collisions once simplified", file=sys.stderr)
        # Comment out the following if you want to see them all, but I don't want to clutter the automated logs
        # print('[WARNING] Dropped strains: ' + ", ".join(["{" + ', '.join([f"{s!r}" for s in strains]) + "}" for strains in dropped]), file=sys.stderr)
    return simplified

type Match = tuple[str, str|None]

def match(titer_strain: str,
          cache: dict[str,str],
          fauna_map: dict[str,str],
          metadata: Metadata,
          epi_isl_to_strain: dict[str, str],
          simplified_metadata: dict[str,str],
          matching_if_simplified_strain_pairs: defaultdict[tuple[str,str], int]
          )->Match:
    """"
    Returns a tuple with the match result ("no", "maybe", "yes") and the strain name.
    If the result is "no" then the strain name wasn't matched and the returnder strain is the input strain.
    If the result is "yes" then the returned strain is (potentially) corrected such that it matches our metadata.
    If the result is "maybe" then the strain name is returned unchanged but there's a potential matching strain - see logs.
    """
    if titer_strain in cache:
        return ("yes", cache[titer_strain])

    if titer_strain in metadata:
        curated_strain = titer_strain # they're the same!
        cache[titer_strain] = curated_strain
        return ("yes", curated_strain)

    if titer_strain in fauna_map:
        epi_isl = fauna_map[titer_strain]
        if epi_isl in epi_isl_to_strain:
            curated_strain = epi_isl_to_strain[epi_isl]
            cache[titer_strain] = curated_strain
            return ("yes", curated_strain)
    
    # TODO XXX add a way to provide a hardcoded mapping of titer-strain-name to new-strain-name
    # and apply those matches here.

    # Don't actually match strains in simplified space, but log them into matching_if_simplified_strain_pairs
    strain_simple = simplify_strain(titer_strain)
    if strain_simple in simplified_metadata:
        curated_strain = simplified_metadata[strain_simple]
        matching_if_simplified_strain_pairs[(titer_strain, curated_strain)]+=1
        return ("maybe", None) # Don't return `curated_strain` - force us to deal with this manually

    return ("no", None)

def get_year(strain:str, metadata:Metadata|None=None, collapse_before:int=1990)->str:
    """
    Find the year of the strain, either by looking it up in metadata (if provided)
    or examining the strain name itself. Unknown years will be *UNKNOWN_YEAR*.
    Years prior to *collapse_before* will be collapsed into f"<={collapse_before}"
    """

    if metadata:
        assert strain in metadata, "Logic error - strain must be in metadata (get_year)"
        year = metadata[strain][1]
    else:
        try:
            year = strain.split('/')[-1]
            year = year.replace('-egg', '')
            if int(year) < 1960 or int(year) > date.today().year:
                # don't bother logging, I've checked and they're not solvable.
                year = UNKNOWN_YEAR
        except:
            print("[warning] Unknown strain year {strain!r}", strain, file=sys.stderr)
            year = UNKNOWN_YEAR

    if year==UNKNOWN_YEAR:
        return UNKNOWN_YEAR
    
    if int(year)<=collapse_before:
        return f"<={str(collapse_before)}"

    return str(year)


def compact_json(obj, indent=0):
    sp = '  ' * indent
    if isinstance(obj, dict):
        items = [f'{sp}  {json.dumps(k)}: {compact_json(v, indent + 1)}' for k, v in obj.items()]
        return '{\n' + ',\n'.join(items) + f'\n{sp}}}'
    elif isinstance(obj, list):
        return json.dumps(obj)  # lists on single line
    else:
        return json.dumps(obj)


def track_matches(matches,
        original_virus_strain: str,
        virus_match: Match,
        original_serum_strain: str,
        serum_match: Match,
        metadata: Metadata)->None:
    """
    Track (in *matches*) which (uncorrected, virus & serum) strains matched and which ones didn't,
    and the year they are from.
    We track (unique) individual strains separately to (unique) measurements (pairs of strains);
    for the latter, the year is taken as the virus_strain not the serum_strain.
    """
    virus_strain_year = get_year(virus_match[1], metadata) if virus_match[1] else get_year(original_virus_strain)
    serum_strain_year = get_year(serum_match[1], metadata) if serum_match[1] else get_year(original_serum_strain)

    # store whether the individual virus & serum strains matched:
    matches['strains'][original_virus_strain] = (virus_match[0], virus_strain_year)
    matches['strains'][original_serum_strain] = (serum_match[1], serum_strain_year)

    # and store whther the pair of strains (i.e. the measurement) matched
    if virus_match[0]=='no' or virus_match[0]=='no':
        measurement_match_result = 'no'
    elif virus_match[0]=='maybe' or virus_match[0]=='maybe':
        measurement_match_result = 'maybe'
    else:
        measurement_match_result = 'yes'
    matches['measurements'][(original_virus_strain, original_serum_strain)] = (measurement_match_result, virus_strain_year)



def calc_stats(matches: dict[str|tuple[str,str], tuple[str, str]])->dict[str,list[int]]:
    years = sorted(set([str(v[1]) for v in matches.values()]))
    by_year = {y: [0,0,0] for y in years} # format: yes matches, maybe matches, total
    by_year['total'] = [0,0,0]
    for match, year in matches.values():
        if match=='yes':
            by_year[str(year)][0] += 1
            by_year['total'][0] += 1
        elif match=='maybe':
            by_year[str(year)][1] += 1
            by_year['total'][1] += 1
        by_year[str(year)][2] += 1
        by_year['total'][2] += 1

    return by_year


def stats_str(a: int, b: int)->str:
    return f"{a:,}/{b:,} ({a/b*100:.1f}%)"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--titers", required=True, metavar='TSV',
                        help="Titers TSV with column 'virus_strain' and 'serum_strain'")
    parser.add_argument("--fauna-strain-map", required=True, metavar='TSV',
                        help="Folder which is expected to have subfolders 'h3n2' etc each with a fauna `fauna-metadata.tsv`")
    parser.add_argument("--metadata", required=True, metavar='TSV',
                        help="Curated metadata TSV")
    parser.add_argument("--output", required=True, metavar='TSV',
                        help="The corrected --titers input TSV")
    parser.add_argument("--stats", required=True, metavar='JSON',
                        help="matching stats")
    parser.add_argument("--stats-metadata", required=True, metavar='JSON_STRING',
                        help="JSON-formatted data to export in the stats JSON")

    args = parser.parse_args()

    # Keep track of strain names seen / remapped to speed things up when we encounter the same
    # strain names
    cache: dict[str, str] = {}

    _strain_matches: dict[str, tuple[str, str]] = {}
    _measurements_matches: dict[tuple[str,str], tuple[str, str]] = {}
    matches = {
        'strains': _strain_matches,
        'measurements': _measurements_matches,
    }

    fauna_map = read_fauna_map(args.fauna_strain_map)
    titers = read_titers(args.titers)
    metadata = parse_metadata(args.metadata)
    epi_isl_to_strain = {v[0]:k for k,v in metadata.items()}
    simplified_metadata = simplify_metadata(metadata)
    # Track pairs of strains (titer-strain, curated-strain) which would match if simplified
    matching_if_simplified_strain_pairs: defaultdict[tuple[str,str], int] = defaultdict(int)

    for row in titers:
        old_virus_strain = row['virus_strain']
        old_serum_strain = row['serum_strain']

        virus_match = match(old_virus_strain, cache, fauna_map, metadata, epi_isl_to_strain, simplified_metadata, matching_if_simplified_strain_pairs)
        serum_match = match(old_serum_strain, cache, fauna_map, metadata, epi_isl_to_strain, simplified_metadata, matching_if_simplified_strain_pairs)

        track_matches(matches, old_virus_strain, virus_match, old_serum_strain, serum_match, metadata)

        row['virus_strain'] = virus_match[1] or old_virus_strain
        row['serum_strain'] = serum_match[1] or old_serum_strain

        row['virus_strain_match'] = virus_match[0]
        row['serum_strain_match'] = serum_match[0]

    write_titers(args.output, titers)

    # Log out titer strains which match curated strains when simplified, as we want to either
    # update the curated metadata to match the titers or remap the titer strain to match.
    # (Remapping titer strains like this is not yet implemented TODO XXX)
    for pair, count in matching_if_simplified_strain_pairs.items():
        print(f"[NEEDS CHECKING] Titer strain {pair[0]!r} would match curated strain {pair[1]!r} if simplified. Strain appears {count:,} times", file=sys.stderr)

    unique_strain_stats = calc_stats(matches['strains'])  # type: ignore
    measurement_stats = calc_stats(matches['measurements'])  # type: ignore

    print(f"Total matching (unique) strains:    {stats_str(unique_strain_stats['total'][0], unique_strain_stats['total'][2])}", file=sys.stderr)
    print(f"Total matching (unique) measurement pairs: {stats_str(measurement_stats['total'][0], measurement_stats['total'][2])}", file=sys.stderr)

    print(f"Potential matching (unique) strains:    {stats_str(unique_strain_stats['total'][1], unique_strain_stats['total'][2])}", file=sys.stderr)
    print(f"Potential matching (unique) measurement pairs: {stats_str(measurement_stats['total'][1], measurement_stats['total'][2])}", file=sys.stderr)


    stats = {
        **json.loads(args.stats_metadata),
        '__help__': "List format is <yes-matches> <maybe-matches> <total>",
        'strains': unique_strain_stats,
        'measurements': measurement_stats,
    }
    with open(args.stats, 'w') as fh:
        print(compact_json(stats), file=fh)
