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
from typing import TypeAlias, cast
import dataclasses
from datetime import date
import re

HARDCODED_CASE_CHANGES: dict[str,str] = {
    # (incorrect) titer strain substring, corrected strain substring
    "/Americansamoa/":  "/AmericanSamoa/",
    "/Angthong/":  "/AngThong/",
    "/AnNahdah/":  "/An_Nahdah/",
    "/BayofPlenty/":  "/BayOfPlenty/",
    "/Chiangmai/":  "/ChiangMai/",
    "/Chiangrai/":  "/ChiangRai/",
    "/ChristChurch/":  "/Christchurch/",
    "/DistrictofColumbia/":  "/DistrictOfColumbia/",
    "/Dominicanrepublic/":  "/DominicanRepublic/",
    "/GansuBaiyin/":  "/Gansu-Baiyin/",
    "/Hawkesbay/":  "/HawkesBay/",
    "/IledeFrance/":  "/IleDeFrance/",
    "/Kyrgyzstan-Bishkek/":  "/Kyrgyzstan_Bishkek/",
    "/Kyrgyzstan-Osh/":  "/Kyrgyzstan_Osh/",
    "/Naknonsithammarat/":  "/NaknonSiThammarat/",
    "/NetherLands/":  "/Netherlands/",
    "/Newcaledonia/":  "/NewCaledonia/",
    "/Newhampshire/":  "/NewHampshire/",
    "/Newjersey/":  "/NewJersey/",
    "/Newmexico/":  "/NewMexico/",
    "/Newyork/":  "/NewYork/",
    "/Northcarolina/":  "/NorthCarolina/",
    "/Northwestauckland/":  "/NorthWestAuckland/",
    "/Papuanewguinea/":  "/PapuaNewGuinea/",
    "/PaysdeLoire/":  "/PaysDeLoire/",
    "/Phranakhonsiayutthaya/":  "/PhraNakhonSiAyutthaya/",
    "/Saopaulo/":  "/SaoPaulo/",
    "/Solomonislands/":  "/SolomonIslands/",
    "/South-Africa/":  "/SouthAfrica/",
    "/Southafrica/":  "/SouthAfrica/",
    "/Southauckland/":  "/SouthAuckland/",
    "/southAuckland/":  "/SouthAuckland/",
    "/Southaustralia/":  "/SouthAustralia/",
    "/Southcarolina/":  "/SouthCarolina/",
    "/Srilanka/":  "/SriLanka/",
    "/Timor-Leste/":  "/TimorLeste/",    
}

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
    
def read_hardcoded_strain_map(fname: str) -> dict[str,str]:
    with open_file(fname) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader) # header
        strain_map = {row[0]:row[1] for row in reader}
    print(f"Parsed hardcoded strain map ({fname}). n={len(strain_map):,}", file=sys.stderr)
    return strain_map

UNKNOWN_YEAR = 'unknown'
def parse_year(date:str) -> int|str:
    try:
        return int(date[0:4])
    except ValueError:
        return UNKNOWN_YEAR

Metadata: TypeAlias = dict[str, tuple[str, int|str]]
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
        print(f"[SIMPLIFIED MATCHING] When creating a lookup of simplified strain name → actual strain name, we dropped {count:,} strains due to collisions once simplified", file=sys.stderr)
        # Comment out the following if you want to see them all, but I don't want to clutter the automated logs
        # print('[WARNING] Dropped strains: ' + ", ".join(["{" + ', '.join([f"{s!r}" for s in strains]) + "}" for strains in dropped]), file=sys.stderr)
    return simplified

@dataclasses.dataclass(frozen=True)
class Match:
    in_metadata: bool # is the (new) strain name found in metadata?    
    original_name: str
    new_name: str
    new_name_matching_when_simplified: bool|str = False # if set implies found_in_metadata=False
    name_changed_via_pattern_match: bool = False

def match(titer_strain: str,
          fauna_map: dict[str,str],
          hardcoded_map: dict[str,str],
          metadata: Metadata,
          epi_isl_to_strain: dict[str, str],
          simplified_metadata: dict[str,str],
          matching_if_simplified_strain_pairs: defaultdict[tuple[str,str], int],
          hardcoded_strains_not_found_in_metadata: defaultdict[tuple[str,str], int],
          )->Match:
    """"
    Matches **titer_strain** (a virus_strain or serum_strain) against metadata using
    a few approaches: Returns a Match object detailing the result, and including the
    (potentially new) strain_name to use in downstream code.
    """
    # hardcoded map takes precedence over everything, and has the special (strange?)
    # case of a hardcoded-corrected strain name which is _not_ in the metadata
    if titer_strain in hardcoded_map:
        curated_strain = hardcoded_map[titer_strain]
        found_in_metadata = curated_strain in metadata
        if not found_in_metadata:
            hardcoded_strains_not_found_in_metadata[(titer_strain, curated_strain)]+=1
        return Match(in_metadata=found_in_metadata, original_name=titer_strain, new_name=curated_strain)

    if titer_strain in metadata:
        return Match(in_metadata=True, original_name=titer_strain, new_name=titer_strain)

    # Check if the titer strain is in the hardcoded fauna map, and if so can we use the EPI_ISL
    # to connect this to the new metadata?
    if titer_strain in fauna_map:
        epi_isl = fauna_map[titer_strain]
        if epi_isl in epi_isl_to_strain: # comes from metadata, and thus curated_strain 
            in_metadata = True # epi_isl_to_strain is from the metadata
            curated_strain = epi_isl_to_strain[epi_isl]
            return Match(in_metadata=in_metadata,  original_name=titer_strain, new_name=curated_strain)

    # Update via hardcoded pattern lookup. One day we'll have a complete pipeline for this, but for
    # now we use a very thin layer to capture some obvious missing strains
    case_corrected_strain = titer_strain
    for case_correction in HARDCODED_CASE_CHANGES.items():
        if case_correction[0] in case_corrected_strain:
            case_corrected_strain = case_corrected_strain.replace(case_correction[0], case_correction[1])
    if case_corrected_strain != titer_strain and case_corrected_strain in metadata:
        return Match(in_metadata=True, original_name=titer_strain, new_name=case_corrected_strain,
            name_changed_via_pattern_match=True)

    # If we've reached this part then we haven't been able to confidently match strains.
    # See if it's possible to match using simplified names, but don't actually
    # change the name using this approach
    strain_simple = simplify_strain(titer_strain)
    if strain_simple in simplified_metadata:
        curated_strain = simplified_metadata[strain_simple]
        matching_if_simplified_strain_pairs[(titer_strain, curated_strain)]+=1
        return Match(in_metadata=False, original_name=titer_strain,
            new_name=titer_strain, new_name_matching_when_simplified=curated_strain)

    # no match :(
    return Match(in_metadata=False, original_name=titer_strain, new_name=titer_strain)


def get_year(match:Match, metadata:Metadata, collapse_before:int=1990)->str:
    """
    Find the year of the strain, either by looking it up in metadata (if provided)
    or examining the strain name itself. Unknown years will be *UNKNOWN_YEAR*.
    Years prior to *collapse_before* will be collapsed into f"<={collapse_before}"
    """

    if match.in_metadata:
        assert match.new_name in metadata
        metadata_year = metadata[match.new_name][1]
        if isinstance(metadata_year, int):
            if metadata_year <= collapse_before:
                return f"<={str(collapse_before)}"
            return str(metadata_year)
        # if the metadata year isn't an integer it means the metadata didn't include it
        # so fall through to strain-name parsing

    strain = cast(str, match.new_name_matching_when_simplified if match.new_name_matching_when_simplified else match.original_name)

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


@dataclasses.dataclass
class SuccessPerYear:
    total: int = 0
    metadata_match: int = 0
    potential_match: int = 0

def process_measurement_pairs(match_pairs: list[tuple[Match,Match]],
        metadata: Metadata,
        additional_stats_data: str,
        stats_fname: str
        ) -> None:
    """Count _unique_ strain names / measurement pairs seen in the data
    and report per-year stats about the number of matches / potential matches
    for downstream visualization.
    Note: uniqueness is based on the original titer strain name(s)
    """
    
    strain_matches_per_year = defaultdict(lambda: SuccessPerYear(total=0))
    measurement_matches_per_year = defaultdict(lambda: SuccessPerYear(total=0))
    strains_seen: set[str] = set()
    measurements_seen: set[tuple[str, str]] = set()

    for pair in match_pairs:
        virus_match, serum_match = pair
        
        if virus_match.original_name not in strains_seen:
            strains_seen.add(virus_match.original_name)
            virus_year = get_year(virus_match, metadata)
            strain_matches_per_year[virus_year].total +=1
            if virus_match.in_metadata:
                strain_matches_per_year[virus_year].metadata_match +=1
            elif virus_match.new_name_matching_when_simplified:
                strain_matches_per_year[virus_year].potential_match +=1
        
        if serum_match.original_name not in strains_seen:
            strains_seen.add(serum_match.original_name)
            serum_year = get_year(serum_match, metadata)
            strain_matches_per_year[serum_year].total +=1
            if serum_match.in_metadata:
                strain_matches_per_year[serum_year].metadata_match +=1
            elif serum_match.new_name_matching_when_simplified:
                strain_matches_per_year[serum_year].potential_match +=1

        # track the pair of matches (i.e. the measurement)
        # using the virus_strain as the source of date information
        if (virus_match.original_name, serum_match.original_name) not in measurements_seen:
            virus_year = get_year(virus_match, metadata)
            measurements_seen.add((virus_match.original_name, serum_match.original_name))
            measurement_matches_per_year[virus_year].total+=1
            if virus_match.in_metadata and serum_match.in_metadata:
                measurement_matches_per_year[virus_year].metadata_match+=1
            elif (virus_match.in_metadata or virus_match.new_name_matching_when_simplified) and \
                    (serum_match.in_metadata or serum_match.new_name_matching_when_simplified):
                measurement_matches_per_year[virus_year].potential_match+=1

    def collect_totals(data):
        """Collect all the per-year data and sum it up to report the totals"""
        total_points = sum([spy.total for spy in data.values()])
        total_metadata_match = sum([spy.metadata_match for spy in data.values()])
        total_potential_match = sum([spy.potential_match for spy in data.values()])
        return SuccessPerYear(total=total_points, metadata_match=total_metadata_match, potential_match=total_potential_match)
    
    strain_matches_per_year['total'] = collect_totals(strain_matches_per_year)
    measurement_matches_per_year['total'] = collect_totals(measurement_matches_per_year)

    print(f"Printing matching stats to {stats_fname!r}", file=sys.stderr)
    with open(stats_fname, 'w') as fh:
        print(compact_json({
            **json.loads(additional_stats_data),
            'strains': dict(sorted([k, dataclasses.asdict(v)] for k,v in strain_matches_per_year.items())),
            'measurements': dict(sorted([k, dataclasses.asdict(v)] for k,v in measurement_matches_per_year.items())),
        }), file=fh)

    def stats_str(a: int, b: int)->str:
        return f"{a:,}/{b:,} ({a/b*100:.1f}%)"

    print(f"Total matching (unique) strains:    {stats_str(strain_matches_per_year['total'].metadata_match, strain_matches_per_year['total'].total)}", file=sys.stderr)
    print(f"Total matching (unique) measurement pairs: {stats_str(measurement_matches_per_year['total'].metadata_match, measurement_matches_per_year['total'].total)}", file=sys.stderr)

    print(f"Potential matching (unique) strains:    {stats_str(strain_matches_per_year['total'].potential_match, strain_matches_per_year['total'].total)}", file=sys.stderr)
    print(f"Potential matching (unique) measurement pairs: {stats_str(measurement_matches_per_year['total'].potential_match, measurement_matches_per_year['total'].total)}", file=sys.stderr)


@dataclasses.dataclass
class StrainCounter:
    n_virus_strain: int = 0
    n_serum_strain: int = 0
    potential_match: str = ""
    year: str = ""

def print_missing_measurements(match_pairs: list[tuple[Match,Match]],
        metadata: Metadata,
        additional_stats_data: str,
        fname: str
        ) -> None:
    print(f"Printing unmatched strains to {fname!r}", file=sys.stderr)
    counts = defaultdict(lambda: StrainCounter(n_virus_strain=0, n_serum_strain=0))
    for pair in match_pairs:
        virus_match, serum_match = pair
        
        if not virus_match.in_metadata:
            counts[virus_match.original_name].n_virus_strain += 1
            counts[virus_match.original_name].year = get_year(virus_match, metadata)
            if virus_match.new_name_matching_when_simplified:
                counts[virus_match.original_name].potential_match = str(virus_match.new_name_matching_when_simplified)
    
        if not serum_match.in_metadata:
            counts[serum_match.original_name].n_serum_strain += 1
            counts[serum_match.original_name].year = get_year(serum_match, metadata)
            if serum_match.new_name_matching_when_simplified:
                counts[serum_match.original_name].potential_match = str(serum_match.new_name_matching_when_simplified)
    
    # hardcoded header information as we expect certain information in --stats-metadata
    header = ['titer_strain', 'year', 'virus_strain_count', 'serum_strain_count', 'potential_matching_strain',
        'lineages', 'centers', 'passages', 'assays']
    
    extra_data = json.loads(additional_stats_data)
    
    with open(fname, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(header)
        for titer_strain, counter in counts.items():
            writer.writerow([
                titer_strain, counter.year,
                counter.n_virus_strain, counter.n_serum_strain,
                counter.potential_match,
                extra_data.get('lineage', 'unknown'),
                extra_data.get('center', 'unknown'),
                extra_data.get('passage', 'unknown'),
                extra_data.get('assay', 'unknown'),
            ])
        
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--titers", required=True, metavar='TSV',
                        help="Titers TSV with column 'virus_strain' and 'serum_strain'")
    parser.add_argument("--fauna-strain-map", required=True, metavar='TSV',
                        help="TSV file which maps (fauna) strain name to EPI_ISL")
    parser.add_argument("--hardcoded-strain-map", required=True, metavar='TSV',
                        help="TSV file which maps (fauna) strain name (curated) strain name. Takes precedence over any other approach to correct names")
    parser.add_argument("--metadata", required=True, metavar='TSV',
                        help="Curated metadata TSV")
    parser.add_argument("--output", required=True, metavar='TSV',
                        help="The corrected --titers input TSV")
    parser.add_argument("--stats", required=True, metavar='JSON',
                        help="matching stats")
    parser.add_argument("--stats-metadata", required=True, metavar='JSON_STRING',
                        help="JSON-formatted data to export in the stats JSON")
    parser.add_argument("--missing-strains", required=False, metavar='TSV',
                        help="Write out information on unmatched strain names")

    args = parser.parse_args()

    _strain_matches: dict[str, tuple[str, str]] = {}
    _measurements_matches: dict[tuple[str,str], tuple[str, str]] = {}
    matches = {
        'strains': _strain_matches,
        'measurements': _measurements_matches,
    }

    hardcoded_map = read_hardcoded_strain_map(args.hardcoded_strain_map)
    fauna_map = read_fauna_map(args.fauna_strain_map)
    titers = read_titers(args.titers)
    metadata = parse_metadata(args.metadata)
    epi_isl_to_strain = {v[0]:k for k,v in metadata.items()}
    simplified_metadata = simplify_metadata(metadata)
    # Track pairs of strains (titer-strain, curated-strain) which would match if simplified
    matching_if_simplified_strain_pairs: defaultdict[tuple[str,str], int] = defaultdict(int)
    hardcoded_strains_not_found_in_metadata: defaultdict[tuple[str,str], int] = defaultdict(int)

    # Keep the match results from all the rows in the titers TSV to report stats later on
    measurement_matches: list[tuple[Match, Match]] = []

    for row in titers:
        virus_match = match(row['virus_strain'], fauna_map, hardcoded_map, metadata, epi_isl_to_strain, simplified_metadata, matching_if_simplified_strain_pairs, hardcoded_strains_not_found_in_metadata)
        serum_match = match(row['serum_strain'], fauna_map, hardcoded_map, metadata, epi_isl_to_strain, simplified_metadata, matching_if_simplified_strain_pairs, hardcoded_strains_not_found_in_metadata)
        measurement_matches.append((virus_match, serum_match))
        # update titers row in-place
        row['virus_strain'] = virus_match.new_name
        row['serum_strain'] = serum_match.new_name
        row['virus_strain_match'] = "yes" if virus_match.in_metadata else "maybe" if virus_match.new_name_matching_when_simplified else "no"
        row['serum_strain_match'] = "yes" if serum_match.in_metadata else "maybe" if serum_match.new_name_matching_when_simplified else "no"

    write_titers(args.output, titers)


    # Log out hardcoded titer strain name changes where the corrected name wasn't found in the metadata
    for pair, count in hardcoded_strains_not_found_in_metadata.items():
        print(f"[WARNING] Hardcoded strain map TSV changed {pair[0]!r} → {pair[1]!r}, but the latter is not found in the metadata. Strain appears {count:,} times", file=sys.stderr)

    # Log out titer strains which match curated strains when simplified, as we want to either
    # update the curated metadata to match the titers or remap the titer strain to match.
    for pair, count in matching_if_simplified_strain_pairs.items():
        print(f"[POTENTIAL MATCH] Titer strain {pair[0]!r} would match curated strain {pair[1]!r} if simplified. Strain appears {count:,} times", file=sys.stderr)

    # Log out how many strain names were changed due to application of pattern matching
    print(f"Changed {sum([m.name_changed_via_pattern_match for pair in measurement_matches for m in pair])} titer strain names via (hardcoded) pattern-matches", file=sys.stderr)

    process_measurement_pairs(measurement_matches, metadata, args.stats_metadata, args.stats)
    
    if args.missing_strains:
        print_missing_measurements(measurement_matches, metadata, args.stats_metadata, args.missing_strains)