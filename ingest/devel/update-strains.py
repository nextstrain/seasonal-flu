#!/usr/bin/env python3

"""
Given a list of hardcoded strain names (e.g. `include.txt`)
Use a number of approaches to update the names to the match the corresponding strain
from the newly curated metadata.
Prints the updated strain list file to STDOUT
"""

import argparse
import difflib
import heapq
import sys
from typing import TypedDict
import json

def read_strains(fname:str)->list[str]:
    strains = list()
    with open(fname, 'r', newline='') as fh:
        for line in fh:
            contents = line.strip()
            if not contents or contents.startswith('#'):
                continue
            strains.append(contents)
    return strains


def read_text_file(fname:str)->list[str]:
    with open(fname, 'r', newline='') as fh:
        # If it's a JSON then we don't even bother to parse the structure and make you do it
        # manually. Sorry.
        if fname.endswith('.json'):
            j = json.load(fh)
            lines = [f"{name}\n" for name in j['nodes'].keys()]
        else:
            lines = [line for line in fh]
    return lines

def normalize_for_matching(s: str) -> str:
    """Normalize string: remove hyphens and convert to lowercase."""
    return s.replace('-', '').replace(' ', '').lower()


def extract_digits(s: str) -> str:
    """Extract all digits from a string in order."""
    return ''.join(c for c in s if c.isdigit())

def no_zeros(s: str) -> str:
    return ''.join(c for c in s if c!='0')


def calculate_similarity(query_norm: str, target_norm: str) -> float:
    """
    Calculate similarity score that:
    - Ignores hyphens and case differences (because inputs are normalized)
    - Heavily penalizes digit mismatches
    """
    # passage mismatch is a non-starter
    if (query_norm.endswith("egg") and not target_norm.endswith("egg")) or (target_norm.endswith("egg") and not query_norm.endswith("egg")):
         return 0.0

    # Extract and compare digits
    query_digits = extract_digits(query_norm)
    target_digits = extract_digits(target_norm)

    # numbers cannot be different _EXCEPT FOR_ dropped zeros
    if no_zeros(query_digits)!=no_zeros(target_digits):
        return 0.0

    # Get base similarity on normalized strings
    base_similarity = difflib.SequenceMatcher(None, query_norm, target_norm).ratio()

    final_score = max(0.0, base_similarity)

    return final_score

class ScoredStrain(TypedDict):
    strain: str
    score: float

def find_highest(similarities: list[tuple[str, float]], n: int) -> list[ScoredStrain]:
    """Return the top n highest scoring matches as a list of dicts with 'strain' and 'score' keys."""
    top_n = heapq.nlargest(n, similarities, key=lambda x: x[1])
    return [{'strain': strain, 'score': score} for strain, score in top_n]

def extract_strain(line: str, is_tsv: bool) -> str|None:
    if line.startswith('#') or line.strip()=='':
        return None
    if not is_tsv:
        # Allow #-comments inline
        return line.split('#')[0].strip()
    # assume it's the first column, and ignore any subtlety in splitting TSVs for now
    return line.split('\t')[0].strip()

def replace_strain(old_strain:str, new_strain:str, line: str) -> str:
    assert old_strain in line, f"old strain not in line! {old_strain!r} {new_strain!r}"
    # assert new_strain not in line, f"new strain already in line! {old_strain!r} {new_strain!r}"
    return line.replace(old_strain, new_strain)

def print_everywhere(*args) -> None:
    print(*args, file=sys.stderr)
    print(*args)

def parse_strain_map(fname:str|None) -> dict[str,str]:
    if not fname:
        return {}
    with open(fname[0]) as fh:
        return dict([line.strip().split() for line in fh])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curated-strains", required=True, help="Text file with curated strain names")
    parser.add_argument("--fauna-strains", required=False, help="Text file with fauna strain names")
    parser.add_argument("--query-strains", required=True, help="Text file with query strain names (e.g. include.txt)")
    parser.add_argument("--strain-map", required=False, nargs='*', help="Hardcoded map of old to new strain names, where known")

    args = parser.parse_args()

    stats: dict[str, int] = {
        'direct_matches': 0,
        'normalized_match': 0,
        'fuzzy_matches': 0,
        'no_matches_and_missing_from_fauna': 0,
        'no_matches': 0,
        'query_strains': 0,
        'strain_map_hit': 0,
    }

    curated_strains = read_strains(args.curated_strains)
    curated_strains_set = set(curated_strains)

    fauna_strains = read_strains(args.fauna_strains) if args.fauna_strains else []

    curated_strains_normalized = {normalize_for_matching(w):w for w in curated_strains}
    # curated_strains_normalized_set = set(curated_strains_normalized.keys())
    query_lines = read_text_file(args.query_strains)
    strain_map = parse_strain_map(args.strain_map)
    is_tsv = args.query_strains.endswith(".tsv")

    for idx, line in enumerate(query_lines):
        if idx and idx%100==0:
            print(f"[update-strains] [progress] {idx:,}/{len(query_lines):,}", file=sys.stderr)

        # strain can be thought of as fauna strain
        strain = extract_strain(line, is_tsv)

        # preserve comments, whitespace etc
        if strain is None or (strain=='strain' and is_tsv):
            print(line, end="")
            continue

        stats['query_strains']+=1

        # Obviously direct hits just get echo'd straight out ;)
        if strain in curated_strains_set:
            print(line, end="")
            stats['direct_matches'] += 1
            continue

        # If it's in the hardcoded strain map then it's similarly trivial
        if new_strain:=strain_map.get(strain, ''):
            print(replace_strain(strain, new_strain, line), end="")
            stats['strain_map_hit']+=1
            continue

        # Now things get a little more tricky...
        # Start by normalizing the string
        query_norm = normalize_for_matching(strain)
        # and checking if it was even in fauna!
        exists_in_fauna = strain in fauna_strains

        # Check if the normalized string exists, and call it a normalized match
        if new_strain:=curated_strains_normalized.get(query_norm, ''):
            print(replace_strain(strain, new_strain, line), end="")
            stats['normalized_match']+=1

        # Loop over the curated strains to calculate similarities
        similarities = [similarity for similarity in 
            [(curated, calculate_similarity(query_norm, curated_norm))
            for curated_norm, curated in curated_strains_normalized.items()]
            if similarity[1] >= 0.8 # TODO XXX
        ]

        if len(similarities):
            stats['fuzzy_matches'] += 1
            best = find_highest(similarities, 3)
            print(f"# [update-strains] BEST {len(best)} FUZZY MATCHES and their scores. Original strain {strain!r}. You must manually choose the correct match.")
            if len(fauna_strains) and not exists_in_fauna:
                print(f"# [update-strains] P.S. Original strain wasn't in fauna!")
            for candidate in best:
                print(f"# [update-strains] [score={candidate['score']:.3f}] " + replace_strain(strain, candidate['strain'], line), end="")
            print(line, end="")
            continue
        
        # Now is the time to admit we failed!
        # ...but wait - was it even in fauna?
        if len(fauna_strains):
            if strain not in fauna_strains:
                print("# [update-strains] NOTE: Following strain unable to be matched, but it also wasn't in fauna to start with:")
                print(line, end="")
                stats['no_matches_and_missing_from_fauna']+=1
            else:
                print("# [update-strains] No match for the following strain, which did exist in fauna:")
                print(line, end="")
                stats['no_matches']+=1
        else:
            print("# [update-strains] No match for the following strain")
            print(line, end="")
            stats['no_matches']+=1

    def _stats_str(a:str, b:str)->str:
        aa = stats[a]
        bb = stats[b]
        return f"{aa:,} / {bb:,} ({100*aa/bb:.1f}%)"

    if stats['query_strains']:
        print_everywhere(f"\n# [update-strains] === Summary {args.query_strains} ===")
        print_everywhere(f"# [update-strains] Direct matches:           {_stats_str('direct_matches', 'query_strains')}")
        print_everywhere(f"# [update-strains] Changed via lookup:       {_stats_str('strain_map_hit', 'query_strains')}")
        print_everywhere(f"# [update-strains] Normalized matches:       {_stats_str('normalized_match', 'query_strains')}")
        print_everywhere(f"# [update-strains] Fuzzy matches:            {_stats_str('fuzzy_matches', 'query_strains')}")
        print_everywhere(f"# [update-strains] Unmatched & not in fauna: {_stats_str('no_matches_and_missing_from_fauna', 'query_strains')}")
        print_everywhere(f"# [update-strains] No matches:               {_stats_str('no_matches', 'query_strains')}")