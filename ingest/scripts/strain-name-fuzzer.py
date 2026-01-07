#!/usr/bin/env python3
"""
Provided with two lists of strain names we want to match all query strains
against the list of curated strains using fuzzy matching.

"""
import argparse
import difflib
import sys

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
        lines = [line for line in fh]
    return lines

def normalize_for_matching(s: str) -> str:
    """Normalize string: remove hyphens and convert to lowercase."""
    return s.replace('-', '').lower()


def extract_digits(s: str) -> str:
    """Extract all digits from a string in order."""
    return ''.join(c for c in s if c.isdigit())

def no_zeros(s: str) -> str:
    return ''.join(c for c in s if c!='0')


def calculate_similarity(query_norm: str, target_norm: str) -> float:
    """
    Calculate similarity score that:
    - Ignores hyphens and case differences
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

def find_highest(similarities: list[tuple[str, float]]) -> tuple[str, float]:
    h:float = 0
    s = ""
    for (strain, score) in similarities:
        if score>h:
            s = strain
            h = score
    return (s, h)



def extract_strain(line: str, is_tsv: bool) -> str|None:
    if line.startswith('#') or line.strip()=='':
        return None
    if not is_tsv:
        return line.strip()
    # assume it's the first column, and ignore any subtlety in splitting TSVs for now
    return line.split('\t')[0].strip()

def replace_strain(new_strain:str, old_strain:str, line: str) -> str:
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
    parser.add_argument("--query-strains", required=True, help="Text file with query strain names (e.g. include.txt)")
    parser.add_argument("--strain-map", required=False, nargs='*', help="Hardcoded map of query strain (fauna) to new curated strain")
    parser.add_argument("--no-fuzz", required=False,action='store_true', help="Skip fuzzing for efficiency reasons")

    args = parser.parse_args()

    direct_matches: int = 0
    normalized_match: int = 0
    fuzzy_matches: int = 0
    no_matches: int = 0
    query_strains: int = 0
    strain_map_hit: int = 0

    curated_strains = read_strains(args.curated_strains)
    curated_strains_set = set(curated_strains)
    curated_strains_normalized = [normalize_for_matching(w) for w in curated_strains]
    curated_strains_normalized_set = set(curated_strains_normalized)
    query_lines = read_text_file(args.query_strains)
    strain_map = parse_strain_map(args.strain_map)
    is_tsv = args.query_strains.endswith(".tsv")

    # STDOUT is the new strains text file
    print(f"# Strains list generated via fuzzy matching")
    print(f"# Original file: {args.query_strains!r}")
    print(f"")


    for idx, line in enumerate(query_lines):
        if idx%100==0:
            print(f"{idx:,}/{len(query_lines):,}", file=sys.stderr)

        # strain can be thought of as fauna strain
        strain = extract_strain(line, is_tsv)

        if strain is None or (strain=='strain' and is_tsv):
            print(line, end="")
            continue

        query_strains+=1

        if strain in strain_map:
            print(replace_strain(strain_map[strain], strain, line), end="")
            strain_map_hit+=1
            continue

        if strain in curated_strains_set: # direct match - just echo back the original line
            print(line, end="")
            direct_matches += 1
            continue
        
        query_norm = normalize_for_matching(strain)

        # Check if the normalized string exists, and call it a normalized match
        if query_norm in curated_strains_normalized_set:
            print(replace_strain(query_norm, strain, line), end="")
            normalized_match+=1
        
        if args.no_fuzz:
            print()
            print("# TODO XXX Couldn't find a matching strain for the following line:")
            print(line, end="")
            no_matches+=1
            continue


        # Loop over the curated strains to calculate similarities
        similarities = [similarity for similarity in 
            [(curated, calculate_similarity(query_norm, curated_norm))
            for (curated, curated_norm) in zip(curated_strains, curated_strains_normalized)]
            if similarity[1] >= 0.8 # TODO XXX
        ]
        if len(similarities)==0:
            print()
            print("# TODO XXX Couldn't find a matching strain for the following line:")
            print(line, end="")
            no_matches += 1
            continue

        best = find_highest(similarities)

        if best[1]==1: 
            print(replace_strain(best[0], strain, line), end="")
            normalized_match+=1
        elif best[1]>0.9:
            print(f"# Following strain was close match ({best[1]:.3f}) to original strain {strain!r}")
            print(replace_strain(best[0], strain, line), end="")
            fuzzy_matches += 1
        else:
            print()
            print("# TODO XXX Couldn't find a matching strain for the following line:")
            print(line, end="")
            no_matches += 1

    if query_strains:
        print_everywhere(f"\n# === Summary {args.query_strains} ===")
        print_everywhere(f"# Direct matches:         {direct_matches:,} / {query_strains:,} ({100*direct_matches/query_strains:.1f}%)")
        print_everywhere(f"# Changed via lookup:     {strain_map_hit:,} / {query_strains:,} ({100*strain_map_hit/query_strains:.1f}%)")
        print_everywhere(f"# Normalized matches:     {normalized_match:,} / {query_strains:,} ({100*normalized_match/query_strains:.1f}%)")
        print_everywhere(f"# Fuzzy matches:          {fuzzy_matches:,} / {query_strains:,} ({100*fuzzy_matches/query_strains:.1f}%)")
        print_everywhere(f"# No matches:             {no_matches:,} / {query_strains:,} ({100*no_matches/query_strains:.1f}%)")
        

    