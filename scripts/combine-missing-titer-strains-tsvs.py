#!/usr/bin/env python3

"""
Combine multiple missing-titer-strains TSV files into a single merged TSV.

Merging logic:
- Unique key: titer_strain (first column)
- Count columns (virus_strain_count, serum_strain_count) are summed
- potential_matching_strain values must be identical across files for a given
  titer_strain (an error is raised if they conflict)
- Categorical columns (lineages, centers, passages, assays) collect unique
  values across files and are written as sorted, comma-separated sets
"""

import argparse
import csv
import sys

ERROR_STRING = "ERROR"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True, metavar="TSV",
                        help="Missing-titer-strains TSV files to merge")
    parser.add_argument("--output", required=True, metavar="TSV",
                        help="Output merged TSV")
    args = parser.parse_args()

    count_cols = ['virus_strain_count', 'serum_strain_count']
    set_cols = ['lineages', 'centers', 'passages', 'assays']

    merged: dict[str, dict] = {}

    for fname in args.input:
        with open(fname) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                strain = row['titer_strain']
                # Initialise empty row each time we see a new strain
                if strain not in merged:
                    merged[strain] = {
                        'titer_strain': strain,
                        'year': row['year'],
                        **{c: 0 for c in count_cols},
                        'potential_matching_strain': row['potential_matching_strain'],
                        **{c: set() for c in set_cols},
                    }
                
                entry = merged[strain]
                
                for c in count_cols:
                    entry[c] += int(row[c])
                
                # Handle "potential_matching_strain" and "year" which we expect
                # to be the same for all occurances of a given strain name. It's
                # conceivable these are different if they are from different
                # lineages and thus are matching against different metadata
                # TSVs. If this ever happens we don't throw an error (because it
                # would bring down the entire pipeline) but rather print a
                # warning and use the ERROR_STRING as the TSV value
                for key in ['year', 'potential_matching_strain']:
                    existing_value = entry[key]
                    new_value = row[key]
                    if existing_value == ERROR_STRING:
                        continue
                    if not existing_value and not new_value:
                        continue
                    if new_value != existing_value:
                        print(f"ERROR: conflicting {key!r} for titer strain {strain!r}:"
                            f" {existing_value!r} ({','.join(entry['lineages'])} {','.join(entry['centers'])} {','.join(entry['passages'])} {','.join(entry['assays'])})"
                            f" vs {new_value!r} ({row['lineages']} {row['centers']} {row['passages']} {row['assays']})",
                            file=sys.stderr)
                        entry[key] = ERROR_STRING

                for c in set_cols:
                    for val in row[c].split(','):
                        val = val.strip()
                        if val:
                            entry[c].add(val)

    def sorted_year(year_str:str) -> int:
        """return an integer representation for sorting purposes"""
        try:
            return int(year_str)
        except ValueError:
            if year_str[0:2] == '<=':
                return 2
            if year_str=='unknown':
                return 1
            return 0

    header = ['titer_strain', 'year', *count_cols, 'potential_matching_strain', *set_cols]
    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        for entry in sorted(merged.values(),
            key=lambda e: (sorted_year(e['year']), sum(e[c] for c in count_cols)),
            reverse=True
        ):
            writer.writerow([
                entry['titer_strain'],
                entry['year'],
                *[entry[c] for c in count_cols],
                entry['potential_matching_strain'],
                *[','.join(sorted(entry[c])) for c in set_cols],
            ])

    print(f"Merged {len(args.input)} files → {len(merged):,} unique strains → {args.output}", file=sys.stderr)
