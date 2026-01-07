#!/usr/bin/env python3
"""
Provides a fauna strain name -> ndjson (curated) strain name mapping
by matching records on EPI ISL.
"""
import argparse
import sys
from augur.io import read_metadata
from augur.io.json import load_ndjson

type EpiIsl = str
type StrainName = str

def parse_existing_metadata(fname: str) -> dict[EpiIsl, StrainName]:
    df = read_metadata(fname)
    return dict(zip(df['gisaid_epi_isl'], df.index))

def parse_ndjson(fname: str) -> dict[EpiIsl, StrainName]:
    with open(fname) as fh:
        records = load_ndjson(fh)
        return {record['gisaid_epi_isl']:record['strain'] for record in records}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-metadata", required=True, help="Original (fauna) metadata TSV with 'strain' and 'gisaid_epi_isl' columns")
    parser.add_argument("--new-metadata", required=True, help="Curated NDJSON")
    parser.add_argument("--changed", required=True, help="TSV of strain names which have changed. Map is FAUNA -> NEW")
    args = parser.parse_args()

    fauna = parse_existing_metadata(args.original_metadata)
    curated = parse_ndjson(args.new_metadata)
    changed_fh = open(args.changed, 'w')
    print("FAUNA_STRAIN\tCURATED_STRAIN", file=changed_fh)

    [unchanged, changed, missing] = [0,0,0]
    for epi_isl, fauna_name in fauna.items():
        new_name = curated.get(epi_isl, None)
        if new_name is None:
            missing+=1
            continue
        if fauna_name == new_name:
            unchanged+=1
        else:
            changed+=1
            print(f"{fauna_name}\t{new_name}", file=changed_fh)

    print(f"Using {args.original_metadata} as source of truth...", file=sys.stderr)
    print(f"{missing=:,}")
    print(f"{changed=:,}")
    print(f"{unchanged=:,}")
    changed_fh.close()