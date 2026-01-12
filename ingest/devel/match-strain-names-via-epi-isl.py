#!/usr/bin/env python3
"""
Provides a fauna strain name -> ndjson (curated) strain name mapping
by matching records on EPI ISL.
"""
import argparse
import sys
import io
import zstandard as zstd
from augur.io import read_metadata
from augur.io.json import load_ndjson

type EpiIsl = str
type StrainName = str

def open_file(file_path):
    """
    Open a file, decompressing with zstd if it has a .zst extension.
    Returns a text-mode file handle.
    """
    if file_path.endswith('.zst'):
        dctx = zstd.ZstdDecompressor()
        fh = open(file_path, 'rb')
        reader = dctx.stream_reader(fh)
        # Wrap in TextIOWrapper to get text mode
        return io.TextIOWrapper(reader, encoding='utf-8')
    else:
        return open(file_path, 'r')

def parse_existing_metadata(fname: str) -> dict[EpiIsl, StrainName]:
    df = read_metadata(fname)
    return dict(zip(df['gisaid_epi_isl'], df.index))

def parse_ndjson(fname: str) -> dict[EpiIsl, StrainName]:
    with open_file(fname) as fh:
        records = load_ndjson(fh)
        return {record['gisaid_epi_isl']:record['strain'] for record in records}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fauna", required=True, metavar='TSV',
                        help="Original (fauna) metadata TSV with 'strain' and 'gisaid_epi_isl' columns")
    parser.add_argument("--curated", required=True, metavar='NDJSON',
                        help="Curated NDJSON")
    parser.add_argument("--output", required=True, metavar='TSV',
                        help="TSV of strain names which have changed. Map is FAUNA STRAIN -> NEWLY CURATED STRAIN")
    args = parser.parse_args()

    fauna = parse_existing_metadata(args.fauna)
    curated = parse_ndjson(args.curated)
    changed_fh = open(args.output, 'w')
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

    print(f"Using {args.fauna} as source of truth...", file=sys.stderr)
    print(f"{missing=:,}")
    print(f"{changed=:,}")
    print(f"{unchanged=:,}")
    changed_fh.close()
