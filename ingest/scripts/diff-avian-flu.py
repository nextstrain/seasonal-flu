#!/usr/bin/env python3

"""
This script is specifically for diffing "avian-flu" NDJSONs [1] against the avian-flu TSV
constructed from Fauna. The aim is to minimise differences, or explain them.

The intention is not for this to be used in day-to-day curation once we are no longer using
Fauna. Instead we'll use `scripts/diff-ndjson` similar to how we do for seasonal-flu.

[1] What an "avian-flu" NDJSON actually means is a little up for grabs, as we're subsampling
by subtypes vs fauna's historical adding of strains via filtering rules which varied over time
"""



import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from typing import Any


type TsvRecord = dict[str, Any]
type NdjsonRecord = dict[str, Any]


def hardcoded_diff_skip(key:str, truth_record: TsvRecord, query_record: NdjsonRecord)->bool:
    if key=='isolate_id':
        return True # skip fauna TSV index field (we used this for matching purposes)

    # Don't report the always-present but not in the NDJSON 'virus="avian_flu"'
    # TODO XXX - do we want to categorise certain subtypes as this instead?
    if key=='virus' and truth_record.get('virus')=='avian_flu':
        return True

    # If the fauna value is '?' and the query is empty string, just ignore this
    if truth_record.get(key, '')=='?' and query_record.get(key, "missing")=='':
        return True
    
    return False

punctuation = ['/','_', '(', ')', ',', '-', ' ', ';']
def remove_punctuation(v:str) -> str:
    vv=v[:]
    for char in punctuation:
        vv = vv.replace(char, '')
    return vv

NULL={'NA', '', '?', 'tbd', 'n/a'}

def is_genbank(accession:str) -> bool:
    genbank_regex = re.compile(r"[A-Z]{2}\d+")
    genbank_regex_versioned = re.compile(r"[A-Z]{2}\d+\.\d+")
    return bool(genbank_regex.fullmatch(accession) or genbank_regex_versioned.fullmatch(accession))


def understand_strain_differences(data: defaultdict[tuple[str,str], list[tuple[TsvRecord, NdjsonRecord]]]) -> None:
    counts = {
        'explained by dashes in new value': 0,
        'explained by removing punctuation': 0,
        'ndjson has trailing /YYYY': 0,
        'fauna has trailing HxNx': 0,
        'differ by field 1 only (type)': 0,
        'differ by field 2 only (host)': 0,
        'differ by field 3 only (location)': 0,
        'differ by field 4 only (id)': 0,
        'differ by field 5 only (year)': 0,
        'unexplained': 0,
    }
    n=0

    def exclude_part(field_num:int, v:str) -> str:
        parts = v.split('/')
        return "/".join([*parts[0:field_num-1], *parts[field_num:]])
    
    for el, records in data.items():
        occurances = len(records)
        n+=occurances
        [fauna_val, ndjson_val] = [el[0].lower(), el[1].lower()]

        ndjson_val = ndjson_val.replace('-', '')
        if fauna_val==ndjson_val:
            counts['explained by dashes in new value'] += occurances
            continue

        if remove_punctuation(fauna_val)==remove_punctuation(ndjson_val):
            counts['explained by removing punctuation'] += occurances
            continue
        
        if groups:=re.match(r'^(.+)/\d{4}$', ndjson_val):
            if remove_punctuation(fauna_val)==remove_punctuation(groups.group(1)):
                counts['ndjson has trailing /YYYY'] += occurances
                continue

        if groups:=re.match(r'^(.+?)_?h\dn\d$', fauna_val):
            if remove_punctuation(groups.group(1))==remove_punctuation(ndjson_val):
                counts['fauna has trailing HxNx'] += occurances
                continue

        if exclude_part(1, fauna_val)==exclude_part(1, ndjson_val):
            counts['differ by field 1 only (type)'] += occurances
            continue
        
        if exclude_part(2, fauna_val)==exclude_part(2, ndjson_val):
            counts['differ by field 2 only (host)'] += occurances
            continue

        if exclude_part(3, fauna_val)==exclude_part(3, ndjson_val):
            counts['differ by field 3 only (location)'] += occurances
            continue

        if exclude_part(4, fauna_val)==exclude_part(4, ndjson_val):
            counts['differ by field 4 only (id)'] += occurances
            continue

        if exclude_part(5, fauna_val)==exclude_part(5, ndjson_val):
            counts['differ by field 5 only (year)'] += occurances
            continue

        counts['unexplained'] +=occurances

    for reason,count in counts.items():
        print(f"\t\t{reason:<35}{count:<7,} ({count*100/n:.1f}% of conflicts)")


def understand_lab_differences(data: defaultdict[tuple[str,str], list[tuple[TsvRecord, NdjsonRecord]]]) -> None:
    counts = {
        'fauna_was_null': 0,
        'punctuation difference': 0,
        'unexplained': 0,
    }
    n=0
    for el, records in data.items():
        occurances = len(records)
        n+=occurances
        [fauna_val, ndjson_val] = [el[0].lower(), el[1].lower()]
        if fauna_val.lower() in NULL:
            counts['fauna_was_null']+=occurances
            continue
        elif remove_punctuation(fauna_val)==remove_punctuation(ndjson_val):
            counts['punctuation difference']+=occurances
            continue
        counts['unexplained']+=occurances

    for reason,count in counts.items():
        if count:
            print(f"\t\t{reason:<35}{count:<7,} ({count*100/n:.1f}% of conflicts)")

def understand_host_differences(data: defaultdict[tuple[str,str], list[tuple[TsvRecord, NdjsonRecord]]]) -> None:
    known_good = {'human', 'avian', 'environment', 'cattle', 'equine', 'nonhuman mammal', 'other'}
    counts = {
        'fauna_was_null': 0,
        'ndjson_is_now_null': 0,
        'both valid, but different!': 0,
        'fauna was too specific': 0,
        'ndjson is too specific!': 0,
        'unexplained': 0,
    }
    n=0
    for el, records in data.items():
        occurances = len(records)
        n+=occurances
        [fauna_val, ndjson_val] = [el[0].lower(), el[1].lower()]
        if fauna_val in NULL:
            counts['fauna_was_null']+=occurances
            continue
        if ndjson_val in NULL:
            counts['ndjson_is_now_null']+=occurances
            continue
        if ndjson_val in known_good and fauna_val in known_good:
            counts['both valid, but different!']+=occurances
            continue
        if ndjson_val in known_good:
            counts['fauna was too specific']+=occurances
            continue
        if fauna_val in known_good:
            counts['ndjson is too specific!']+=occurances
            continue
        counts['unexplained']+=occurances

    for reason,count in counts.items():
        if count:
            print(f"\t\t{reason:<35}{count:<7,} ({count*100/n:.1f}% of conflicts)")


def understand_location_differences(key: str, data: defaultdict[tuple[str,str], list[tuple[TsvRecord, NdjsonRecord]]]) -> None:
    counts = {
        'fauna_was_null': 0,
        'ndjson_is_now_null': 0,
        'explained by fauna bugs': 0,
        'unexplained': 0,
    }
    n=0
    for el, records in data.items():
        occurances = len(records)
        n+=occurances
        [fauna_val, ndjson_val] = [el[0].lower(), el[1].lower()]
        if fauna_val in NULL:
            counts['fauna_was_null']+=occurances
            continue
        if ndjson_val in NULL:
            counts['ndjson_is_now_null']+=occurances
            continue

        # see 2 bugs in fauna https://github.com/nextstrain/fauna/blob/c38b1b44cb119973793438529769ccb17250419d/vdb/avian_flu_upload.py#L612-L616
        # just look at the first Fauna record, assuming the country is the same...
        if (key=='division' or key=='location') and fauna_val==records[0][0]['country'].lower():
            counts['explained by fauna bugs']+=occurances
            continue
        counts['unexplained']+=occurances

    for reason,count in counts.items():
        if count:
            print(f"\t\t{reason:<35}{count:<7,} ({count*100/n:.1f}% of conflicts)")

def understand_subtype_differences(data: defaultdict[tuple[str,str], list[tuple[TsvRecord, NdjsonRecord]]]) -> None:
    counts = {
        'fauna_was_null': 0,
        'unexplained': 0,
    }
    n=0
    for el, records in data.items():
        occurances = len(records)
        n+=occurances
        [fauna_val, ndjson_val] = [el[0].lower(), el[1].lower()]
        if fauna_val in NULL:
            counts['fauna_was_null']+=occurances
            continue
        counts['unexplained']+=occurances

    for reason,count in counts.items():
        if count:
            print(f"\t\t{reason:<35}{count:<7,} ({count*100/n:.1f}% of conflicts)")

def understand_date_differences(data: defaultdict[tuple[str,str], list[tuple[TsvRecord, NdjsonRecord]]]) -> None:
    counts = {
        'genbank accessions': 0,
        'fauna date pre 1950': 0,
    }
    n=0
    for el, records in data.items():
        occurances = len(records)
        n+=occurances
        counts['genbank accessions']+=sum([is_genbank(r[0]['isolate_id']) for r in records])
        
        [fauna_val, ndjson_val] = [el[0], el[1]]
        if int(fauna_val.split('-')[0])<1950:
            counts['fauna date pre 1950']+=occurances

    for reason,count in counts.items():
        if count:
            print(f"\t\t{reason:<35}{count:<7,} ({count*100/n:.1f}% of conflicts)")


def understand_generic_differences(data: defaultdict[tuple[str,str], list[tuple[TsvRecord, NdjsonRecord]]]) -> None:
    counts = {
        'fauna_was_null': 0,
        'ndjson_is_now_null': 0,
        'unexplained': 0,
    }
    n=0
    for el, records in data.items():
        occurances = len(records)
        n+=occurances
        [fauna_val, ndjson_val] = [el[0].lower(), el[1].lower()]
        if fauna_val in NULL:
            counts['fauna_was_null']+=occurances
            continue
        if ndjson_val in NULL:
            counts['ndjson_is_now_null']+=occurances
            continue
        counts['unexplained']+=occurances

    for reason,count in counts.items():
        if count:
            print(f"\t\t{reason:<35}{count:<7,} ({count*100/n:.1f}% of conflicts)")




def diff(truth_record: TsvRecord, query_record: NdjsonRecord) -> list[tuple[str, str, Any, Any]]:
    """
    Compare two dictionaries and return a list of differences.
    Each difference is a tuple of (change_type, key, truth_value, query_value).

    Allow the query record to have extra keys and don't report these
    (as the NDJSON has heaps of fields not in the TSV)
    """
    differences = []

    for key, truth_value in truth_record.items():
        if hardcoded_diff_skip(key, truth_record, query_record):
            continue

        if key not in query_record:
            differences.append(("key_removed", key, truth_value, None))
            continue
        
        query_value = query_record[key]

        if truth_value==query_value:
            continue
        
        if truth_value.lower()==query_value.lower():
            differences.append(("conflict_case_only", key, truth_value, query_value))
        else:
            differences.append(("conflict", key, truth_value, query_value))

    return differences


def load_ndjson(filepath):
    """Load NDJSON file and return a dictionary keyed by strain."""
    records = {}

    #   "segment": "HA",
    #   "accession": "EPI289873",

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            try:
                record = json.loads(line.strip())
                sequences = record['sequences']
                try:
                    key = [s['accession'] for s in sequences if s['segment']=='HA'][0] # isolate_id
                except IndexError:
                    print(f"Skipping strain {record.get('strain')} as no HA sequence")
                    continue
                # key = record.get('strain')
                # if key is None:
                    # raise Exception(f"Warning: Line {line_num} in {filepath} missing 'strain' key")
                records[key] = record
            except json.JSONDecodeError as e:
                raise Exception("JSON ERROR!", e)
    return records


def load_tsv(filepath):
    """Load TSV file and return a dictionary keyed by 'isolate_id' (HA)."""
    records = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for line_num, row in enumerate(reader, 2):  # Start at 2 to account for header
            key = row.get('isolate_id')
            if key is None:
                raise Exception(f"Warning: Line {line_num} in {filepath} missing 'isolate_id' column")
            # Convert empty strings to None for consistency
            record = {k: (None if v == '' else v) for k, v in row.items()}
            records[key] = record
    return records

def format_value(val):
    """Format a value for human-readable output."""
    if val is None:
        return "<missing>"
    if isinstance(val, str):
        return f'"{val}"'
    if isinstance(val, (dict, list)):
        # For complex types, show a truncated JSON representation
        s = json.dumps(val)
        if len(s) > 80:
            return s[:77] + "..."
        return s
    return str(val)



def remap_via_simplified_strain(missing_keys: set[str], added_keys: set[str], truth_records: dict[str, TsvRecord], query_records: dict[str, NdjsonRecord]) -> None:
    def simplify(name:str, records: dict[str, Any])->dict[str,str]:
        lookup:defaultdict[str,list[tuple[str, dict[str,Any]]]] = defaultdict(list)
        for accession, record in records.items():
            strain = record['strain'].lower().replace("/",'').replace(' ','').replace('.','')
            lookup[strain].append((accession, record))
        lookup_no_dups: dict[str,str] = {k:v[0][0] for k,v in lookup.items() if len(v)==1}
        print(f"\t{name} n={len(lookup_no_dups):,} unambiguous strain names, n={len(lookup)-len(lookup_no_dups):,} (simplified) strain names were duplicates")
        return lookup_no_dups
    
    truths = simplify("TRUTH", truth_records)
    queries = simplify("QUERY", query_records)
    count = 0

    for simplified_strain_name, truth_accession in truths.items():
        if simplified_strain_name not in queries:
            continue
        query_accession = queries[simplified_strain_name]
        if truth_accession in missing_keys:
            if query_accession not in added_keys:
                print(f"[WARNING] Query accession {query_accession} matches an accession in the truth, but via strain matching matches another record (accession: {truth_accession}). We will use the accession match.")
                continue
            count+=1
            query_record = query_records[query_accession]
            query_records[truth_accession] = query_record
            del query_records[query_accession]
        
    print(f"\tremapped {count:,} by comparing (unambiguous & simplified) strain names")


def remap_via_epi_epi(missing_keys: set[str], added_keys: set[str], truth_records: dict[str, TsvRecord], query_records: dict[str, NdjsonRecord]) -> None:
    # Sometimes the truth (TSV) accession has a duplicated "EPI" prefix, e.g. EPIEPI2743546
    count = 0
    for truth_accession in missing_keys:
        if truth_accession.startswith("EPIEPI"):
            if (query_accession:=truth_accession.replace("EPIEPI", "EPI")) in query_records:
                count+=1
                query_record = query_records[query_accession]
                query_records[truth_accession] = query_record
                del query_records[query_accession]
                continue
        
    print(f"\tremapped {count:,} via EPIEPI -> EPI comparison")

def compare_records(truth_records, query_records):
    truth_keys = set(truth_records.keys())
    query_keys = set(query_records.keys())
    missing_keys = truth_keys - query_keys
    added_keys = query_keys - truth_keys
    common_keys = truth_keys & query_keys

    looks_like_genbank = set()
    for accession in missing_keys:
        if is_genbank(accession):
            looks_like_genbank.add(accession)

    print(f"\t{len(missing_keys):,} MISSING ROWS (in truth, not in query)")
    print(f"\t\tNOTE: {len(looks_like_genbank)} of these look like GenBank accessions")
    print(f"\t{len(added_keys):,} ADDED ROWS (in query, not in truth)")
    print(f"\t{len(common_keys):,} COMMON ROWS (based on HA isolate id)")


    return (truth_keys, query_keys, missing_keys, added_keys, common_keys)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True, metavar="FILE", help="Source of truth file (TSV)")
    parser.add_argument("--query", required=True, metavar="FILE", help="Query file (NDJSON )")

    args = parser.parse_args()

    # Detect and load files
    truth_records: dict[str, dict[str, Any]] = load_tsv(args.truth) # key=isolate_id
    query_records: dict[str, dict[str, Any]] = load_ndjson(args.query) # key=HA accession

    # Find missing and added rows based on the HA accession to start with
    print("\n\nNDJSON - TSV linkage via HA isolate id / (HA) isolate_id")
    (truth_keys, query_keys, missing_keys, added_keys, common_keys) = compare_records(truth_records, query_records)

    print("\nConsidering duplicated EPIEPI accessions in the truth TSV...")
    remap_via_epi_epi(missing_keys, added_keys, truth_records, query_records)
    (truth_keys, query_keys, missing_keys, added_keys, common_keys) = compare_records(truth_records, query_records)

    print("\nConsidering strain name matches (lowecase, no spaces comparisons)")
    remap_via_simplified_strain(missing_keys, added_keys, truth_records, query_records)
    (truth_keys, query_keys, missing_keys, added_keys, common_keys) = compare_records(truth_records, query_records)


    print()
    # Allow missing "truth" keys IF --extra-truth, otherwise report the rows
    missing_keys_by_type: defaultdict[str, int] = defaultdict(int)
    if missing_keys:
        found_differences = True
        print(f"{len(missing_keys):,} MISSING ROWS (in truth, not in query):")
        for isolate_id in sorted(missing_keys):
            type_info = f"{truth_records[isolate_id]['subtype']}"
            host = truth_records[isolate_id]['host']
            strain = truth_records[isolate_id]['strain']
            print(f"  - {isolate_id:<12} {type_info:<8} {host:<10} ({strain})")
            missing_keys_by_type[type_info]+=1
        print()

    # Report added rows with some helpful info about the record
    added_keys_by_type: defaultdict[str, int] = defaultdict(int)
    if added_keys:
        found_differences = True
        print(f"{len(added_keys):,} ADDED ROWS (in query, not in truth):")
        for isl in sorted(added_keys): # not actually ISL - it's the HA accession!!!!
            type_info = f"{query_records[isl]['vtype']}/{query_records[isl]['subtype']}/{query_records[isl]['lineage']}"
            host = query_records[isl]['host']
            strain = query_records[isl]['strain']
            print(f"  - {isl:<15} {type_info:<12} {host:<12} ({strain})")
            added_keys_by_type[type_info]+=1
        print()

    n_common_rows_which_are_different = 0
    total_differences_added_key = 0.0
    total_differences_removed_key = 0.0
    total_differences_value_casing = 0.0
    total_differences_value_contents = 0.0
    differences_value_casing:defaultdict[str,int]=defaultdict(int)

    # Store the value differences (for rows we can match) so that we can run heuristics
    # on them later on with the intention of understanding differences
    value_differences: defaultdict[str, defaultdict[tuple[str,str], list[tuple[TsvRecord, NdjsonRecord]]]] = defaultdict(lambda: defaultdict(list))


    # Compare common rows
    for key in sorted(common_keys):
        truth_record = truth_records[key]
        query_record = query_records[key]

        differences = diff(truth_record, query_record)

        if differences:
            n_common_rows_which_are_different+=1
            print(f"DIFFERENCES FOR {key}\t{query_record['gisaid_epi_isl']}\t{query_record['gisaid_strain']}")
            for conflict_type, key, truth_val, query_val in differences:
                if conflict_type=='key_removed':
                    total_differences_removed_key+=1
                    print(f"  Key removed: {key} (was: {format_value(truth_val)})")
                elif conflict_type=='conflict_case_only':
                    total_differences_value_casing+=1
                    differences_value_casing[key]+=1
                    print(f"  Value differs BY CASE ONLY at {key}: {format_value(truth_val)} → {format_value(query_val)}")
                elif conflict_type=='conflict':
                    total_differences_value_contents+=1
                    value_differences[key][(truth_val, query_val)].append((truth_record, query_record))
                    print(f"  Value differs at {key}: {format_value(truth_val)} → {format_value(query_val)}")
                else:
                    raise Exception(f"Unknown conflict type {conflict_type=}")


            print()

    print(f"\n\n------ SUMMARY -----")
    print(f"n(truth records, TSV) = {len(truth_records.keys()):,}")
    print(f"n(query records, NDJSON) = {len(query_records.keys()):,}")
    print(f"\n{len(missing_keys):,} MISSING ROWS (in truth, not in query)")
    if (len(missing_keys)):
        for key,value in sorted(list(missing_keys_by_type.items()), key=lambda el: el[1], reverse=True):
            print(f"\t{key:<15}{value:,}")
    print(f"\n{len(added_keys):,} ADDED ROWS (in query, not in truth)")
    if (len(added_keys)):
        for key,value in sorted(list(added_keys_by_type.items()), key=lambda el: el[1], reverse=True):
            print(f"\t{key:<15}{value:,}")
    print(f"\n{len(common_keys):,} COMMON ROWS (where we were able to link up records)")
    if n_common_rows_which_are_different:
        total_differences_added_key/=n_common_rows_which_are_different
        total_differences_removed_key/=n_common_rows_which_are_different
        total_differences_value_casing/=n_common_rows_which_are_different
        total_differences_value_contents/=n_common_rows_which_are_different
        print(f"\t{n_common_rows_which_are_different:,} of these have disagreements. On average:")
        print(f"\t{total_differences_added_key:.2f} added keys (i.e. new in query)")
        print(f"\t{total_differences_removed_key:.2f} removed keys")
        print(f"\t{total_differences_value_casing:.2f} value differences due to case only")
        print(f"\t{total_differences_value_contents:.2f} real value differences")
    # now try to understand value changes

    print("\n\nAmong records in common, values which differ by case only:")
    for key,value in sorted(list(differences_value_casing.items()), key=lambda el: el[1], reverse=True):
        print(f"\t{key:<18}{value:<7,} ({value*100/n_common_rows_which_are_different:.1f}% of records)")


    print("\n\nAmong records in common, values which differ by more than case:")
    for key,data in sorted(list(value_differences.items()), key=lambda el: sum([len(x) for x in el[1].values()]), reverse=True):    
        count = sum([len(x) for x in data.values()])
        print(f"\t{key:<18}{count:<7,} ({count*100/n_common_rows_which_are_different:.1f}% of records)")
        if key=='strain':
            understand_strain_differences(data)
        elif key in ['authors', 'originating_lab', 'submitting_lab']:
            understand_lab_differences(data)
        elif key in ['region', 'country', 'division', 'location']:
            understand_location_differences(key, data)
        elif key == 'subtype':
            understand_subtype_differences(data)
        elif key == 'host':
            understand_host_differences(data)
        elif key == 'date':
            understand_date_differences(data)
        else:
            understand_generic_differences(data)