#!/usr/bin/env python3
"""
Write a GFF3 annotation whose coding regions are CDS rows.

TreeTime reads the coding regions from CDS rows with a strand and a phase. The genemaps
of this workflow describe them as gene rows in the style of augur, with the name in
`gene_name`. Each gene row becomes a CDS row with phase 0 and the name as ID and Name;
an annotation that has CDS rows already is copied unchanged.
"""
import argparse
import re


def gene_name(attributes):
    for key in ('Name', 'gene_name', 'gene', 'ID'):
        match = re.search(rf'(?:^|;)\s*{key}\s*=\s*"?([^";]+)"?', attributes)
        if match:
            return match.group(1)
    raise ValueError(f"No name in the attributes '{attributes}'")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--annotation', required=True, help='GFF3 annotation')
    parser.add_argument('--output', required=True, help='GFF3 annotation with CDS rows')
    args = parser.parse_args()

    lines = open(args.annotation).read().splitlines()
    rows = [line.split('\t') for line in lines if line and not line.startswith('#')]
    has_cds = any(len(row) > 2 and row[2] == 'CDS' for row in rows)

    with open(args.output, 'w') as fh:
        for line in lines:
            fields = line.split('\t')
            if not has_cds and len(fields) == 9 and fields[2] == 'gene':
                name = gene_name(fields[8])
                fields[2], fields[7] = 'CDS', '0'
                fields[8] = f"ID={name};Name={name};gene_name={name}"
                line = '\t'.join(fields)
            fh.write(line + '\n')
