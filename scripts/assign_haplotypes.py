#!/usr/bin/env python3
import argparse
from collections import defaultdict
from functools import partial
import sys

from augur.io import read_metadata
from augur.io.file import PANDAS_READ_CSV_OPTIONS
from augur.utils import write_json
import pandas as pd


def nucleotide_substitutions_match(substitutions_string, substitutions_list):
    """Returns True/False based on whether the given comma-delimited string of
    nucleotide substitutions matches all substitutions in the given list.

    substitutions_string looks like T291C,A566C,G615A,A1011G,A1257G
    substitutions_list looks like ['291C', '566C']
    """
    if len(substitutions_string) == 0 or len(substitutions_list) == 0:
        return True

    substitutions_string_list = [
        sub[1:]
        for sub in substitutions_string.split(",")
    ]
    return all(
        sub in substitutions_string_list
        for sub in substitutions_list
    )


def aa_substitutions_match(substitutions_string, substitutions_list):
    """Returns True/False based on whether the given comma-delimited string of
    amino acid substitutions matches all substitutions in the given list.

    substitutions_string looks like HA1:G78S,HA1:T135A,HA1:R150K,HA1:V223I
    substitutions_list looks like [('HA1', '122D'), ('HA1', '276E')]
    """
    if len(substitutions_string) == 0 or len(substitutions_list) == 0:
        return True

    substitutions_string_list = [
        (sub.split(":")[0], sub.split(":")[1][1:])
        for sub in substitutions_string.split(",")
    ]
    return all(
        sub in substitutions_string_list
        for sub in substitutions_list
    )


def assign_haplotype(record, haplotype_definitions, clade_column, default_haplotype):
    """Assign the most precise haplotype to the given record based on the given
    haplotype definitions.

    # subclade
    # substitutions
    # aaSubstitutions
    # founderMuts['subclade'].substitutions
    # founderMuts['subclade'].aaSubstitutions

    Example haplotype definitions with a single haplotype:

    {'J.2:135A': {'aa': [('HA1', '135A')], 'clade': 'J.2'}}
    """
    assigned_name = default_haplotype
    for name, definition in haplotype_definitions.items():
        if "clade" in definition:
            # Try to assign this haplotype based on its clade and clade-specific
            # substitutions.
            clade_match = (record[clade_column] == definition["clade"])
            nucleotide_column = f"founderMuts[\'{clade_column}\'].substitutions"
            aa_column = f"founderMuts[\'{clade_column}\'].aaSubstitutions"
        else:
            # Try to assign this haplotype based on all substitutions.
            clade_match = True
            nucleotide_column = "substitutions"
            aa_column = "aaSubstitutions"

        if (
            clade_match and
            nucleotide_substitutions_match(record[nucleotide_column], definition.get("nuc", [])) and
            aa_substitutions_match(record[aa_column], definition.get("aa", []))
        ):
            assigned_name = name

    return assigned_name


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--substitutions", required=True, help="TSV file with clades and substitutions from Nextclade")
    parser.add_argument("--haplotypes", required=True, help="TSV file of haplotype definitions in 'augur clades' format")
    parser.add_argument("--metadata-id-columns", default=["strain", "seqName"], help="names of possible columns in the substitutions table to use as the record id")
    parser.add_argument("--clade-column", default="subclade", help="name of the column in the substitutions table corresponding to clades used in the haplotype definitions")
    parser.add_argument("--haplotype-column-name", default="haplotype", help="name of the column or attribute to store the annotated haplotype in the output")
    parser.add_argument("--default-haplotype", default="unassigned", help="value to assign to records without any match to the given haplotypes")
    parser.add_argument("--output-table", required=True, help="TSV file of substitutions annotated by haplotype")
    parser.add_argument("--output-node-data", help="JSON in Nextstrain's node data format with haplotypes annotated per record id")

    args = parser.parse_args()

    substitutions = read_metadata(
        args.substitutions,
        id_columns=args.metadata_id_columns,
    )

    if args.haplotype_column_name in substitutions.columns:
        print(
            f"ERROR: The requested column name for haplotype annotations, '{args.haplotype_column_name}', already exists in the substitutions table.",
            file=sys.stderr,
        )
        sys.exit(1)

    haplotype_definitions = pd.read_csv(
        args.haplotypes,
        sep='\t' if args.haplotypes.endswith('.tsv') else ',',
        comment='#',
        na_filter=False,
        **PANDAS_READ_CSV_OPTIONS,
    )

    haplotype_definition_by_name = {}
    for haplotype_name, haplotype_definition in haplotype_definitions.groupby("clade", sort=False):
        definition = {}
        for record in haplotype_definition.to_dict(orient="records"):
            if record["gene"] == "clade":
                # When a haplotype has a clade definition, all substitutions in
                # the definition will be relative to that clade.
                if "clade" in definition:
                    print(
                        f"ERROR: The haplotype '{haplotype_name}' has multiple clades in its definition which is not possible.",
                        file=sys.stderr,
                    )
                    sys.exit(1)

                definition["clade"] = record["site"]
            elif record["gene"] == "nuc":
                # Nucleotide substitutions look like "A7G", so we only need to
                # compare each haplotype defining substitution to the position
                # and derived allele.
                if "nuc" not in definition:
                    definition["nuc"] = []

                definition["nuc"].append(record["site"] + record["alt"])
            else:
                # Amino acid substitutions look like "HA1:Q173P", so we need to
                # compare each haplotype defining substitution's gene and
                # position/derived allele.
                if "aa" not in definition:
                    definition["aa"] = []

                definition["aa"].append(
                    (
                        record["gene"],
                        record["site"] + record["alt"],
                    ),
                )

        haplotype_definition_by_name[haplotype_name] = definition

    assign_haplotype_per_record = partial(
        assign_haplotype,
        haplotype_definitions=haplotype_definition_by_name,
        clade_column=args.clade_column,
        default_haplotype=args.default_haplotype,
    )

    # Assign haplotypes to each row.
    substitutions[args.haplotype_column_name] = substitutions.apply(
        assign_haplotype_per_record,
        axis=1,
    )

    substitutions.to_csv(
        args.output_table,
        sep="\t",
        index=True,
        header=True,
    )

    if args.output_node_data:
        node_data = {
            strain: {args.haplotype_column_name: haplotype}
            for strain, haplotype in substitutions["haplotype"].to_dict().items()
        }
        write_json({"nodes": node_data}, args.output_node_data)
