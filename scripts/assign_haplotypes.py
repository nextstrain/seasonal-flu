#!/usr/bin/env python3
import argparse
from collections import defaultdict
from functools import partial
import sys

from augur.io import read_metadata, write_json
from augur.io.file import PANDAS_READ_CSV_OPTIONS
import pandas as pd


def nucleotide_substitutions_match(record_substitutions, required_substitutions):
    """Returns True/False based on whether the given comma-delimited string of
    nucleotide substitutions matches all substitutions in the given list.

    record_substitutions looks like T291C,A566C,G615A,A1011G,A1257G
    required_substitutions looks like ['291C', '566C']

    When required substitutions don't exist in the record's substitution list,
    there can't be a match.

    >>> nucleotide_substitutions_match("", ['291C', '566C'])
    False

    When there aren't any required substitutions, we call this a match.

    >>> nucleotide_substitutions_match("T291C,A566C,G615A,A1011G,A1257G", [])
    True

    When both record and required substitutions exist, all required
    substitutions must be in the record substitutions list.

    >>> nucleotide_substitutions_match("T291C,A566C,G615A,A1011G,A1257G", ['291C', '566C'])
    True
    >>> nucleotide_substitutions_match("T291C,A566C,G615A,A1011G,A1257G", ['291C', '566C', '123C'])
    False

    """
    if len(record_substitutions) > 0:
        record_substitutions = [
            sub[1:]
            for sub in record_substitutions.split(",")
        ]
    else:
        record_substitutions = []

    return all(
        sub in record_substitutions
        for sub in required_substitutions
    )


def aa_substitutions_match(record_substitutions, required_substitutions):
    """Returns True/False based on whether the given comma-delimited string of
    amino acid substitutions matches all substitutions in the given list.

    >>> aa_substitutions_match("", [('HA1', '122D'), ('HA1', '276E')])
    False
    >>> aa_substitutions_match("HA1:G78S,HA1:T135A,HA1:R150K,HA1:V223I", [])
    True
    >>> aa_substitutions_match("HA1:G78S,HA1:T135A,HA1:R150K,HA1:V223I", [('HA1', '122D'), ('HA1', '276E')])
    False
    >>> aa_substitutions_match("HA1:G78S,HA1:T135A,HA1:R150K,HA1:V223I", [('HA1', '135A'), ('HA1', '223I')])
    True

    """
    if len(record_substitutions) > 0:
        record_substitutions = [
            (sub.split(":")[0], sub.split(":")[1][1:])
            for sub in record_substitutions.split(",")
        ]
    else:
        record_substitutions = []

    return all(
        sub in record_substitutions
        for sub in required_substitutions
    )


def assign_haplotype(record, haplotype_definitions, clade_column, default_haplotype, use_clade_as_default_haplotype=False):
    """Assign the most precise haplotype to the given record based on the given
    haplotype definitions.

    >>> haplotype_definitions = {'J.2:135A': {'aa': [('HA1', '135A')], 'clade': 'J.2'}, 'J.3': {'clade': 'J.3'}, 'J.1': {'nuc': ['135C'], 'aa': [('HA1', '123R')]}}
    >>> assign_haplotype({'subclade': 'J.2', "founderMuts['subclade'].aaSubstitutions": "HA1:A135T"}, haplotype_definitions, "subclade", "unassigned")
    'unassigned'
    >>> assign_haplotype({'subclade': 'J.2', "founderMuts['subclade'].aaSubstitutions": "HA1:A135T"}, haplotype_definitions, "subclade", "unassigned", use_clade_as_default_haplotype=True)
    'J.2'
    >>> assign_haplotype({'subclade': 'J.2', "founderMuts['subclade'].aaSubstitutions": "HA1:T135A"}, haplotype_definitions, "subclade", "unassigned")
    'J.2:135A'
    >>> assign_haplotype({'subclade': 'J', 'aaSubstitutions': 'HA1:K123R,HA1:T135A', 'substitutions': 'A100T,T110C,G135C,T200G'}, haplotype_definitions, "subclade", "unassigned")
    'J.1'
    >>> assign_haplotype({'subclade': 'J.3'}, haplotype_definitions, "subclade", "unassigned")
    'J.3'

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

    # Allow unassigned records to default to their original clade annotation
    # instead of a hardcoded default value.
    if assigned_name == default_haplotype and use_clade_as_default_haplotype:
        assigned_name = record[clade_column]

    return assigned_name


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--substitutions", required=True, help="TSV file with clades and substitutions from Nextclade")
    parser.add_argument("--haplotypes", required=True, help="""
    TSV file of haplotype definitions in 'augur clades' format except with the 'clade' column name replaced with 'haplotype'.
    Haplotypes will be assigned to each input record in the order they appear in this file.
    Records matching multiple haplotypes will receive the haplotype that appears latest in the file.
    Define haplotypes that derive from existing clades by specifying 'clade' in the 'gene' field and the clade name in the 'site' field.
    All defining substitutions for derived haplotypes will be checked against the column specified with the `--clade-column` argument and corresponding 'founderMuts' column of the Nextclade annotations.
    For example, if a haplotype is defined relative the 'subclade' column, its amino acid substitutions will be checked against the "founderMuts['subclade'].aaSubstitutions" column.
    """)
    parser.add_argument("--metadata-id-columns", default=["strain", "seqName"], help="names of possible columns in the substitutions table to use as the record id")
    parser.add_argument("--clade-column", default="subclade", help="name of the column in the substitutions table corresponding to clades used in the haplotype definitions")
    parser.add_argument("--haplotype-column-name", default="haplotype", help="name of the column or attribute to store the annotated haplotype in the output")
    parser.add_argument("--default-haplotype", default="unassigned", help="value to assign to records without any match to the given haplotypes")
    parser.add_argument("--use-clade-as-default-haplotype", action="store_true", help="use the existing clade annotation for records without assigned haplotypes instead of using the hardcoded default value")
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

    if "haplotype" not in haplotype_definitions.columns:
        print(
            f"ERROR: The column 'haplotype' is missing from the given haplotype definitions file, '{args.haplotypes}'.",
            file=sys.stderr,
        )
        sys.exit(1)

    haplotype_definition_by_name = {}
    for haplotype_name, haplotype_definition in haplotype_definitions.groupby("haplotype", sort=False):
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

                definition["nuc"].append(str(record["site"]) + record["alt"])
            else:
                # Amino acid substitutions look like "HA1:Q173P", so we need to
                # compare each haplotype defining substitution's gene and
                # position/derived allele.
                if "aa" not in definition:
                    definition["aa"] = []

                definition["aa"].append(
                    (
                        record["gene"],
                        str(record["site"]) + record["alt"],
                    ),
                )

        haplotype_definition_by_name[haplotype_name] = definition

    assign_haplotype_per_record = partial(
        assign_haplotype,
        haplotype_definitions=haplotype_definition_by_name,
        clade_column=args.clade_column,
        default_haplotype=args.default_haplotype,
        use_clade_as_default_haplotype=args.use_clade_as_default_haplotype,
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
            for strain, haplotype in substitutions[args.haplotype_column_name].to_dict().items()
        }
        write_json({"nodes": node_data}, args.output_node_data)
