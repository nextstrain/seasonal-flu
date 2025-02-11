"""
Links the records in the downloaded files from GISAID EpiFlu.
Excel metadata records are linked with segment sequences in the FASTA and
output as NDJSON to stdout.

Each record represents a single GISAID record, formatted as:

    {
        “Isolate_Id”: “...”,
        “Isolate_Name”: “...”,
        “Collection_Date”: “...”,
        “Passage_History”: “...”,
        […other metadata fields…],
        “sequences”: [
            {
                “segment”: “HA”,
                “accession”: “...”,
                “sequence”: “...”
            },
            {
                “segment”: “NA”,
                “accession”: “...”,
                “sequence”: “...”
            },
            […other 8 segments…]
        ]
    }
"""
import argparse
import os
import pyfastx
import re
from textwrap import dedent
from typing import Iterable
from augur.io.json import dump_ndjson
from augur.io.metadata import DEFAULT_DELIMITERS, read_table_to_dict
from augur.io.print import print_err


# Expected segments from GISAID EpiFlu
# The last two segments (HE and P3) are unused, but keeping for completion
#   -Jover, 11 February 2025
SEGMENTS = [
    "PB2",
    "PB1",
    "PA",
    "HA",
    "NP",
    "NA",
    "MP",
    "NS",
    "HE",
    "P3",
]

SEGMENT_ACCESSION_PATTERN = r'^(?P<accession>EPI\d+)\|'
DEFAULT_UNKNOWN_VALUE = "?"

# Expected columns in the GISAID xls file
# Hard-coding here as I don't expect them to change, but if they get updated
# often enough then we can consider making them CLI options.
#   -Jover, 11 February 2025
RECORD_ID_COLUMN = "Isolate_Id"
SEGMENT_COLUMNS = {
    segment: f"{segment} Segment_Id"
    for segment in SEGMENTS
}


def link_metadata_and_sequences(metadata: Iterable[dict],
                                sequences: pyfastx.Fasta) -> Iterable[dict]:
    """
    Link records in the provided *metadata* with segment sequences in the
    provided *sequences*. Drops the segment fields and adds a "sequences"
    field, which is an array of all of the segments.
    """
    for record in metadata:
        linked_record = record.copy()
        record_id = linked_record[RECORD_ID_COLUMN]
        record_seqs = []
        unmatched_accessions = {}
        unmatched_sequences = {}
        for segment, segment_column in SEGMENT_COLUMNS.items():
            segment_id = linked_record.pop(segment_column)
            accession = parse_segment_accession(segment_id)
            sequence = get_segment_sequence(sequences, record_id,
                                            segment, accession)

            if segment_id and accession == DEFAULT_UNKNOWN_VALUE:
                unmatched_accessions[segment] = segment_id

            if accession != DEFAULT_UNKNOWN_VALUE and sequence == DEFAULT_UNKNOWN_VALUE:
                unmatched_sequences[segment] = accession

            record_seqs.append({
                "segment": segment,
                "accession": accession,
                "sequence": sequence
            })

        if len(unmatched_accessions):
            print_err(
                "WARNING: Could not match segment accessions for record",
                f"{record_id!r} for the following segments: {unmatched_accessions!r}"
            )

        if len(unmatched_sequences):
            print_err(
                f"WARNING: Could not find sequences for record {record_id!r}",
                f"for the following segment accessions: {unmatched_sequences}"
            )

        linked_record["sequences"] = record_seqs
        yield linked_record


def parse_segment_accession(segment_id: str) -> str:
    """
    Parses the segment accession from the provided *segment_id*. If unable to
    parse the accession, then returns the `DEFAULT_UNKNOWN_VALUE`.
    """
    accession = DEFAULT_UNKNOWN_VALUE

    matches = re.search(SEGMENT_ACCESSION_PATTERN, segment_id)
    if matches is not None:
        accession = matches["accession"]

    return accession


def get_segment_sequence(sequences: pyfastx.Fasta, record_id: str,
                         segment: str, accession: str) -> str:
    """
    Gets the sequence matching the provided *accession* from the indexed
    *sequences. If there is not matching sequence, returns the
    `DEFAULT_UNKNOWN_VALUE`.
    """
    sequence = DEFAULT_UNKNOWN_VALUE

    # try/except was consistently faster than checking accession is in sequences
    #   -Jover, 11 February 2025
    try:
        sequence_record = sequences[accession]
    except KeyError:
        pass
    else:
        sequence = str(sequence_record.seq).upper()

    return sequence


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--metadata", metavar="<xls>",
        help=dedent(f"""\
            GISAID EpiFlu metadata xls file, which is expected to have the
            record id column {RECORD_ID_COLUMN!r} and the segment ID columns
            {[SEGMENT_COLUMNS.values()]!r}. Each segment ID column is expected
            to contain segment IDs with the segment accession matching
            {SEGMENT_ACCESSION_PATTERN}.
            """))

    parser.add_argument("--sequences", metavar="<fasta>",
        help=dedent(f"""\
            GISAID EpiFlu FASTA file, where the headers should only be the
            sequence accession. This is the “DNA Accession no.” field in the
            GISAID "Sequences (DNA) as FASTA" download options.
            """))

    args = parser.parse_args()

    metadata = read_table_to_dict(
        table=args.metadata,
        delimiters=DEFAULT_DELIMITERS,
        id_column=RECORD_ID_COLUMN
    )

    # Remove the old Pyfastx index to force rebuild of index
    # so we don't have to worry about a stale cached index
    #   -Jover, 11 February 2025
    try:
        os.remove(f"{args.sequences}.fxi")
    except FileNotFoundError:
        pass

    sequences = pyfastx.Fasta(args.sequences)
    linked_records = link_metadata_and_sequences(metadata, sequences)
    dump_ndjson(linked_records)
