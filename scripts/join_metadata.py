"""Merge two tables.
"""
import argparse
import sys

from augur.io import read_metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", nargs="+", required=True, help="metadata files to join")
    parser.add_argument("--segments", nargs="+", required=True, help="names of segments corresponding to each metadata file")
    parser.add_argument("--segment-columns", default=["accession"], nargs="+", help="segment-specific columns of the metadata to copy from each segment's file")
    parser.add_argument("--how", default="outer", choices=["inner", "outer"], help="how to join metadata files")
    parser.add_argument("--output", required=True, help="joined metadata")

    args = parser.parse_args()

    if len(args.metadata) != len(args.segments):
        print(
            "ERROR: You must provide as many segment names as metadata files.",
            file=sys.stderr
        )
        sys.exit(1)

    metadata = None
    for segment, segment_metadata_file in zip(args.segments, args.metadata):
        segment_metadata = read_metadata(segment_metadata_file)
        segment_metadata[segment] = True
        segment_columns = {
            column: f'{column}_{segment}'
            for column in args.segment_columns
        }

        if metadata is None:
            metadata = segment_metadata.rename(columns=segment_columns)
        else:
            # Augur's `read_metadata` function indexes by the first valid strain
            # id it finds, so we do not need to specify a key to join on.
            metadata = metadata.join(
                segment_metadata.loc[:, [segment] + args.segment_columns].rename(columns=segment_columns),
                how=args.how,
            )
            metadata[segment] = metadata[segment].fillna(False)

    metadata.to_csv(
        args.output,
        sep='\t',
        na_rep="N/A",
    )
