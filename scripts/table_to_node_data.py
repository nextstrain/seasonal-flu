"""Create Augur-compatible node data JSON from a pandas data frame.
"""
import argparse
import pandas as pd
from augur.utils import write_json


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", help="table to convert to a node data JSON")
    parser.add_argument("--index-column", default="strain", help="name of the column to use as an index")
    parser.add_argument("--delimiter", default=",", help="separator between columns in the given table")
    parser.add_argument("--node-name", default="nodes", help="name of the node data attribute in the JSON output")
    parser.add_argument("--output", help="node data JSON file")

    args = parser.parse_args()

    if args.output is not None:
        table = pd.read_csv(
            args.table,
            sep=args.delimiter,
            index_col=args.index_column,
            dtype=str,
        )

        # # Convert columns that aren't strain names or labels to floats.
        # for column in table.columns:
        #     if column != "strain" and not "label" in column:
        #         table[column] = table[column].astype(float)

        table_dict = table.transpose().to_dict()
        write_json({args.node_name: table_dict}, args.output)
