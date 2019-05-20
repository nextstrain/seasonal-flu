import argparse
import pandas as pd

def select(file, mergeby, fields):
    entries = pd.read_csv(file, sep='\t')
    mapping = {}
    for index, row in entries.iterrows():
        if row[mergeby]:
            values = []
            for field in fields:
                if row[field]:
                    values.append(row[field])
                else:
                    values.append("?")
            mapping[row[mergeby]] = values
    for key, values in mapping.items():
        print(key + "\t" + "\t".join(values))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Concatenate multiple tsvs, merging specified columns",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--files', nargs='+', type=str, required=True, help="list of tsv files")
    parser.add_argument('--mergeby', type=str, help="column name to merge on")
    parser.add_argument('--fields', nargs='+', type=str, help="column names to include")
    args = parser.parse_args()

    print(args.mergeby + "\t" + "\t".join(args.fields))
    for file in args.files:
        select(file, args.mergeby, args.fields)
