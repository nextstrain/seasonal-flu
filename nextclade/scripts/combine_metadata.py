import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine metadata with nextclade output")
    parser.add_argument("--nextclade", required=True, help="Path to nextclade output tsv file")
    parser.add_argument("--metadata-raw", required=True, help="Path to metadata tsv file")
    parser.add_argument("--nextclade-ha", required=True, help="Path to nextclade tsv file")
    parser.add_argument("--nextclade-na", required=True, help="Path to nextclade tsv file")
    parser.add_argument("--nextclade-columns", nargs='+', required=True, help="Columns to select from nextclade output")
    parser.add_argument("--lineage", required=True, help="Lineage name")
    parser.add_argument("--segment", required=True, help="Segment name")
    parser.add_argument("--output", required=True, help="Path to output tsv file")
    args = parser.parse_args()

    # read nextclade output and metadata
    nextclade_df = pd.read_csv(args.nextclade, sep='\t', usecols=args.nextclade_columns, index_col='seqName')
    metadata_df = pd.read_csv(args.metadata_raw, sep='\t')


    # merge on strain column
    combined_df = pd.merge(metadata_df, nextclade_df, left_on='strain', right_index=True, how='left')

    if args.lineage in ['h3n2', 'h1n1pdm', 'vic', 'b']:
        if args.segment != 'ha':
            nextclade_ha_df = pd.read_csv(args.nextclade_ha, sep='\t', index_col='seqName')
            combined_df = pd.merge(combined_df, nextclade_ha_df[['clade']], left_on='strain', right_index=True, how='left', suffixes=('', '_ha'))

        if args.segment != 'na':
            nextclade_na_df = pd.read_csv(args.nextclade_na, sep='\t', index_col='seqName')
            combined_df = pd.merge(combined_df, nextclade_na_df[['clade']], left_on='strain', right_index=True, how='left', suffixes=('', '_na'))
    # write combined metadata to output file
    combined_df.to_csv(args.output, sep='\t', index=False)
