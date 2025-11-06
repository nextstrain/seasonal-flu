"""
This part of the workflow handles fetching the metadata and sequences from
open data sources.

OUTPUTS:

    metadata      = data/{lineage}/metadata.tsv
    sequences     = data/{lineage}/{segment}.fasta

"""

# Fetch metadata from GenSpectrum
# Fetch sequences from GenSpectrum
# Pull out GenBank accessions from GenSpectrum metadata
# Fetch from Entrez
# Merge and collapse segment Entrez metadata with GenSpectrum metadata
# Produce 1 metadata TSV and 8 segment FASTA
