"""
This part of the workflow handles the curation of the open data

REQUIRED INPUTS:

    metadata      = data/{lineage}/metadata.tsv
    sequences     = data/{lineage}/{segment}.fasta

OUTPUTS:

    metadata      = results/{lineage}/metadata.tsv
    sequences     = results/{lineage}/{segment}.fasta

"""

# Curate metadata
# Deduplicate by strain
# Filter FASTAs by deduped metadata
# Replace seq ids in FASTAs with strain

