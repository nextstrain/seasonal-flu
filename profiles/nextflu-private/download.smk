from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

rule download_sequences:
    input:
        sequences=S3.remote("nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{segment}/raw_sequences.fasta.xz")
    output:
        sequences="data/{lineage}/raw_{segment}.fasta"
    shell:
        """
        xz -c -d {input.sequences} > {output.sequences}
        """

rule download_titers:
    input:
        titers=S3.remote("nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{titer_collection}_titers.tsv.gz")
    output:
        titers="data/{lineage}/{titer_collection}_titers.tsv"
    shell:
        """
        gzip -c -d {input.titers} > {output.titers}
        """
