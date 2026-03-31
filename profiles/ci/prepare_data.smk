rule prepare_sequences:
    input:
        sequences="example_data/h3n2_{segment}.fasta",
    output:
        sequences="data/h3n2/{segment}.fasta",
    shell:
        """
        cp -f {input.sequences} {output.sequences}
        """


rule prepare_metadata:
    input:
        metadata="example_data/h3n2_metadata.tsv",
    output:
        metadata="data/h3n2/metadata.tsv",
    shell:
        """
        cp -f {input.metadata} {output.metadata}
        """
