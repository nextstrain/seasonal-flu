rule prepare_sequences:
    input:
        sequences="example_data/h3n2_{segment}_dms.fasta",
    output:
        sequences="data/h3n2/raw_{segment}.fasta",
    shell:
        """
        cp -f {input.sequences} {output.sequences}
        """
