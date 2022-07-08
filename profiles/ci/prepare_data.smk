rule prepare_sequences:
    input:
        sequences="example_data/h3n2_{segment}.fasta",
    output:
        sequences="data/h3n2/raw_{segment}.fasta",
    shell:
        """
        cp -f {input.sequences} {output.sequences}
        """

rule prepare_titers:
    input:
        titers="example_data/cdc_h3n2_cell_fra_titers.tsv",
    output:
        titers="data/h3n2/cdc_cell_fra_titers.tsv",
    shell:
        """
        cp -f {input.titers} {output.titers}
        """
