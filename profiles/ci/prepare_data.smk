rule prepare_sequences:
    input:
        sequences="example_data/h3n2_{segment}.fasta",
    output:
        sequences="data/h3n2/raw_{segment}.fasta",
    shell:
        """
        cp -f {input.sequences} {output.sequences}
        """

rule prepare_nextclade:
    input:
        nextclade="example_data/nextclade_{segment}.tsv.xz",
    output:
        nextclade="data/h3n2/{segment}/nextclade.tsv.xz",
    shell:
        """
        cp -f {input.nextclade} {output.nextclade}
        """
