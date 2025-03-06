SEGMENT_MAP = {
    "ha": "seg4",
    "na": "seg6",
    "mp": "seg7",
    "ns": "seg8",
    "np": "seg5",
    "pa": "seg3",
    "pb1": "seg2",
    "pb2": "seg1"
}

COLUMN_MAPPING = {
    "accessionVersion": "strain",
    "sampleCollectionDate": "date",
    "country": "region",
    "ncbiReleaseDate": "date_submitted"
}

def rename_columns(input_file, output_file, mapping=COLUMN_MAPPING):
    with open(input_file, "r") as f:
        header = f.readline().strip().split("\t")
        header = [mapping.get(h, h) for h in header]
        with open(output_file, "w") as g:
            g.write("\t".join(header) + "\n")
            for line in f:
                g.write(line)


rule download_metadata:
    output:
        metadata="data/{lineage}/raw_metadata_{segment}.tsv",
    conda:
        "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        curl -s "https://api.loculus.genspectrum.org/influenza-a/sample/details?downloadAsFile=true&versionStatus=LATEST_VERSION&isRevocation=false&dataFormat=tsv&subtypeHA=H3&subtypeNA=N2" > {output.metadata}
        """


rule download_sequences:
    output:
        sequences="data/{lineage}/{segment}.fasta",
    conda:
        "../../workflow/envs/nextstrain.yaml"
    params:
        segment_number=lambda wildcards: SEGMENT_MAP[wildcards.segment],
    shell:
        """
        curl -s "https://api.loculus.genspectrum.org/influenza-a/sample/unalignedNucleotideSequences/{params.segment_number}?downloadAsFile=true&versionStatus=LATEST_VERSION&isRevocation=false&dataFormat=fasta&subtypeHA=H3&subtypeNA=N2" > {output.sequences}
        """
    
rule rename_columns:
    input:
        metadata="data/{lineage}/raw_metadata_{segment}.tsv",
    output:
        renamed_metadata="data/{lineage}/metadata_{segment}.tsv",
    params:
        mapping=COLUMN_MAPPING,
    run:
        rename_columns(
            input.metadata, output.renamed_metadata, mapping=params.mapping
        )