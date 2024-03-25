ruleorder: download_parsed_sequences > parse
ruleorder: download_parsed_metadata > annotate_metadata_with_reference_strains

rule download_sequences:
    output:
        sequences="data/{lineage}/raw_{segment}.fasta"
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{segment}/raw_sequences.fasta.xz"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - | xz -c -d > {output.sequences}
        """

rule download_titers:
    output:
        titers="data/{lineage}/{titer_collection}_titers.tsv"
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{titer_collection}_titers.tsv.gz"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - | gzip -c -d > {output.titers}
        """

rule download_parsed_sequences:
    output:
        sequences="data/{lineage}/{segment}.fasta"
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{segment}/sequences.fasta.xz"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - | xz -c -d > {output.sequences}
        """

rule download_parsed_metadata:
    output:
        metadata="data/{lineage}/metadata.tsv",
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/metadata.tsv.xz"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - | xz -c -d > {output.metadata}
        """
