ruleorder: download_parsed_sequences > parse
ruleorder: download_parsed_metadata > annotate_metadata_with_reference_strains
ruleorder: download_nextclade > run_nextclade


DEFAULT_S3_PATHS = {
    "raw_sequences": "s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{segment}/raw_sequences.fasta.xz",
    "titers": "s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{titer_collection}_titers.tsv.gz",
    "sequences": "s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{segment}/sequences.fasta.xz",
    "metadata": "s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/metadata.tsv.xz",
}


def _get_s3_path(file_type: str) -> str:
    """
    Returns `config[{files_type}_url]` if it's provided,
    otherwise use the DEFAULT_S3_PATHS.

    Will raise Exception if the provided *file_type* is not included in
    the DEFAULT_S3_PATHS.
    """
    if file_type not in DEFAULT_S3_PATHS:
        raise Exception(f"Encountered unsupported file_type {file_type!r} for S3 path.")

    return config.get(f"{file_type}_url", DEFAULT_S3_PATHS[file_type])


rule download_sequences:
    output:
        sequences="data/{lineage}/raw_{segment}.fasta"
    params:
        s3_path=_get_s3_path("raw_sequences"),
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - | xz -c -d > {output.sequences}
        """

rule download_titers:
    output:
        titers="data/{lineage}/{titer_collection}_titers.tsv"
    params:
        s3_path=_get_s3_path("titers")
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - | gzip -c -d > {output.titers}
        """

rule download_parsed_sequences:
    output:
        sequences="data/{lineage}/{segment}.fasta"
    params:
        s3_path=_get_s3_path("sequences")
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - | xz -c -d > {output.sequences}
        """

rule download_parsed_metadata:
    output:
        metadata="data/{lineage}/metadata.tsv",
    params:
        s3_path=_get_s3_path("metadata")
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} - | xz -c -d > {output.metadata}
        """

rule download_nextclade:
    output:
        nextclade="data/{lineage}/{segment}/nextclade.tsv.xz",
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{segment}/nextclade.tsv.xz"
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp {params.s3_path} {output.nextclade}
        """
