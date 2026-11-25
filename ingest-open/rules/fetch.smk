"""
This part of the workflow handles fetching the metadata and sequences from
open data sources.

OUTPUTS:

    ndjson = "data/{lineage}/open.ndjson.zst"

"""
from urllib.parse import urlencode


def _genspectrum_lapis_url(wildcards):
    """
    Returns the URL for the GenSpectrum LAPIS query engine for
    the provided *wildcards.lineage*.

    See <https://loculus.genspectrum.org/api-documentation>
    """
    # Map wildcards.lineage to GenSpectrum lineages
    LINEAGE_MAP = {
        "h1n1pdm": "h1n1pdm",
        "h3n2":    "h3n2",
        "vic":     "b-victoria",
    }

    if (lineage := LINEAGE_MAP.get(wildcards.lineage)) is None:
        raise InvalidConfigError(
            f"Encountered unsupported lineage {wildcards.lineage!r}. "
            f"Lineage must be one of {list(LINEAGE_MAP.keys())}")

    return f"https://api.loculus.genspectrum.org/{lineage}"


def _genspectrum_metadata_url(wildcards):
    """
    Returns the URL for downloading the metadata TSV.

    See <https://api.loculus.genspectrum.org/b-victoria/swagger-ui/index.html#/lapis-controller/getDetailsAsCsv>
    """
    endpoint = f"{_genspectrum_lapis_url(wildcards)}/sample/details"
    params = {
        "dataFormat": "TSV",
        "versionStatus": "LATEST_VERSION",
        "compression": "zstd",
        "fields": config["genspectrum_metadata_fields"],
        **config["genspectrum_filters"],
    }
    query = urlencode(params, doseq=True, encoding="utf-8")
    return f"{endpoint}?{query}"


def _genspectrum_sequences_url(wildcards):
    """
    Returns the URL for downloading the sequences FASTA per *wildcards.segment*.
    Only includes the LATEST_VERSION of GenSpectrum records and uses the
    GenSpectrum accession field as the FASTA header.

    See <https://api.loculus.genspectrum.org/b-victoria/swagger-ui/index.html#/multi-segmented-sequence-controller/getUnalignedNucleotideSequence>
    """
    endpoint = f"{_genspectrum_lapis_url(wildcards)}/sample/unalignedNucleotideSequences"
    params = {
        "dataFormat": "FASTA",
        "fastaHeaderTemplate": config["genspectrum_fastaHeaderTemplate"],
        "versionStatus": "LATEST_VERSION",
        "compression": "zstd",
        **config["genspectrum_filters"],
    }
    query = urlencode(params, doseq=True, encoding="utf-8")
    return f"{endpoint}?{query}"


# Fetch metadata from GenSpectrum
rule download_genspectrum_metadata:
    output:
        metadata = "data/{lineage}/metadata.tsv.zst",
    retries: 5
    benchmark:
        "benchmarks/{lineage}/download_genspectrum_metadata.txt",
    log:
        "logs/{lineage}/download_genspectrum_metadata.txt",
    params:
        metadata_url = lambda w: _genspectrum_metadata_url(w)
    shell:
        r"""
        exec &> >(tee {log:q})

        curl -fsSL {params.metadata_url:q} -o {output.metadata:q}
        """


# Fetch sequences from GenSpectrum
rule download_genspectrum_sequences:
    output:
        sequences = "data/{lineage}/sequences.fasta.zst",
    retries: 5
    benchmark:
        "benchmarks/{lineage}/download_genspectrum_sequences.txt",
    log:
        "logs/{lineage}/download_genspectrum_sequences.txt",
    params:
        sequences_url = lambda w: _genspectrum_sequences_url(w)
    shell:
        r"""
        exec &> >(tee {log:q})

        curl -fsSL {params.sequences_url:q} -o {output.sequences:q}
        """


rule decompress_sequences:
    input:
        sequences = "data/{lineage}/sequences.fasta.zst"
    output:
        sequences = temp("data/{lineage}/sequences.fasta")
    benchmark:
        "benchmarks/{lineage}/decompress_sequences.txt",
    log:
        "logs/{lineage}/decompress_sequences.txt",
    shell:
        r"""
        zstd -dcq {input.sequences:q} > {output.sequences:q}
        """


rule link_metadata_and_sequences:
    input:
        metadata = "data/{lineage}/metadata.tsv.zst",
        sequences = "data/{lineage}/sequences.fasta",
    output:
        ndjson = "data/{lineage}/open.ndjson.zst",
    benchmark:
        "benchmarks/{lineage}/link_metadata_and_sequences.txt",
    log:
        "logs/{lineage}/link_metadata_and_sequences.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        {workflow.basedir}/scripts/link-metadata-and-sequences \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            | zstd -T0 -c > {output.ndjson:q}
        """
