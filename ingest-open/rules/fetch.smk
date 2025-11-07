"""
This part of the workflow handles fetching the metadata and sequences from
open data sources.

OUTPUTS:

    metadata      = data/{lineage}/metadata.tsv
    sequences     = data/{lineage}/{segment}.fasta

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
    # Map wildcards.segment to GenSpectrum segments
    SEGMENT_MAP = {
        "pb2": "seg1",
        "pb1": "seg2",
        "pa":  "seg3",
        "ha":  "seg4",
        "np":  "seg5",
        "na":  "seg6",
        "mp":  "seg7",
        "ns":  "seg8",
    }

    if (segment := SEGMENT_MAP.get(wildcards.segment)) is None:
        raise InvalidConfigError(
            f"Encountered unsupported segment {wildcards.segment!r}. "
            f"Segment must be one of {list(SEGMENT_MAP.keys())}")

    endpoint = f"{_genspectrum_lapis_url(wildcards)}/sample/unalignedNucleotideSequences/{segment}"
    params = {
        "dataFormat": "FASTA",
        "fastaHeaderTemplate": "{accession}",
        "versionStatus": "LATEST_VERSION"
    }
    query = urlencode(params, doseq=True, encoding="utf-8")
    return f"{endpoint}?{query}"


# Fetch metadata from GenSpectrum
rule download_genspectrum_metadata:
    output:
        metadata = "data/{lineage}/metadata.tsv",
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
        sequences = "data/{lineage}/{segment}.fasta",
    retries: 5
    benchmark:
        "benchmarks/{lineage}/{segment}/download_genspectrum_sequences.txt",
    log:
        "logs/{lineage}/{segment}/download_genspectrum_sequences.txt",
    params:
        sequences_url = lambda w: _genspectrum_sequences_url(w)
    shell:
        r"""
        exec &> >(tee {log:q})

        curl -fsSL {params.sequences_url:q} -o {output.sequences:q}
        """


# Pull out GenBank accessions from GenSpectrum metadata
# Fetch from Entrez
# Merge and collapse segment Entrez metadata with GenSpectrum metadata
# Produce 1 metadata TSV and 8 segment FASTA
