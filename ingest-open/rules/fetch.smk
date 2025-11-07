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
rule genspectrum_to_genbank:
    input:
        genspectrum_metadata = "data/{lineage}/genspectrum/metadata.tsv",
    output:
        genbank_accessions = "data/{lineage}/genspectrum_to_genbank.tsv"
    benchmark:
        "benchmarks/{lineage}/genspectrum_to_genbank.txt"
    log:
        "logs/{lineage}/genspectrum_to_genbank.txt",
    params:
        accession_columns = ",".join([
            "accession",
            "insdcAccessionBase_seg1",
            "insdcAccessionBase_seg2",
            "insdcAccessionBase_seg3",
            "insdcAccessionBase_seg4",
            "insdcAccessionBase_seg5",
            "insdcAccessionBase_seg6",
            "insdcAccessionBase_seg7",
            "insdcAccessionBase_seg8",
        ]),
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk cut -t -f {params.accession_columns:q} \
            {input.genspectrum_metadata:q} \
            > {output.genbank_accessions:q}
        """


rule gather_genbank_accessions:
    input:
        genbank_accessions = "data/{lineage}/genspectrum_to_genbank.tsv",
    output:
        genbank_accessions = "data/{lineage}/genbank_accessions.txt",
    benchmark:
        "benchmarks/{lineage}/gather_genbank_accessions.txt"
    log:
        "logs/{lineage}/gather_genbank_accessions.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.genbank_accessions:q} \
            | csvtk gather -t -k item -v value -f -1 \
            | csvtk filter2 -t -f 'len($value) > 0' \
            | csvtk cut -t -f value \
            | csvtk del-header -t > {output.genbank_accessions:q}
        """


# Fetch from Entrez
rule fetch_from_ncbi_entrez:
    input:
        genbank_accessions = "data/{lineage}/genbank_accessions.txt",
    output:
        genbank="data/{lineage}/genbank.gb",
    # Allow retries in case of network errors
    retries: 5
    benchmark:
        "benchmarks/{lineage}/fetch_from_ncbi_entrez.txt"
    log:
        "logs/{lineage}/fetch_from_ncbi_entrez.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        {workflow.basedir}/scripts/fetch-from-ncbi-entrez-with-accessions \
            --accessions {input.genbank_accessions:q} \
            --output {output.genbank:q}
        """
# Merge and collapse segment Entrez metadata with GenSpectrum metadata
# Produce 1 metadata TSV and 8 segment FASTA
