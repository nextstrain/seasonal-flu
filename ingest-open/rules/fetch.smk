"""
This part of the workflow handles fetching the metadata and sequences from
open data sources.

OUTPUTS:

    metadata      = data/{lineage}/metadata.tsv
    sequences     = data/{lineage}/{segment}.fasta

"""
from urllib.parse import urlencode


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
        metadata = "data/{lineage}/genspectrum_metadata.tsv",
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
        genspectrum_metadata = "data/{lineage}/genspectrum_metadata.tsv",
    output:
        genspectrum_to_genbank = "data/{lineage}/{segment}/genspectrum_to_genbank.tsv",
        genbank_accessions = temp("data/{lineage}/{segment}/genbank_accessions.txt"),
    benchmark:
        "benchmarks/{lineage}/{segment}/genspectrum_to_genbank.txt"
    log:
        "logs/{lineage}/{segment}/genspectrum_to_genbank.txt",
    params:
        accession_columns = lambda w: ",".join([
            "accession",
            f"insdcAccessionBase_{SEGMENT_MAP[w.segment]}",
        ]),
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk cut -t -f {params.accession_columns:q} \
            {input.genspectrum_metadata:q} \
            | csvtk filter2 -t -f 'len($2) > 0' \
            | tee {output.genspectrum_to_genbank:q} \
            | csvtk cut -t -f 2 \
            | csvtk del-header -t \
            > {output.genbank_accessions:q}
        """


# Fetch from Entrez
# Limit concurrent connections via Entrez to avoid hitting Too Many Requests error
workflow.global_resources.setdefault("concurrent_entrez", 2)
rule fetch_from_ncbi_entrez:
    input:
        genbank_accessions = "data/{lineage}/{segment}/genbank_accessions.txt",
    output:
        genbank="data/{lineage}/{segment}/genbank.gb",
    resources:
        concurrent_entrez = 1
    # Allow retries in case of network errors
    retries: 5
    benchmark:
        "benchmarks/{lineage}/{segment}/fetch_from_ncbi_entrez.txt"
    log:
        "logs/{lineage}/{segment}/fetch_from_ncbi_entrez.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        {workflow.basedir}/scripts/fetch-from-ncbi-entrez-with-accessions \
            --accessions {input.genbank_accessions:q} \
            --output {output.genbank:q}
        """


rule parse_genbank_to_tsv:
    """
    Parse relevant fields from GenBank and output as a TSV
    """
    input:
        genbank="data/{lineage}/{segment}/genbank.gb",
    output:
        tsv="data/{lineage}/{segment}/ncbi_entrez.tsv",
    benchmark:
        "benchmarks/{lineage}/{segment}/parse_genbank_to_ndjson.txt"
    log:
        "logs/{lineage}/{segment}/parse_genbank_to_ndjson.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        bio json --lines {input.genbank:q} \
            | jq -c --arg segment {wildcards.segment} '
                {{
                  "accession_\($segment)":        .record.accessions[0],
                  "strain_\($segment)":           .record.strain[0],
                  "date_\($segment)":             .record.collection_date[0],
                  "organism_\($segment)":         .record.organism[0],
                  "isolation_source_\($segment)": .record.isolation_source[0],
                  "lab_host_\($segment)":         .record.lab_host[0],
                  "note_\($segment)":             .record.note[0],
                }}
              ' \
            | augur curate passthru \
                --output-metadata {output.tsv:q}
        """


# Merge and collapse segment Entrez metadata with GenSpectrum metadata
rule link_entrez_to_genspectrum_accession:
    input:
        genspectrum_to_genbank = "data/{lineage}/{segment}/genspectrum_to_genbank.tsv",
        entrez_tsv="data/{lineage}/{segment}/ncbi_entrez.tsv",
    output:
        entrez_with_genspectrum_accession = "data/{lineage}/{segment}/entrez_with_genspectrum_accession.tsv",
    benchmark:
        "benchmarks/{lineage}/{segment}/link_entrez_to_genspectrum_accession.txt"
    log:
        "logs/{lineage}/{segment}/link_entrez_to_genspectrum_accession.txt"
    params:
        genspectrum_col = lambda w: f"insdcAccessionBase_{SEGMENT_MAP[w.segment]}",
        entrez_col = lambda w: f"accession_{w.segment}",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur merge \
            --metadata genspectrum={input.genspectrum_to_genbank:q} entrez={input.entrez_tsv:q} \
            --metadata-id-columns genspectrum={params.genspectrum_col:q} entrez={params.entrez_col:q} \
            --output-metadata {output.entrez_with_genspectrum_accession:q}
        """


rule merge_segment_metadata:
    # Using SEGMENT_MAP.keys instead of config["segments"] because we want to
    # merge the full metadata of segment records regardless of which segments
    # are requested as final output segment FASTAs.
    input:
        **{
            segment: f"data/{{lineage}}/{segment}/entrez_with_genspectrum_accession.tsv"
            for segment in SEGMENT_MAP.keys()
        },
    output:
        all_segment_metadata = "data/{lineage}/entrez_all_segment_metadata.tsv",
    benchmark:
        "benchmarks/{lineage}/merge_segment_metadata.txt"
    log:
        "logs/{lineage}/merge_segment_metadata.txt"
    params:
        metadata = lambda _, input: list(map("=".join, input.items())),
        id_field = "accession",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur merge \
            --metadata {params.metadata:q} \
            --metadata-id-columns {params.id_field:q} \
            --output-metadata {output.all_segment_metadata:q}
        """


rule collapse_segment_metadata:
    input:
        all_segment_metadata = "data/{lineage}/entrez_all_segment_metadata.tsv",
    output:
        entrez_metadata = "data/{lineage}/entrez_metadata.tsv",
    benchmark:
        "benchmarks/{lineage}/collapse_segment_metadata.txt"
    log:
        "logs/{lineage}/collapse_segment_metadata.txt"
    params:
        segments = list(SEGMENT_MAP.keys()),
        columns = ["strain", "date", "isolation_source", "note"],
    shell:
        r"""
        exec &> {log:q}

        {workflow.basedir}/scripts/collapse-segment-metadata \
            --metadata {input.all_segment_metadata:q} \
            --segments {params.segments:q} \
            --columns {params.columns:q} \
            --output-metadata {output.entrez_metadata:q}
        """

# Produce 1 metadata TSV and 8 segment FASTA
