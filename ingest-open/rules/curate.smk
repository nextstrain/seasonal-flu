"""
This part of the workflow handles the curation of the open data

REQUIRED INPUTS:

    metadata      = data/{lineage}/metadata.tsv
    sequences     = data/{lineage}/{segment}.fasta

OUTPUTS:

    metadata      = results/{lineage}/metadata.tsv
    sequences     = results/{lineage}/{segment}.fasta

"""

# Curate metadata
def format_field_map(field_map: dict[str, str]) -> list[str]:
    """
    Format entries to the format expected by `augur curate --field-map`.
    When used in a Snakemake shell block, the list is automatically expanded and
    spaces are handled by quoted interpolation.
    """
    return [f'{key}={value}' for key, value in field_map.items()]

            
rule curate:
    input:
        metadata="data/{lineage}/metadata.tsv",
        geolocation_rules=config["curate"]["local_geolocation_rules"],
        annotations=config["curate"]["annotations"],
    output:
        metadata="data/{lineage}/curated_metadata.tsv",
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        division_field=config["curate"]["genspectrum_division_field"],
        location_field=config["curate"]["location_field"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        authors_field=config["curate"]["authors_field"],
        authors_default_value=config["curate"]["authors_default_value"],
        abbr_authors_field=config["curate"]["abbr_authors_field"],
        annotations_id=config["curate"]["annotations_id"],
    benchmark:
        "benchmarks/{lineage}/curate.txt"
    log:
        "logs/{lineage}/curate.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur curate rename \
            --metadata {input.metadata:q} \
            --field-map {params.field_map:q} \
            | augur curate normalize-strings \
            | augur curate format-dates \
                --date-fields {params.date_fields:q} \
                --expected-date-formats {params.expected_date_formats:q} \
            | {workflow.basedir}/scripts/parse-genspectrum-division \
                --division-field {params.division_field:q} \
                --location-field {params.location_field:q} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields:q} \
                --articles {params.articles:q} \
                --abbreviations {params.abbreviations:q} \
            | augur curate abbreviate-authors \
                --authors-field {params.authors_field:q} \
                --default-value {params.authors_default_value:q} \
                --abbr-authors-field {params.abbr_authors_field:q} \
            | augur curate apply-geolocation-rules \
            | augur curate apply-geolocation-rules \
                --no-default-rules \
                --geolocation-rules {input.geolocation_rules:q} \
            | augur curate apply-record-annotations \
                --annotations {input.annotations:q} \
                --id-field {params.annotations_id:q} \
                --output-metadata {output.metadata:q}
        """


rule subset_metadata:
    input:
        metadata="data/{lineage}/curated_metadata.tsv",
    output:
        subset_metadata="results/{lineage}/metadata.tsv",
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    benchmark:
        "benchmarks/{lineage}/subset_metadata.txt"
    log:
        "logs/{lineage}/subset_metadata.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk cut -t -f {params.metadata_fields:q} \
            {input.metadata:q} > {output.subset_metadata:q}
        """

# Deduplicate by strain
# Filter FASTAs by deduped metadata
# Replace seq ids in FASTAs with strain
