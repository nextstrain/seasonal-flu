"""
This part of the workflow handles the curation of data from GISAID.
Requires the `SEGMENTS` variable to be set upstream in the workflow.

REQUIRED INPUTS:

    ndjson      = data/gisaid.ndjson

OUTPUTS:

    metadata    = data/subset_metadata.tsv
    sequences   = results/sequences.fasta

"""


def format_field_map(field_map: dict[str, str]) -> str:
    """
    Format dict to `"key1"="value1" "key2"="value2"...` for use in shell commands.
    """
    return " ".join([f'"{key}"="{value}"' for key, value in field_map.items()])


rule curate:
    input:
        sequences_ndjson="data/gisaid.ndjson",
        geolocation_rules=config["curate"]["local_geolocation_rules"],
        annotations=config["curate"]["annotations"],
    output:
        curated_ndjson="data/curated_gisaid.ndjson",
    log:
        "logs/curate.txt",
    benchmark:
        "benchmarks/curate.txt"
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        annotations_id=config["curate"]["annotations_id"],
    shell:
        r"""
        (cat {input.sequences_ndjson:q} \
            | augur curate rename \
                --field-map {params.field_map} \
            | augur curate normalize-strings \
            | augur curate format-dates \
                --date-fields {params.date_fields:q} \
                --expected-date-formats {params.expected_date_formats:q} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields:q} \
                --articles {params.articles:q} \
                --abbreviations {params.abbreviations:q} \
            | augur curate apply-geolocation-rules \
                --geolocation-rules {input.geolocation_rules:q} \
            | augur curate apply-record-annotations \
                --annotations {input.annotations:q} \
                --id-field {params.annotations_id:q} \
                > {output.curated_ndjson}) 2>> {log}
        """


# TODO: Split ndjson by subtype before splitting into metadata + FASTAs
rule split_ndjson_by_segment:
    input:
        curated_ndjson="data/curated_gisaid.ndjson",
    output:
        metadata="data/all_metadata.tsv",
        sequences=expand("results/{segment}/sequences.fasta", segment=SEGMENTS),
    params:
        segments=SEGMENTS,
        seq_output_dir=lambda w, output: output.sequences[0].split("/")[0],
        id_field=config["curate"]["output_id_field"],
    shell:
        r"""
        cat {input.curated_ndjson:q} \
            | ./scripts/split-gisaid-ndjson-by-segment \
                --segments {params.segments:q} \
                --output-metadata {output.metadata:q} \
                --sequences-output-dir {params.seq_output_dir:q} \
                --output-id-field {params.id_field:q}
        """


rule subset_metadata:
    input:
        metadata="data/all_metadata.tsv",
    output:
        subset_metadata="results/metadata.tsv",
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    shell:
        r"""
        csvtk cut -t -f {params.metadata_fields:q} \
            {input.metadata:q} > {output.subset_metadata:q}
        """
