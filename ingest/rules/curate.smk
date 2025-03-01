"""
This part of the workflow handles the curation of data from GISAID.
Requires the `SEGMENTS` variable to be set upstream in the workflow.

REQUIRED INPUTS:

    ndjson      = data/gisaid.ndjson

OUTPUTS:

    metadata    = data/subset_metadata.tsv
    sequences   = results/sequences.fasta

"""


rule fetch_from_fauna:
    """
    Fetch files from nextstrain/fauna/source-data that we need to ensure the
    strain names of the sequences are in-sync with strain names in titer data.

    Currently used to fetch the following files for rule curate:
    - flu_strain_name_fix.tsv
    - flu_fix_location_label.tsv
    """
    output: "data/fauna-source-data/{source_data_file}",
    params:
        url=lambda w: f"https://github.com/nextstrain/fauna/raw/@/source-data/{w.source_data_file}",
    shell:
        """
        curl -fsSL {params.url} > {output}
        """


def format_field_map(field_map: dict[str, str]) -> str:
    """
    Format dict to `"key1"="value1" "key2"="value2"...` for use in shell commands.
    """
    return " ".join([f'"{key}"="{value}"' for key, value in field_map.items()])


rule curate:
    input:
        sequences_ndjson="data/gisaid.ndjson",
        strain_replacements="data/fauna-source-data/flu_strain_name_fix.tsv",
        strain_location_replacements="data/fauna-source-data/flu_fix_location_label.tsv",
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
        gisaid_location_field=config["curate"]["gisaid_location_field"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        passage_field=config["curate"]["passage_field"],
        passage_category_field=config["curate"]["passage_category_field"],
        gisaid_subtype_field=config["curate"]["gisaid_subtype_field"],
        gisaid_lineage_field=config["curate"]["gisaid_lineage_field"],
        new_type_field=config["curate"]["new_type_field"],
        new_subtype_field=config["curate"]["new_subtype_field"],
        new_lineage_field=config["curate"]["new_lineage_field"],
        gisaid_strain_field=config["curate"]["gisaid_strain_field"],
        new_strain_field=config["curate"]["new_strain_field"],
        gihsn_field=config["curate"]["gihsn_field"],
        age_field=config["curate"]["age_field"],
        age_unit_field=config["curate"]["age_unit_field"],
        new_age_field=config["curate"]["new_age_field"],
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
            | ./scripts/parse-gisaid-location \
                --location-field {params.gisaid_location_field:q} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields:q} \
                --articles {params.articles:q} \
                --abbreviations {params.abbreviations:q} \
            | augur curate apply-geolocation-rules \
                --geolocation-rules {input.geolocation_rules:q} \
            | ./scripts/annotate-with-passage-category \
                --passage-field {params.passage_field:q} \
                --passage-category-field {params.passage_category_field:q} \
            | ./scripts/standardize-lineage \
                --subtype-field {params.gisaid_subtype_field:q} \
                --lineage-field {params.gisaid_lineage_field:q} \
                --new-type-field {params.new_type_field:q} \
                --new-subtype-field {params.new_subtype_field:q} \
                --new-lineage-field {params.new_lineage_field:q} \
            | ./scripts/standardize-strain-names \
                --strain-field {params.gisaid_strain_field:q} \
                --passage-field {params.passage_category_field:q} \
                --type-field {params.new_type_field:q} \
                --new-strain-field {params.new_strain_field:q} \
                --strain-replacements {input.strain_replacements:q} \
                --location-replacements {input.strain_location_replacements:q} \
            | ./scripts/annotate-with-gihsn \
                --strain-field {params.gisaid_strain_field:q} \
                --gihsn-field {params.gihsn_field:q} \
            | ./scripts/curate-age \
                --age-field {params.age_field:q} \
                --age-unit-field {params.age_unit_field:q} \
                --new-age-field {params.new_age_field:q} \
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
