"""
This part of the workflow handles the curation of the open data

REQUIRED INPUTS:

    ndjson = "data/{lineage}/open.ndjson.zst"

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
        ndjson="data/{lineage}/open.ndjson.zst",
        geolocation_rules=resolve_config_path(config["curate"]["local_geolocation_rules"]),
        annotations=resolve_config_path(config["curate"]["annotations"]),
    output:
        ndjson="data/{lineage}/curated.ndjson.zst",
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        strain_regex=config["curate"]["strain_regex"],
        strain_backup_fields=config["curate"]["strain_backup_fields"],
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

        zstdcat {input.ndjson:q} \
            | augur curate rename \
                --field-map {params.field_map:q} \
            | augur curate normalize-strings \
            | augur curate transform-strain-name \
                --strain-regex {params.strain_regex:q} \
                --backup-fields {params.strain_backup_fields:q} \
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
            | zstd -T0 -c > {output.ndjson:q}
        """


# Copied from
# <https://github.com/nextstrain/zika/blob/4d8d8c452b8f3f50b771d95a17662c88818b6048/phylogenetic/rules/annotate_phylogeny.smk#L81-L97>
def conditional(option, argument):
    """Used for config-defined arguments whose presence necessitates a command-line option
    (e.g. --foo) prepended and whose absence should result in no option/arguments in the CLI command.
    *argument* can be falsey, in which case an empty string is returned (i.e. "don't pass anything
    to the CLI"), or a *list* or *string* or *number* in which case a flat list of options/args is returned,
    or *True* in which case a list of a single element (the option) is returned.
    Any other argument type is a WorkflowError
    """
    if not argument:
        return ""
    if argument is True: # must come before `isinstance(argument, int)` as bool is a subclass of int
        return [option]
    if isinstance(argument, list):
        return [option, *argument]
    if isinstance(argument, int) or isinstance(argument, float) or isinstance(argument, str):
        return [option, argument]
    raise WorkflowError(f"Workflow function conditional() received an argument value of unexpected type: {type(argument).__name__}")


def get_prioritized_id_file(wildcards):
    if prioritized_id_file := config["curate"].get("prioritized_strain_ids", {}).get(wildcards.lineage):
        return resolve_config_path(prioritized_id_file)
    return []

rule prioritize_id_per_strain:
    input:
        curated_ndjson="data/{lineage}/curated.ndjson.zst",
        prioritized_strain_ids=get_prioritized_id_file,
    output:
        prioritized_ids="data/{lineage}/prioritized_ids.txt",
    benchmark:
        "benchmarks/{lineage}/prioritize_id_per_strain.txt"
    log:
        "logs/{lineage}/prioritize_id_per_strain.txt"
    params:
        strain_field=config["curate"]["strain_field"],
        id_field=config["curate"]["record_id_field"],
        seq_field="sequences",
        prioritized_strain_ids=lambda _, input: conditional('--prioritized-ids', input.prioritized_strain_ids),
    shell:
        r"""
        exec &> >(tee {log:q})

        zstdcat {input.curated_ndjson:q} \
            | {workflow.basedir}/../ingest/scripts/prioritize-id-per-strain \
                --strain-field {params.strain_field:q} \
                --id-field {params.id_field:q} \
                --seq-field {params.seq_field:q} \
                {params.prioritized_strain_ids:q} \
                --output {output.prioritized_ids:q}
        """


rule deduplicate_ndjson_by_strain:
    input:
        curated_ndjson="data/{lineage}/curated.ndjson.zst",
        prioritized_ids="data/{lineage}/prioritized_ids.txt",
    output:
        deduped_ndjson=temp("data/{lineage}/deduped_curated.ndjson.zst"),
    benchmark:
        "benchmarks/{lineage}/deduplicate_ndjson_by_strain.txt"
    log:
        "logs/{lineage}/deduplicate_ndjson_by_strain.txt"
    params:
        id_field=config["curate"]["record_id_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

        zstdcat {input.curated_ndjson:q} \
            | {workflow.basedir}/../ingest/scripts/filter-ndjson-by-value \
                --field {params.id_field:q} \
                --include {input.prioritized_ids:q} \
            | zstd -T0 -c > {output.deduped_ndjson:q}
        """


rule split_ndjson_by_segment:
    input:
        deduped_ndjson="data/{lineage}/deduped_curated.ndjson.zst",
    output:
        metadata="data/{lineage}/curated_metadata.tsv",
        sequences=expand("results/{{lineage}}/{segment}.fasta", segment=config["segments"]),
    benchmark:
        "benchmarks/{lineage}/split_ndjson_by_segment.txt"
    log:
        "logs/{lineage}/split_ndjson_by_segment.txt"
    params:
        segments=config["segments"],
        seq_output_dir=lambda w, output: Path(output.sequences[0]).parent,
        id_field=config["curate"]["record_id_field"],
        select_seq="error",
    shell:
        r"""
        exec &> >(tee {log:q})

        zstdcat {input.deduped_ndjson:q} \
            | {workflow.basedir}/../ingest/scripts/split-ndjson-by-segment \
                --segments {params.segments:q} \
                --select-seq {params.select_seq:q} \
                --output-metadata {output.metadata:q} \
                --sequences-output-dir {params.seq_output_dir:q} \
                --output-id-field {params.id_field:q}
        """


def metadata_fields(wildcards) -> str:
    """
    Returns config defined columns and any additional segment
    columns added by ./scripts/split-ndjson-by-segment
    """
    metadata_columns = config["curate"]["metadata_columns"].copy()
    for segment in config["segments"]:
        metadata_columns.extend([
            segment,
            f"accession_{segment}",
        ])
    metadata_columns.append("n_segments")
    return ",".join(metadata_columns)


rule subset_metadata:
    input:
        metadata="data/{lineage}/curated_metadata.tsv",
    output:
        subset_metadata="results/{lineage}/metadata.tsv",
    params:
        metadata_fields=metadata_fields,
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
