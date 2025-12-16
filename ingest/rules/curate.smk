"""
This part of the workflow handles the curation of data from GISAID.

REQUIRED INPUTS:

    ndjson      = data/gisaid.ndjson

OUTPUTS:

    metadata    = results/{dataset}/metadata.tsv
    sequences   = results/{dataset}/{segment}.fasta

"""
from pathlib import Path


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
        lineage_annotations=config["curate"]["lineage_annotations"],
        strain_replacements_seasonal="data/fauna-source-data/flu_strain_name_fix.tsv",
        strain_replacements_avian="data/fauna-source-data/avian_flu_strain_name_fix.tsv",
        strain_location_replacements="data/fauna-source-data/flu_fix_location_label.tsv",
        strain_name_fixes=config['curate']['strain_name_fixes'],
        geolocation_rules=config["curate"]["local_geolocation_rules"],
        location_annotations=config["curate"]["location_annotations"],
        final_annotations=config["curate"]["final_annotations"],
        gisaid_location_rules=config["curate"]["gisaid_location_rules"],
    output:
        curated_ndjson=temp("data/curated_gisaid.ndjson"),
    log:
        "logs/curate.txt",
    benchmark:
        "benchmarks/curate.txt"
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        gisaid_subtype_field=config["curate"]["gisaid_subtype_field"],
        gisaid_lineage_field=config["curate"]["gisaid_lineage_field"],
        gisaid_note_field=config["curate"]["gisaid_note_field"],
        new_type_field=config["curate"]["new_type_field"],
        new_subtype_field=config["curate"]["new_subtype_field"],
        new_lineage_field=config["curate"]["new_lineage_field"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        gisaid_location_field=config["curate"]["gisaid_location_field"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        lowercase_fields=config["curate"]["lowercase"]["fields"],
        passage_field=config["curate"]["passage_field"],
        passage_category_field=config["curate"]["passage_category_field"],
        gisaid_strain_field=config["curate"]["gisaid_strain_field"],
        new_strain_field=config["curate"]["new_strain_field"],
        gihsn_field=config["curate"]["gihsn_field"],
        age_field=config["curate"]["age_field"],
        age_unit_field=config["curate"]["age_unit_field"],
        new_age_field=config["curate"]["new_age_field"],
        gender_field=config["curate"]["gender_field"],
        new_gender_field=config["curate"]["new_gender_field"],
        annotations_id=config["curate"]["annotations_id"],
        host_field=config["curate"]["host_field"],
    shell:
        r"""
        (cat {input.sequences_ndjson:q} \
            | augur curate rename \
                --field-map {params.field_map} \
            | augur curate normalize-strings \
            | ./scripts/lowercase-fields \
                --fields {params.lowercase_fields:q} \
            | ./scripts/standardize-lineage \
                --subtype-field {params.gisaid_subtype_field:q} \
                --lineage-field {params.gisaid_lineage_field:q} \
                --note-field {params.gisaid_note_field:q} \
                --new-type-field {params.new_type_field:q} \
                --new-subtype-field {params.new_subtype_field:q} \
                --new-lineage-field {params.new_lineage_field:q} \
                --annotations {input.lineage_annotations:q} \
            | augur curate format-dates \
                --date-fields {params.date_fields:q} \
                --expected-date-formats {params.expected_date_formats:q} \
            | ./scripts/parse-gisaid-location \
                --location-field {params.gisaid_location_field:q} \
                --strain-field {params.gisaid_strain_field:q} \
                --annotations {input.location_annotations:q} \
                --rules {input.gisaid_location_rules:q} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields:q} \
                --articles {params.articles:q} \
                --abbreviations {params.abbreviations:q} \
            | ./scripts/curate-host \
                --host-field {params.host_field:q} \
                --strain-field {params.gisaid_strain_field:q} \
            | augur curate apply-geolocation-rules \
            | augur curate apply-geolocation-rules \
                --no-default-rules \
                --geolocation-rules {input.geolocation_rules:q} \
            | ./scripts/annotate-with-passage-category \
                --passage-field {params.passage_field:q} \
                --passage-category-field {params.passage_category_field:q} \
            | ./scripts/standardize-strain-names \
                --strain-field {params.gisaid_strain_field:q} \
                --passage-field {params.passage_category_field:q} \
                --type-field {params.new_type_field:q} \
                --new-strain-field {params.new_strain_field:q} \
                --strain-replacements {input.strain_replacements_avian:q} {input.strain_replacements_seasonal:q} {input.strain_name_fixes:q} \
                --location-replacements {input.strain_location_replacements:q} \
            | ./scripts/annotate-with-gihsn \
                --strain-field {params.gisaid_strain_field:q} \
                --gihsn-field {params.gihsn_field:q} \
            | ./scripts/curate-age \
                --age-field {params.age_field:q} \
                --age-unit-field {params.age_unit_field:q} \
                --new-age-field {params.new_age_field:q} \
            | ./scripts/curate-gender \
                --gender-field {params.gender_field:q} \
                --new-gender-field {params.new_gender_field:q} \
            | augur curate apply-record-annotations \
                --annotations {input.final_annotations:q} \
                --id-field {params.annotations_id:q} \
                > {output.curated_ndjson:q}) 2> {log:q}
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


rule filter_curated_data:
    input:
        curated_ndjson="data/curated_gisaid.ndjson",
    output:
        filtered_ndjson=temp("data/{dataset}/curated_gisaid.ndjson"),
    log:
        "logs/{dataset}/filter_curated_data.txt",
    benchmark:
        "benchmarks/{dataset}/filter_curated_data.txt",
    params:
        gisaid_id_field=config["curate"]["gisaid_id_field"],
        optional_lineage_field=lambda w: conditional('--lineage-field', config["filtering"][w.dataset].get('lineage_field', None)),
        # TODO XXX: There's a gotcha where the lineages have to be a list if there's multiple - we should remove the :q as we know all lineage values are a single word
        lineages_to_include=lambda w: config["filtering"][w.dataset].get('lineages', None),
        optional_additional_field=lambda w:conditional('--additional-field', config["filtering"][w.dataset].get('additional_field', None)),
        optional_additional_field_values=lambda w:conditional('--additional-field-values', config["filtering"][w.dataset].get('additional_field_values', None)),
    shell:
        r"""
        cat {input.curated_ndjson:q} \
            | ./scripts/filter-ndjson \
                --id-field {params.gisaid_id_field:q} \
                {params.optional_lineage_field:q} \
                --lineages {params.lineages_to_include:q} \
                {params.optional_additional_field:q} \
                {params.optional_additional_field_values:q} \
                > {output.filtered_ndjson:q} 2> {log:q}
        """


rule deduplicate_ndjson_by_strain:
    input:
        curated_ndjson="data/{dataset}/curated_gisaid.ndjson",
    output:
        deduped_ndjson=temp("data/{dataset}/deduped_curated.ndjson"),
    log:
        "logs/{dataset}/deduplicate_ndjson_by_strain.txt"
    params:
        strain_field=config["curate"]["new_strain_field"],
        id_field=config["curate"]["gisaid_id_field"],
        # TODO XXX - can we make this an input to get snakemake's file checking?
        prioritized_strain_ids=lambda w: conditional('--prioritized-ids', config["filtering"][w.dataset].get('prioritized_strain_ids', None)),
    shell:
        r"""
         cat {input.curated_ndjson:q} \
            | ./scripts/dedup-by-strain \
                --strain-field {params.strain_field:q} \
                --id-field {params.id_field:q} \
                {params.prioritized_strain_ids:q} \
                > {output.deduped_ndjson:q} 2> {log:q}
        """


rule split_ndjson_by_segment:
    input:
        deduped_ndjson="data/{dataset}/deduped_curated.ndjson",
    output:
        metadata="data/{dataset}/curated_metadata.tsv",
        sequences=expand("results/{{dataset}}/{segment}.fasta", segment=config["segments"]),
    log:
        "logs/{dataset}/split_ndjson_by_segment.txt"
    params:
        segments=config["segments"],
        seq_output_dir=lambda w, output: Path(output.sequences[0]).parent,
        id_field=config["curate"]["output_id_field"],
    shell:
        r"""
        cat {input.deduped_ndjson:q} \
            | ./scripts/split-gisaid-ndjson-by-segment \
                --segments {params.segments:q} \
                --output-metadata {output.metadata:q} \
                --sequences-output-dir {params.seq_output_dir:q} \
                --output-id-field {params.id_field:q} 2> {log:q}
        """


# Modified from top level rule
# <https://github.com/nextstrain/seasonal-flu/blob/f073d3e055ab6efaae6d4c91a06efa621b6247d9/workflow/snakemake_rules/select_strains.smk#L94-L113>
rule build_reference_strains_table:
    input:
        references=lambda w: config['filtering'][w.dataset]['reference_strains'],
    output:
        references="data/{dataset}/reference_strains.tsv",
    params:
        reference_column=config["curate"]["reference_column"],
        id_field=config["curate"]["output_id_field"],
    shell:
        r"""
        csvtk add-header \
            --names {params.id_field:q} \
            {input.references:q} \
            | csvtk uniq \
            | csvtk --out-tabs mutate2 \
                --name {params.reference_column:q} \
                --expression "'True'" > {output.references:q}
        """


# Modified from top level rule
# https://github.com/nextstrain/seasonal-flu/blob/f073d3e055ab6efaae6d4c91a06efa621b6247d9/workflow/snakemake_rules/select_strains.smk#L115-L137
rule annotate_metadata_with_reference_strains:
    input:
        references="data/{dataset}/reference_strains.tsv",
        metadata="data/{dataset}/curated_metadata.tsv",
    output:
        metadata="data/{dataset}/curated_metadata_with_references.tsv",
    params:
        id_field=config["curate"]["output_id_field"],
    shell:
        r"""
        csvtk -t join \
            --left-join \
            --na "False" \
            -f {params.id_field:q} \
            {input.metadata:q} \
            {input.references:q} > {output.metadata:q}
        """


def metadata_selector(wildcards):
    """
    If the config sets a "reference_strains" TSV file for the dataset then this
    function returns the metadata TSV file which has the reference strain appended.
    Otherwise this points to an earlier TSV file, thus skipping the rules which
    add the reference strains.
    """
    if config['filtering'][wildcards.dataset].get('reference_strains', False):
        return "data/{dataset}/curated_metadata_with_references.tsv",
    return "data/{dataset}/curated_metadata.tsv",


def metadata_fields(wildcards) -> str:
    """
    Returns config defined columns and any additional segment
    columns added by ./scripts/split-gisaid-ndjson-by-segment
    """
    metadata_columns = config["filtering"][wildcards.dataset]["metadata_columns"].copy()
    for segment in config["segments"]:
        metadata_columns.extend([
            segment,
            f"accession_{segment}",
        ])
    metadata_columns.append("n_segments")
    return ",".join(metadata_columns)


rule subset_metadata:
    input:
        metadata=metadata_selector,
    output:
        subset_metadata="results/{dataset}/metadata.tsv",
    params:
        metadata_fields=metadata_fields,
    shell:
        r"""
        csvtk cut -t -f {params.metadata_fields:q} \
            {input.metadata:q} > {output.subset_metadata:q}
        """
