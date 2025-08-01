'''
This file contains rules that select strains for a build, extracts the sequences
and the meta data subset.

input:
 - "data/{lineage}/raw_{segment}.fasta" # e.g., from fauna or GISAID download with metadata in FASTA headers.
 - "data/{lineage}/{center}_{passage}_{assay}_titers.tsv" # e.g., from fauna

intermediate files:
 - "data/{lineage}/metadata.tsv"
 - "data/{lineage}/{segment}.fasta"

output:
 - builds/{build_name}/strains.txt
 - builds/{build_name}/metadata.tsv
 - builds/{build_name}/titers.tsv
 - builds/{build_name}/{segment}/sequences.fasta
'''

localrules: titer_priorities, select_titers

build_dir = config.get("build_dir", "builds")

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = "data/{lineage}/raw_{segment}.fasta",
    output:
        sequences = "data/{lineage}/{segment}.fasta",
        metadata = "data/{lineage}/metadata_{segment}.tsv",
    params:
        fasta_fields=config.get("fasta_fields", ""),
        prettify_fields_arg=lambda wildcards: f"--prettify-fields {' '.join(config['prettify_fields'])}" if "prettify_fields" in config else "",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/parse_{lineage}_{segment}.txt"
    log:
        "logs/parse_{lineage}_{segment}.txt"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields} \
            {params.prettify_fields_arg} 2>&1 | tee {log}
        """

rule join_metadata:
    input:
        segment_metadata=lambda w: [f"data/{w.lineage}/metadata_{segment}.tsv" for segment in config['segments']],
    output:
        metadata="data/{lineage}/metadata_joined.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/join_metadata_{lineage}.txt"
    log:
        "logs/join_metadata_{lineage}.txt"
    params:
        segments=lambda w: config['segments'],
        segment_columns=["accession"],
        how="outer",
    shell:
        """
        python3 scripts/join_metadata.py \
            --metadata {input.segment_metadata:q} \
            --segments {params.segments:q} \
            --segment-columns {params.segment_columns:q} \
            --how {params.how:q} \
            --output {output.metadata:q} 2>&1 | tee {log}
        """

# Annotate strains in the metadata that have "GIHSN" in the strain name to
# indicate whether it was collected as part of the Global Influenza Hospital
# Surveillance Network (GIHSN)
rule annotate_metadata_with_gihsn:
    input:
        metadata="data/{lineage}/metadata_joined.tsv",
    output:
        metadata=temp("data/{lineage}/metadata_with_gihsn.tsv"),
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/annotate_metadata_with_gihsn_{lineage}.txt"
    log:
        "logs/annotate_metadata_with_gihsn_{lineage}.txt"
    shell:
        """
        csvtk --tabs mutate2 \
            --expression '${{strain}}=~"(GIHSN)" ? "True" : "False"' \
            --name gihsn_sample \
            {input.metadata} > {output.metadata}
        """

rule build_reference_strains_table:
    input:
        references="config/{lineage}/reference_strains.txt",
    output:
        references="data/{lineage}/reference_strains.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/build_reference_strains_table_{lineage}.txt"
    log:
        "logs/build_reference_strains_table_{lineage}.txt"
    shell:
        """
        csvtk add-header \
            --names strain \
            {input.references} \
            | csvtk uniq \
            | csvtk --out-tabs mutate2 \
                --name is_reference \
                --expression "'True'" > {output.references}
        """

# Annotate strains in the metadata based on whether they are reference strains
# or not, so we can subsample these strains by attribute from augur filter
# later.
rule annotate_metadata_with_reference_strains:
    input:
        metadata="data/{lineage}/metadata_with_gihsn.tsv",
        references="data/{lineage}/reference_strains.tsv",
    output:
        metadata="data/{lineage}/metadata.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/annotate_metadata_with_reference_strains_{lineage}.txt"
    log:
        "logs/annotate_metadata_with_reference_strains_{lineage}.txt"
    shell:
        """
        csvtk --tabs join \
            --left-join \
            --na "False" \
            -f "strain" \
            {input.metadata} \
            {input.references} > {output.metadata}
        """

rule concat_titers_for_build:
    input:
        titers=lambda wildcards: [collection["data"] for collection in config["builds"][wildcards.build_name]["titer_collections"]],
    output:
        titers="builds/{build_name}/all_titers.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/concat_titers_for_build_{build_name}.txt"
    log:
        "logs/concat_titers_for_build_{build_name}.txt"
    shell:
        """
        tsv-append -H {input.titers} | tsv-select -H -f virus_strain,serum_strain > {output.titers} 2> {log}
        """

rule titer_priorities:
    input:
        titers = "builds/{build_name}/all_titers.tsv",
    output:
        priorities = build_dir + "/{build_name}/titer_priorities.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/titer_priorities_{build_name}.txt"
    log:
        "logs/titer_priorities_{build_name}.txt"
    shell:
        """
        tsv-summarize -H --group-by virus_strain --count {input.titers} | sed 1d > {output.priorities} 2> {log}
        """

rule build_titer_strains_table:
    input:
        titers="builds/{build_name}/all_titers.tsv",
    output:
        titer_strains=build_dir + "/{build_name}/titer_strains.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/build_titer_strains_table_{build_name}.txt"
    log:
        "logs/build_titer_strains_table_{build_name}.txt"
    shell:
        """
        csvtk --tabs cut \
            --fields virus_strain \
            {input.titers} \
            | csvtk rename \
              --fields virus_strain \
              --names strain \
            | csvtk uniq \
            | csvtk --out-tabs mutate2 \
                --name is_titer_strain \
                --expression "'True'" > {output.titer_strains}
        """

rule build_titer_reference_strains_table:
    input:
        titers="builds/{build_name}/all_titers.tsv",
    output:
        titer_strains=build_dir + "/{build_name}/titer_reference_strains.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/build_titer_reference_strains_table_{build_name}.txt"
    log:
        "logs/build_titer_reference_strains_table_{build_name}.txt"
    shell:
        """
        csvtk --tabs cut \
            --fields serum_strain \
            {input.titers} \
            | csvtk rename \
              --fields serum_strain \
              --names strain \
            | csvtk uniq \
            | csvtk --out-tabs mutate2 \
                --name is_titer_reference_strain \
                --expression "'True'" > {output.titer_strains}
        """

# Annotate strains in the metadata based on whether they have titer data or not,
# so we can include these strains by attribute from augur filter later.
rule annotate_metadata_with_titer_strains:
    input:
        metadata=lambda wildcards: f"data/{config['builds'][wildcards.build_name]['lineage']}/metadata.tsv",
        titer_strains=build_dir + "/{build_name}/titer_strains.tsv",
        titer_reference_strains=build_dir + "/{build_name}/titer_reference_strains.tsv",
    output:
        metadata=build_dir + "/{build_name}/full_metadata_with_titer_annotations.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/annotate_metadata_with_titer_strains_{build_name}.txt"
    log:
        "logs/annotate_metadata_with_titer_strains_{build_name}.txt"
    shell:
        """
        csvtk --tabs join --left-join --na "False" -f "strain" {input.metadata} {input.titer_strains} \
            | csvtk --tabs join --left-join --na "False" -f "strain" /dev/stdin {input.titer_reference_strains} > {output.metadata}
        """

# Run Nextclade, if we don't already have access to Nextclade annotations from
# elsewhere (e.g., S3).
rule get_nextclade_dataset_for_lineage_and_segment:
    output:
        nextclade_dir=directory("nextclade_dataset/{lineage}_{segment}/"),
    params:
        nextclade_server_arg=lambda wildcards: f"--server={config['nextclade_server']}" if config.get("nextclade_server") else "",
    shell:
        r"""
        nextclade3 dataset get \
            -n 'nextstrain/flu/{wildcards.lineage}/{wildcards.segment}' \
            {params.nextclade_server_arg} \
            --output-dir {output.nextclade_dir}
        """

rule run_nextclade:
    input:
        nextclade_dir="nextclade_dataset/{lineage}_{segment}/",
        sequences="data/{lineage}/{segment}.fasta",
    output:
        annotations="data/{lineage}/{segment}/nextclade.tsv.xz",
    log:
        "logs/run_nextclade_{lineage}_{segment}.txt"
    threads: 8
    shell:
        r"""
        nextclade3 run \
            -j {threads} \
            -D {input.nextclade_dir} \
            --output-tsv {output.annotations} \
            {input.sequences}
        """

def get_metadata_for_nextclade_merge(wildcards):
    # Use metadata annotated with a given build's titer strains, if we are
    # building the measurements panel or running titer models.
    if config['builds'][wildcards.build_name].get("enable_measurements") or config['builds'][wildcards.build_name].get("enable_titer_models"):
        return f"{build_dir}/{wildcards.build_name}/full_metadata_with_titer_annotations.tsv"
    else:
        return f"data/{config['builds'][wildcards.build_name]['lineage']}/metadata.tsv"

def get_nextclade_for_subsampling(wildcards):
    return f"data/{config['builds'][wildcards.build_name]['lineage']}/ha/nextclade.tsv.xz"

rule merge_nextclade_with_metadata:
    """
    Nextclade data is either merged with regular metadata (1) or titered metadata (2), see functions above
    """
    input:
        metadata=get_metadata_for_nextclade_merge,
        nextclade=get_nextclade_for_subsampling,
    output:
        merged = "{build_dir}/{build_name}/metadata_with_nextclade.tsv"
    params:
        metadata_id="strain",
        nextclade_id="seqName",
    conda: "../envs/nextstrain.yaml"
    log:
        "logs/{build_dir}/{build_name}/merge_nextclade_with_metadata.txt"
    shell:
        r"""
        augur merge \
           --metadata \
             metadata={input.metadata} \
             nextclade={input.nextclade} \
           --metadata-id-columns \
             metadata={params.metadata_id} \
             nextclade={params.nextclade_id} \
           --output-metadata {output.merged} 2>&1 | tee {log}
        """

def get_subsample_input(w):
    files = {"metadata": f"{build_dir}/{w.build_name}/metadata_with_nextclade.tsv"}
    if config['builds'][w.build_name]['subsamples'][w.subsample].get('priorities', '')=='titers':
        files['titers']=build_dir + f"/{w.build_name}/titer_priorities.tsv"
    return files

rule subsample:
    input:
        unpack(get_subsample_input)
    output:
        strains = build_dir + "/{build_name}/strains_{subsample}.txt",
        filter_log = build_dir + "/{build_name}/strains_{subsample}_filter_log.tsv",
    params:
        filters =  lambda w: config["builds"][w.build_name]["subsamples"][w.subsample]["filters"],
        priorities = lambda w: f"--priority {build_dir}/{w.build_name}/titer_priorities.tsv" \
                               if config['builds'][w.build_name]['subsamples'][w.subsample].get('priorities', '')=='titers' else ''
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/subsample_{build_name}_{subsample}.txt"
    log:
        "logs/subsample_{build_name}_{subsample}.txt"
    shell:
        r"""
        augur filter \
            --metadata {input.metadata} \
            {params.filters} \
            {params.priorities} \
            --output-strains {output.strains} \
            --output-log {output.filter_log} 2>&1 | tee {log}
        """

rule select_strains:
    input:
        metadata = build_dir + "/{build_name}/metadata_with_nextclade.tsv",
        subsamples = lambda w: [f"{build_dir}/{w.build_name}/strains_{s}.txt" for s in config['builds'][w.build_name]['subsamples']],
    output:
        metadata = build_dir + "/{build_name}/metadata.tsv",
        strains = build_dir + "/{build_name}/strains.txt",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/select_strains_{build_name}.txt"
    log:
        "logs/select_strains_{build_name}.txt"
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --exclude-all \
            --include {input.subsamples} \
            --output-metadata {output.metadata} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """

rule select_sequences:
    input:
        sequences = lambda w: f"data/{config['builds'][w.build_name]['lineage']}/{w.segment}.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv",
        strains = build_dir + "/{build_name}/strains.txt",
    output:
        sequences = build_dir + "/{build_name}/{segment}/sequences.fasta",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/select_sequences_{build_name}_{segment}.txt"
    log:
        "logs/select_sequences_{build_name}_{segment}.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude-all \
            --include {input.strains} \
            --output-sequences {output.sequences} 2>&1 | tee {log}
        """

def get_titer_collection_data(wildcards):
    return [
        collection["data"]
        for collection in config["builds"][wildcards.build_name]["titer_collections"]
        if collection["name"] == wildcards.titer_collection
    ][0]

rule select_titers:
    input:
        strains = build_dir + "/{build_name}/strains.txt",
        titers = get_titer_collection_data,
    output:
        titers = build_dir + "/{build_name}/titers/{titer_collection}.tsv",
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/select_titers_{build_name}_{titer_collection}.txt"
    log:
        "logs/select_titers_{build_name}_{titer_collection}.txt"
    shell:
        """
        head -n 1 {input.titers} > {output.titers};
        tsv-join \
            --key-fields 1 \
            --filter-file {input.strains} \
            {input.titers} >> {output.titers} 2> {log}
        """
