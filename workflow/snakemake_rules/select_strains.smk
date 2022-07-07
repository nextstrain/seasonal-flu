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
        fasta_fields=config["fasta_fields"],
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
        metadata="data/{lineage}/metadata.tsv",
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

def get_titers_for_build(w):
    return "data/{lineage}/{center}_{passage}_{assay}_titers.tsv".format(**config['builds'][w.build_name])

rule titer_priorities:
    input:
        titers = get_titers_for_build,
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

def get_subsample_input(w):
    files = {"metadata": f"data/{config['builds'][w.build_name]['lineage']}/metadata.tsv"}
    if config['builds'][w.build_name]['subsamples'][w.subsample].get('priorities', '')=='titers':
        files['titers']=build_dir + f"/{w.build_name}/titer_priorities.tsv"
    return files

rule subsample:
    input:
        unpack(get_subsample_input)
    output:
        strains = build_dir + "/{build_name}/strains_{subsample}.txt",
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
        """
        augur filter \
            --metadata {input.metadata} \
            {params.filters} \
            {params.priorities} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """

rule select_strains:
    input:
        metadata = lambda w: f"data/{config['builds'][w.build_name]['lineage']}/metadata.tsv",
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

rule select_titers:
    input:
        strains = build_dir + "/{build_name}/strains.txt",
        titers = get_titers_for_build
    output:
        titers = build_dir + "/{build_name}/titers.tsv"
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/select_titers_{build_name}.txt"
    log:
        "logs/select_titers_{build_name}.txt"
    shell:
        """
        tsv-join \
            --key-fields 1 \
            --filter-file {input.strains} \
            {input.titers} > {output.titers} 2> {log}
        """
