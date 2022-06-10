'''
This file contains rules that select strains for a build, extracts the sequences
and the meta data subset.

input:
 - "data/lineage/metadata.tsv"
 - "data/lineage/{segment}.fasta"
 - "data/{lineage}/{center}_{passage}_{assay}_titers.tsv"

output:
 - builds/{build_name}/strains.txt
 - builds/{build_name}/metadata.tsv
 - builds/{build_name}/titers.tsv
 - builds/{build_name}/{segment}/sequences.fasta
'''

localrules: titer_priorities, select_titers

build_dir = config.get("build_dir", "builds")


def get_titers_for_build(w):
    return "data/{lineage}/{center}_{passage}_{assay}_titers.tsv".format(**config['builds'][w.build_name])

rule titer_priorities:
    input:
        titers = get_titers_for_build,
    output:
        priorities = build_dir + "/{build_name}/titer_priorities.tsv",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        tsv-summarize -H --group-by virus_strain --count {input.titers} | sed 1d > {output.priorities}
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
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --exclude-all \
            --include {input.subsamples} \
            --output-metadata {output.metadata} \
            --output-strains {output.strains}
        """

rule select_sequences:
    input:
        sequences = lambda w: f"data/{config['builds'][w.build_name]['lineage']}/{w.segment}.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv",
        strains = build_dir + "/{build_name}/strains.txt",
    output:
        sequences = build_dir + "/{build_name}/{segment}/sequences.fasta",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude-all \
            --include {input.strains} \
            --output-sequences {output.sequences}
        """

rule select_titers:
    input:
        strains = build_dir + "/{build_name}/strains.txt",
        titers = get_titers_for_build
    output:
        titers = build_dir + "/{build_name}/titers.tsv"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        tsv-join \
            --key-fields 1 \
            --filter-file {input.strains} \
            {input.titers} > {output.titers}
        """
