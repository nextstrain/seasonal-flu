segments = ['ha', 'na']
lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam']
resolutions = ['6m', '2y'] #, '3y', '6y', '12y']

passages = ['cell']
centers = ['cdc']
assays = ['hi']

localrules: download_all, simplify_auspice_names, targets, clean, clobber
include: "Snakefile_base"

rule all_live:
    input:
        auspice_main = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json",
                              lineage=lineages, segment=segments, resolution=resolutions),
        auspice_tip_frequencies = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json",
                              lineage=lineages, segment=segments, resolution=resolutions)

# separate rule for interaction with fauna
rule download_all:
    input:
        titers = expand("data/{center}_{lineage}_{passage}_{assay}_titers.tsv",
                         center=centers, lineage=lineages, passage=passages, assay=assays),
        sequences = expand("data/{lineage}_{segment}.fasta", lineage=lineages, segment=segments)


def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model,
        rules.clades.output.clades,
        rules.traits.output.node_data,
        rules.lbi.output.lbi
    ]

    # Only request a distance file for builds that have distance map
    # configurations defined.
    if _get_build_distance_map_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = files.auspice_config,
        node_data = _get_node_data_for_export
    output:
        auspice = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_tree.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice} \
            --minify-json
        """

rule simplify_auspice_names:
    input:
        tree = "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_hi_tree.json",
        frequencies = "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_hi_tip-frequencies.json"
    output:
        tree = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json"
    shell:
        '''
        mv {input.tree} {output.tree} &
        mv {input.frequencies} {output.frequencies} &
        '''

rule targets:
    input:
        tree = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tree.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json"
    output:
        target = "targets/flu_seasonal_{lineage}_{segment}_{resolution}"
    shell:
        '''
        touch {output.target}
        '''

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "targets ",
        "auspice ",
        "auspice-who ",
        "logs"
    shell:
        "rm -rfv {params}"

rule clobber:
    message: "Removing directories: {params}"
    params:
        "results ",
        "targets ",
        "auspice ",
        "auspice-who ",
        "logs ",
        "data"
    shell:
        "rm -rfv {params}"
