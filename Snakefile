from datetime import date
from treetime.utils import numeric_date

path_to_fauna = '../fauna'
if os.environ.get('FAUNA_PATH'):
    path_to_fauna = os.environ.get('FAUNA_PATH')

min_length = 800
segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns'] # ordering is used by scripts/reassort
lineages = ['h3n2', 'h1n1pdm']
resolutions = ['2y']

def reference_strain(v):
    references = {'h3n2':"A/Beijing/32/1992",
                  'h1n1pdm':"A/California/07/2009",
                  'vic':"B/HongKong/02/1993",
                  'yam':"B/Singapore/11/1994"
                  }
    return references[v.lineage]

def gene_names(w):
    genes_to_translate = {'ha':['SigPep', 'HA1', 'HA2'], 'na':['NA']}
    return genes_to_translate[w.segment]

def translations(w):
    genes = gene_names(w)
    return ["results/aa-seq_%s_%s_%s_%s.fasta"%(w.lineage, w.segment, w.resolution, g)
            for g in genes]

def clock_rate(w):
    rate = {
        ('h3n2', 'ha'): 0.00356, ('h3n2', 'na'): 0.00298, ('h3n2', 'mp'): 0.00082,
        ('h3n2', 'np'): 0.00140, ('h3n2', 'ns'): 0.00193, ('h3n2', 'pa'): 0.00226,
        ('h3n2', 'pb1'): 0.00177, ('h3n2', 'pb2'): 0.00230,
        ('h1n1pdm', 'ha'): 0.00329, ('h1n1pdm', 'na'): 0.00342, ('h1n1pdm', 'mp'): 0.00209,
        ('h1n1pdm', 'np'): 0.00196, ('h1n1pdm', 'ns'): 0.00278, ('h1n1pdm', 'pa'): 0.00235,
        ('h1n1pdm', 'pb1'): 0.00188, ('h1n1pdm', 'pb2'): 0.00224,
        ('vic', 'ha'): 0.0024, ('vic', 'na'): 0.0015,
        ('yam', 'ha'): 0.0019, ('yam', 'na'): 0.0013
    }
    return rate[(w.lineage, w.segment)]

#
# Define clades functions
#
def _get_clades_file_for_wildcards(wildcards):
    if wildcards.segment == "ha":
        return "config/clades_%s_ha.tsv"%(wildcards.lineage)
    else:
        return "results/clades_%s_ha_%s.json"%(wildcards.lineage, wildcards.resolution)

#
# Define rules.
#

rule all:
    input:
        auspice_tree = expand("auspice/seattle_flu_seasonal_{lineage}_{segment}_{resolution}_global_tree.json", lineage=lineages, segment=segments, resolution=resolutions),
        auspice_meta = expand("auspice/seattle_flu_seasonal_{lineage}_{segment}_{resolution}_global_meta.json", lineage=lineages, segment=segments, resolution=resolutions)

rule files:
    params:
        seattle_metadata = "metadata/seattle_metadata.tsv",
        outliers = "config/outliers_{lineage}.txt",
        references = "config/references_{lineage}.txt",
        reference = "config/reference_{lineage}_{segment}.gb",
        colors = "config/colors.tsv",
        lat_longs = "config/lat_longs.tsv",
        auspice_config = "config/auspice_config_{lineage}.json",

files = rules.files.params

rule download_background_seqmeta:
    message: "Downloading background sequences from fauna"
    output:
        seqmeta = "data/background_seqmeta_{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus flu \
            --fasta_fields {params.fasta_fields} \
            --resolve_method split_passage \
            --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
            --path data \
            --fstem background_seqmeta_{wildcards.lineage}_{wildcards.segment}
        """

rule parse_background_seqmeta:
    message: "Parsing fasta into sequences and metadata"
    input:
        seqmeta = rules.download_background_seqmeta.output.seqmeta
    output:
        sequences = "data/background_sequences_{lineage}_{segment}.fasta",
        metadata = "data/background_metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields = "strain virus isolate_id date region country division location passage authors age sex"
    shell:
        """
        augur parse \
            --sequences {input.seqmeta} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule download_seattle_sequences:
    message: "Downloading Seattle sequences from fauna"
    output:
        sequences = "data/seattle_sequences_{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus seattle \
            --fasta_fields {params.fasta_fields} \
            --select segment:{wildcards.segment} type:{wildcards.lineage} \
            --path data \
            --fstem seattle_sequences_{wildcards.lineage}_{wildcards.segment}
        """

rule concat_sequences:
    message: "Concatenating background and Seattle sequences"
    input:
        background_sequences = rules.parse_background_seqmeta.output.sequences,
        seattle_sequences = rules.download_seattle_sequences.output.sequences
    output:
        sequences = "data/sequences_{lineage}_{segment}.fasta"
    shell:
        """
        cat {input.background_sequences} {input.seattle_sequences} > {output.sequences}
        """

rule concat_metadata:
    message: "Concatenating background and Seattle metadata"
    input:
        background_metadata = rules.parse_background_seqmeta.output.metadata,
        seattle_metadata = files.seattle_metadata
    output:
        metadata = "data/metadata_{lineage}_{segment}.tsv"
    shell:
        """
        python3 scripts/concat_metadata.py \
            --files {input.background_metadata} {input.seattle_metadata} \
            --mergeby strain \
            --fields date region site_type flu_shot gender \
            > {output.metadata}
        """

rule filter:
    message:
        """
        Filtering {wildcards.lineage} {wildcards.segment} sequences:
          - less than {params.min_length} bases
          - outliers
          - samples with missing region and country metadata
        """
    input:
        sequences = rules.concat_sequences.output.sequences,
        metadata = rules.concat_metadata.output.metadata,
        exclude = files.outliers
    output:
        sequences = 'results/filtered_{lineage}_{segment}.fasta'
    params:
        min_length = min_length
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --non-nucleotide \
            --exclude {input.exclude} \
            --exclude-where country=? region=? passage=egg \
            --output {output}
        """

rule select_strains:
    message:
        """
        Selecting Strains (scripts/select_strains.py)
        This automatically includes all strains with region=seattle
        """
    input:
        sequences = expand("results/filtered_{{lineage}}_{segment}.fasta", segment=segments),
        metadata = expand("data/metadata_{{lineage}}_{segment}.tsv", segment=segments),
        include = files.references
    output:
        strains = "results/strains_{lineage}_{resolution}.txt",
    params:
        viruses_per_month = 90
    shell:
        """
        python3 scripts/select_strains.py \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --segments {segments} \
            --include {input.include} \
            --lineage {wildcards.lineage} \
            --resolution {wildcards.resolution} \
            --viruses_per_month {params.viruses_per_month} \
            --output {output.strains}
        """

rule extract:
    message:
        """
        Extract sequences from a given FASTA file that match the given list of sample names
        """
    input:
        sequences = rules.filter.output.sequences,
        strains = rules.select_strains.output.strains
    output:
        sequences = 'results/extracted_{lineage}_{segment}_{resolution}.fasta'
    shell:
        """
        python3 scripts/extract_sequences.py \
            --sequences {input.sequences} \
            --samples {input.strains} \
            --output {output}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.extract.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{lineage}_{segment}_{resolution}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference \
            --nthreads 1
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree-raw_{lineage}_{segment}_{resolution}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads 1
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.concat_metadata.output.metadata
    output:
        tree = "results/tree_{lineage}_{segment}_{resolution}.nwk",
        node_data = "results/branch-lengths_{lineage}_{segment}_{resolution}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = clock_rate
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-rate {params.clock_rate} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt-muts_{lineage}_{segment}_{resolution}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa-muts_{lineage}_{segment}_{resolution}.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule reconstruct_translations:
    message: "Reconstructing translations required for titer models and frequencies"
    input:
        tree = rules.refine.output.tree,
        node_data = "results/aa-muts_{lineage}_{segment}_{resolution}.json",
    output:
        aa_alignment = "results/aa-seq_{lineage}_{segment}_{resolution}_{gene}.fasta"
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene {wildcards.gene} \
            --output {output.aa_alignment} \
            --internal-nodes
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.concat_metadata.output.metadata
    output:
        node_data = "results/traits_{lineage}_{segment}_{resolution}.json",
    params:
        columns = "region"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule clades:
    message: "Annotating clades"
    input:
        tree = "results/tree_{lineage}_ha_{resolution}.nwk",
        nt_muts = rules.ancestral.output,
        aa_muts = rules.translate.output,
        clades = _get_clades_file_for_wildcards
    output:
        node_data = "results/clades_{lineage}_{segment}_{resolution}.json"
    run:
        if wildcards.segment == 'ha':
            shell("""
                augur clades \
                    --tree {input.tree} \
                    --mutations {input.nt_muts} {input.aa_muts} \
                    --clades {input.clades} \
                    --output {output.node_data}
            """)
        else:
            shell("""
                python3 scripts/import_tip_clades.py \
                    --tree {input.tree} \
                    --clades {input.clades} \
                    --output {output.node_data}
            """)

rule lbi:
    message: "Calculating LBI"
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = 0.3,
        window = 0.5,
        names = "lbi"
    output:
        node_data = "results/lbi_{lineage}_{segment}_{resolution}.json"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window} \
            --output {output.node_data}
        """

rule hamming_distance:
    message:
        """
        Calculating hamming distance for {wildcards.lineage} {wildcards.segment} {wildcards.resolution}
        """
    input:
        aligned = rules.align.output.alignment
    output:
        data = "results/distance_{lineage}_{segment}_{resolution}.json"
    shell:
        """
        python3 scripts/reassort/hamming.py --fasta {input.aligned} --output {output.data}
        """

rule clustering:
    message: "Identifying clusters based on connected components"
    input:
        nt_muts = expand("results/nt-muts_{{lineage}}_{segment}_{{resolution}}.json", segment=segments),
    params:
        cutoff = 30
    output:
        node_data = "results/clustering_{lineage}_{resolution}.json"
    shell:
        """
        python3 scripts/connected_components.py \
            --nt-muts {input.nt_muts} \
            --cutoff {params.cutoff} \
            --output {output.node_data}
        """

# def _get_trees_for_all_segments(wildcards):
#     trees = []
#     for seg in segments:
#         trees.append(rules.refine.output.tree.format(**wildcards, **{"segment": seg}))
#     return trees
#
# rule identify_non_reassorting_tips:
#     message: "Identifying sets of tips which have not reassorted"
#     input:
#         trees = _get_trees_for_all_segments,
#         mutations = lambda wildcards: [rules.ancestral.output.node_data.format(**wildcards, **{"segment": seg}) for seg in segments]
#     output:
#         data = "results/reassort_{lineage}_{resolution}.json"
#     shell:
#         """
#         python3 scripts/reassort --trees {input.trees} --mutations {input.mutations} --output {output.data}
#         """

def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.clades.output.node_data,
        rules.traits.output.node_data,
        rules.lbi.output.node_data,
        rules.clustering.output.node_data
    ]

    # HA gets the reassortant information
    # if wildcards["segment"] == "ha":
    #     inputs.append(rules.identify_non_reassorting_tips.output.data)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.concat_metadata.output.metadata,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config,
        node_data = _get_node_data_for_export
    output:
        auspice_tree = "auspice/seattle_flu_seasonal_{lineage}_{segment}_{resolution}_global_tree.json",
        auspice_meta = "auspice/seattle_flu_seasonal_{lineage}_{segment}_{resolution}_global_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
