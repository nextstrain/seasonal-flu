from datetime import date
from treetime.utils import numeric_date

path_to_fauna = '../fauna'
segments = ['ha', 'na']
lineages = ['h3n2']
resolutions = ['2y', '3y', '6y', '12y']
frequency_regions = ['north_america', 'south_america', 'europe', 'china',
                     'southeast_asia', 'japan_korea', 'south_asia', 'africa']


def reference_strain(v):
    references = {'h3n2':"A/Beijing/32/1992",
                  'h1n1pdm':"A/California/07/2009",
                  'vic':"B/HongKong/02/1993",
                  'yam':"B/Singapore/11/1994"
                  }
    return references[v.lineage]

genes_to_translate = {'ha':['HA1', 'HA2'], 'na':['NA']}
def gene_names(w):
    return genes_to_translate[w.segment]

def translations(w):
    genes = gene_names(w)
    return ["results/aa-seq_%s_%s_%s_%s.fasta"%(w.lineage, w.segment, w.resolution, g)
            for g in genes]

def pivots_per_year(w):
    pivots_per_year = {'2y':12, '3y':6, '6y':4, '12y':2}
    return pivots_per_year[w.resolution]

def min_date(w):
    now = numeric_date(date.today())
    return now - int(w.resolution[:-1])

def max_date(w):
    return numeric_date(date.today())

def substitution_rates(w):
    references = {('h3n2', 'ha'): 0.0038, ('h3n2', 'na'):0.0028}
    return references[(w.lineage, w.segment)]

def vpm(v):
    vpm = {'2y':2, '3y':2, '6y':2, '12y':1}
    return vpm[v.resolution] if v.resolution in vpm else 5


rule all:
    input:
        auspice_tree = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tree.json", lineage=lineages, segment=segments, resolution=resolutions),
        auspice_meta = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_meta.json", lineage=lineages, segment=segments, resolution=resolutions)

rule files:
    params:
        outliers = "config/outliers_{lineage}.txt",
        references = "config/references_{lineage}.txt",
        reference = "config/{lineage}_{segment}_outgroup.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json",

files = rules.files.params

rule download_sequences:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    shell:
        """
        env PYTHONPATH={path_to_fauna} \
            python2 {path_to_fauna}/vdb/download.py \
                --database vdb \
                --virus flu \
                --fasta_fields {params.fasta_fields} \
                --resolve_method split_passage \
                --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
                --path data \
                --fstem {wildcards.lineage}_{wildcards.segment}
        """

rule download_titers:
    message: "Downloading titers from fauna"
    output:
        titers = "data/{lineage}_hi_titers.tsv"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    shell:
        """
        env PYTHONPATH={path_to_fauna} \
            python2 {path_to_fauna}/tdb/download.py \
                --database cdc_tdb \
                --virus flu \
                --subtype {wildcards.lineage} \
                --select assay_type:hi \
                --path data \
                --fstem {wildcards.lineage}_hi
        """

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.download_sequences.output.sequences
    output:
        sequences = "results/sequences_{lineage}_{segment}.fasta",
        metadata = "results/metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields =  "strain virus isolate_id date region country division location passage authors age gender"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule select_strains:
    input:
        metadata = lambda w:expand("results/metadata_{lineage}_{segment}.tsv", segment=segments, lineage=w.lineage),
        titers = rules.download_titers.output.titers
    output:
        strains = "results/strains_{lineage}_{resolution}.txt",
    params:
        viruses_per_month = vpm,
        exclude = files.outliers,
        include = files.references
    shell:
        """
        python scripts/select_strains.py \
            --metadata {input.metadata} \
            --segments {segments} \
            --exclude {params.exclude} \
            --include {params.include} \
            --lineage {wildcards.lineage} \
            --resolution {wildcards.resolution} \
            --viruses_per_month {params.viruses_per_month} \
            --titers {input.titers} \
            --output {output.strains}
        """

rule filter:
    input:
        metadata = rules.parse.output.metadata,
        sequences = 'results/sequences_{lineage}_{segment}.fasta',
        strains = rules.select_strains.output.strains
    output:
        sequences = 'results/filtered_{lineage}_{segment}_{resolution}.fasta'
    run:
        from Bio import SeqIO
        with open(input.strains) as infile:
            strains = set(map(lambda x:x.strip(), infile.readlines()))
        with open(output.sequences, 'w') as outfile:
            for seq in SeqIO.parse(input.sequences, 'fasta'):
                if seq.name in strains:
                    SeqIO.write(seq, outfile, 'fasta')

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
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
            --remove-reference
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
            --output {output.tree}
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
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree_{lineage}_{segment}_{resolution}.nwk",
        node_data = "results/branch-lengths_{lineage}_{segment}_{resolution}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
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
    params:
        gene = "{gene}",
    output:
        aa_alignment = "results/aa-seq_{lineage}_{segment}_{resolution}_{gene}.fasta"
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene {params.gene} \
            --output {output.aa_alignment}
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
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

rule titers_sub:
    input:
        titers = rules.download_titers.output.titers,
        aa_muts = rules.translate.output,
        alignments = translations
    params:
        genes = gene_names
    output:
        titers_model = "results/titers-sub-model_{lineage}_{segment}_{resolution}.json",
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --output {output.titers_model}
        """

rule titers_tree:
    input:
        titers = rules.download_titers.output.titers,
        tree = rules.refine.output.tree
    output:
        titers_model = "results/titers-tree-model_{lineage}_{segment}_{resolution}.json",
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --output {output.titers_model}
        """

rule mutation_frequencies:
    input:
        metadata = rules.parse.output.metadata,
        alignment = translations
    params:
        genes = gene_names,
        min_date = min_date,
        max_date = max_date,
        pivots_per_year = pivots_per_year
    output:
        mut_freq = "results/mutation-frequencies_{lineage}_{segment}_{resolution}.json"
    shell:
        """
        augur frequencies \
            --alignments {input.alignment} \
            --metadata {input.metadata} \
            --gene-names {params.genes} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --pivots-per-year {params.pivots_per_year} \
            --output {output.mut_freq}
        """

rule tree_frequencies:
    input:
        metadata = rules.parse.output.metadata,
        tree = rules.refine.output.tree
    params:
        regions = frequency_regions + ['global'],
        min_date = min_date,
        max_date = max_date,
        pivots_per_year = pivots_per_year
    output:
        tree_freq = "results/tree-frequencies_{lineage}_{segment}_{resolution}.json",
    shell:
        """
        augur frequencies \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --pivots-per-year {params.pivots_per_year} \
            --regions {params.regions} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output.tree_freq}
        """

rule clades:
    message: "Annotating clades"
    input:
        tree = rules.refine.output.tree,
        nt_muts = rules.ancestral.output,
        aa_muts = rules.translate.output,
        clade_definitions = "config/clades_{lineage}_{segment}.tsv"
    output:
        clades = "results/clades_{lineage}_{segment}_{resolution}.json"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nt_muts} {input.aa_muts} \
            --clades {input.clade_definitions} \
            --output {output.clades}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        node_data = rules.refine.output.node_data,
        metadata = rules.parse.output.metadata,
        nt_muts = rules.ancestral.output,
        aa_muts = rules.translate.output,
        titers_tree_model = rules.titers_tree.output.titers_model,
        clades = rules.clades.output.clades,
        traits = rules.traits.output.node_data,
        auspice_config = files.auspice_config
    output:
        auspice_tree = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tree.json",
        auspice_meta = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} {input.nt_muts} {input.aa_muts} {input.titers_tree_model} {input.clades} {input.traits} \
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
