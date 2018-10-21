path_to_fauna = '../fauna/data'
segments = ['ha', 'na']

def reference_strain(v):
    references = {'h3n2':"A/Beijing/32/1992",
                  'h1n1pdm':"A/California/07/2009",
                  'vic':"B/HongKong/02/1993",
                  'yam':"B/Singapore/11/1994"}
    return references[v.lineage]



def vpm(v):
    vpm = {'3y':15, '6y':12, '12y':5}
    return vpm[v.resolution] if v.resolution in vpm else 5


rule all:
    input:
        auspice_tree = "auspice/flu_{build}_{lineage}_{segment}_{resolution}_tree.json",
        auspice_meta = "auspice/flu_{build}_{lineage}_{segment}_{resolution}_meta.json"

rule files:
    params:
        input_fasta = path_to_fauna+"/{lineage}_{segment}.fasta",
        # dropped_strains = "config/dropped_strains_{lineage}.txt",
        reference = "config/outgroup_{lineage}_{segment}.gb"

files = rules.files.params


rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = files.input_fasta
    output:
        sequences = "results/sequences_{lineage}_{segment}.fasta",
        metadata = "results/metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields =  "strain virus isolate_id date region country division passage authors age gender"

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
        metadata = lambda w:expand("results/metadata_{lineage}_{segment}.tsv", segment=segments, lineage=w.lineage)
    output:
        strains = "results/strains_{build}_{lineage}_{resolution}.txt",
    params:
        virus_per_month = vpm
    shell:
        """
        python scripts/prepare.py --metadata {input.metadata} \
                                  --segments {segments} \
                                  --resolution {wildcards.resolution} --lineage {wildcards.lineage} \
                                  --output {output.strains}
        """


rule filter:
    input:
        metadata = rules.parse.output.metadata,
        sequences = 'results/sequences_{lineage}_{segment}.fasta',
        strains = rules.select_strains.output.strains
    output:
        sequences = 'results/sequences_{build}_{lineage}_{segment}_{resolution}.fasta'
    run:
        from Bio import SeqIO
        with open(input.strains) as infile:
            strains = set(map(lambda x:x.strip(), infile.readlines()))
        with open(output.sequences, 'w') as outfile:
            for seq in SeqIO.parse(input.sequences, 'fasta'):
                if seq.name in strains:
                    SeqIO.write(seq, outfile, 'fasta')

