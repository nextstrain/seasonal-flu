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


rule select_strains:
    input:
        sequences = expand("data/{{lineage}}/{segment}.fasta", segment=segments),
        metadata = "data/{lineage}/metadata.tsv",
        titers = expand("data/{{lineage}}/{center}_{passage}_{assay}_titers.tsv",
                        center=config[wildcards.build_name]["center"],
                        passage=config[wildcards.build_name]["passage"],
                        assay=config[wildcards.build_name]["assay"])
        include = files.references
    output:
        strains = "builds/{build_name}/strains.txt",
    params:
        viruses_per_month = vpm
    conda: "environment.yaml"
    shell:
        """
        python3 scripts/select_strains.py \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --segments {segments} \
            --include {input.include} \
            --lineage {wildcards.lineage} \
            --resolution {wildcards.resolution} \
            --viruses-per-month {params.viruses_per_month} \
            --titers {input.titers} \
            --output {output.strains}
        """

rule select_sequences:
    input:
        strains = "builds/{build_name}/strains.txt",
        sequences = "data/{lineage}/{segment}.fasta"
    output:
        sequences = "builds/{build_name}/{segment}/sequences.fasta"
    shell:
        """
        seqkit grep -f {input.strains} {input.sequences} -o {output.sequences}
        """

rule select_metadata:
    input:
        strains = "builds/{build_name}/strains.txt",
        metadata = "data/{lineage}/metadata.tsv"
    output:
        metadata = "builds/{build_name}/metadata.tsv"
    run:
        import pandas as pd
        with open(input.strains) as fh:
            strains = [x.strip() for x in fh]

        d = pd.read_csv(input.metadata, sep='\t', index_col=0).loc[strains]

        d.to_csv(output.metadata, sep='\t')

rule select_titers:
    input:
        strains = "builds/{build_name}/strains.txt",
        titers = expand("data/{{lineage}}/{center}_{passage}_{assay}_titers.tsv",
                        center=config[wildcards.build_name]["center"],
                        passage=config[wildcards.build_name]["passage"],
                        assay=config[wildcards.build_name]["assay"])
    output:
        titers = "builds/{build_name}/titers.tsv"
    run:
        import pandas as pd
        with open(input.strains) as fh:
            strains = set([x.strip() for x in fh])

        d = pd.read_csv(input.titers, sep='\t')
        relevant_titers = []
        for ri, row in d.rowiter():
            if row[0] in strains and row[1] in strains:
                relevant_titers.append(row)

        pd.DataFrame(relevant_titers).to_csv(output.titers, sep='\t')



