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

build_dir = config.get("build_dir", "builds")


rule subsample:
    input:
        metadata = lambda w: f"data/{config['builds'][w.build_name]['lineage']}/metadata.tsv"
    output:
        strains = build_dir + "/{build_name}/strains_{subsample}.txt",
    params:
        filters =  lambda w: config["builds"][w.build_name]["subsamples"][w.subsample]["filters"]
    conda: "environment.yaml"
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            {params.filters} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """

rule select_strains:
    input:
        metadata = lambda w: f"data/{config['builds'][w.build_name]['lineage']}/metadata.tsv",
        subsamples = lambda w: [f"{build_dir}/{w.build_name}/strains_{s}.txt" for s in config['builds'][w.build_name]['subsamples']]
    output:
        strains = build_dir + "/{build_name}/strains.txt",
    run:
        strains = set()
        for fname in input.subsamples:
            with open(fname) as fh:
                for line in fh:
                    strains.add(line.strip())

        with open(output.strains, 'w') as fh:
            fh.write('\n'.join(strains)+'\n')

rule select_sequences:
    input:
        strains = build_dir + "/{build_name}/strains.txt",
        sequences = lambda w: f"data/{config['builds'][w.build_name]['lineage']}/{w.segment}.fasta"
    output:
        sequences = build_dir + "/{build_name}/{segment}/sequences.fasta"
    shell:
        """
        seqkit grep -f {input.strains} {input.sequences} | seqkit replace -p "[^acgtACGT-]|N" -s -o {output.sequences}
        """

rule select_metadata:
    input:
        strains = build_dir + "/{build_name}/strains.txt",
        metadata = lambda w: f"data/{config['builds'][w.build_name]['lineage']}/metadata.tsv"
    output:
        metadata = build_dir + "/{build_name}/metadata.tsv"
    run:
        import pandas as pd
        with open(input.strains) as fh:
            strains = [x.strip() for x in fh]

        d = pd.read_csv(input.metadata, sep='\t', index_col=0).loc[strains]

        d.to_csv(output.metadata, sep='\t')



