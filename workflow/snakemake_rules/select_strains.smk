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

localrules: titer_priorities, select_strains, select_metadata, select_titers

build_dir = config.get("build_dir", "builds")


def get_titers_for_build(w):
    return "data/{lineage}/{center}_{passage}_{assay}_titers.tsv".format(**config['builds'][w.build_name])

rule titer_priorities:
    input:
        titers = get_titers_for_build,
        metadata = lambda w: f"data/{config['builds'][w.build_name]['lineage']}/metadata.tsv"
    output:
        priorities = build_dir + "/{build_name}/titer_priorities.tsv"
    run:
        import pandas as pd
        metadata = pd.read_csv(input.metadata, sep='\t')
        titer_counts = {s:0 for s in metadata.strain}

        try:
            titers = pd.read_csv(input.titers, sep='\t')
            for ri, row in titers.iterrows():
                test_strain = row.iloc[0]
                if test_strain in titer_counts:
                    titer_counts[test_strain] += 1
        except:
            pass

        with open(output.priorities, 'w') as fh:
            for s,p in titer_counts.items():
                fh.write(f"{s}\t{p}\n")

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


rule select_titers:
    input:
        strains = build_dir + "/{build_name}/strains.txt",
        titers = get_titers_for_build
    output:
        titers = build_dir + "/{build_name}/titers.tsv"
    run:
        with open(input.strains) as fh:
            strains = set([x.strip() for x in fh.readlines()])

        with open(input.titers) as fh:
            with open(output.titers, 'w') as out_fh:
                for line in fh:
                    if line.split('\t')[0] in strains:
                        out_fh.write(line)




