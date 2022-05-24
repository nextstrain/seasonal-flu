
build_dir = config.get("build_dir", "builds")
localrules: glyc, lbi

glyc_gene = {'ha':'HA1', 'na':'NA'}

rule glyc:
    input:
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    output:
        glyc = build_dir + "/{build_name}/{segment}/glyc.json"
    params:
        alignment = lambda w: f"{build_dir}/{w.build_name}/{w.segment}/nextalign/sequences.gene.{glyc_gene.get(w.segment)}_withInternalNodes.fasta",
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/glyc.py \
            --tree {input.tree} \
            --alignment {params.alignment} \
            --output {output.glyc}
        """

rule lbi:
    message: "Calculating LBI"
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = 0.5,
        window = 0.5,
        names = "lbi"
    output:
        lbi =  build_dir + "/{build_name}/{segment}/lbi.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --output {output} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window}
        """
