'''
This file contains rules infer titer models

input:
 - builds/{build_name}/metadata.tsv
 - builds/{build_name}/titers.tsv
 - builds/{build_name}/{segment}/tree.nwk

output:

'''
build_dir = config.get("build_dir", "builds")

rule titers_sub:
    input:
        titers = build_dir +"/{build_name}/titers.tsv",
        tree = rules.refine.output.tree,
        translations_done = build_dir + "/{build_name}/{segment}/translations.done"
    params:
        translations = lambda w: [f"{build_dir}/{w.build_name}/{w.segment}/nextalign/sequences.gene.{gene}_withInternalNodes.fasta" for gene in GENES[w.segment]],
        genes = lambda w: GENES[w.segment]
    output:
        titers_model = build_dir + "/{build_name}/{segment}/titers-sub-model.json",
    conda: "environment.yaml"
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {params.translations} \
            --gene-names {params.genes} \
            --tree {input.tree} \
            --allow-empty-model \
            --output {output.titers_model}
        """

rule titers_tree:
    input:
        titers = build_dir +"/{build_name}/titers.tsv",
        tree = rules.refine.output.tree
    output:
        titers_model = build_dir + "/{build_name}/{segment}/titers-tree-model.json",
    conda: "environment.yaml"
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --allow-empty-model \
            --output {output.titers_model}
        """
