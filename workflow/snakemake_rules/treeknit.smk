ruleorder: refine_ARG>refine


rule treeknit:
    input:
        trees = rules.sanitize_trees.output.trees,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        trees = expand("{build_dir}/{{build_name}}/TreeKnit/tree_common_{segment}.resolved.nwk",
                       segment=config['segments'][:2], build_dir=[build_dir]),
        mccs = build_dir + "/{build_name}/TreeKnit/MCCs.dat"
    params:
        treetime_tmpdir = build_dir + "/{build_name}/TreeTime_tmp",
        tmp_trees = expand("{build_dir}/{{build_name}}/TreeKnit/tree_{segment}.nwk",
                            segment=config['segments'][:2], build_dir=[build_dir]),
        treeknit_tmpdir = build_dir + "/{build_name}/TreeKnit",
        clock_filter=4
    shell:
        """
        treeknit {input.trees[0]} {input.trees[1]} --outdir {params.treeknit_tmpdir}
        """

rule treetime_arg:
    message:
        """
        Refining tree
          - estimate timetree
        """
    input:
        trees = rules.treeknit.output.trees,
        mccs = rules.treeknit.output.mccs,
        alignments = expand("{build_dir}/{{build_name}}/{segment}/aligned.fasta",
                            segment=config['segments'][:2], build_dir=[build_dir]),
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        directory(build_dir + "/{build_name}/treetime_arg")
    params:
        coalescent = "const",
        date_inference = "always"
    resources:
        mem_mb=16000
    shell:
        """
        treetime arg --trees {input.trees} \
                      --mccs {input.mccs} \
                      --alignments {input.alignments} \
                      --confidence --clock-std-dev 0.001 \
                      --clock-filter 0 \
                      --outdir  {output} \
                      --time-marginal {params.date_inference} \
                      --dates {input.metadata} --keep-root --keep-polytomies
        """


rule refine_ARG:
    message:
        """
        Refining tree
          - estimate timetree
        """
    input:
        treetime_arg = rules.treetime_arg.output,
        mccs = rules.treeknit.output.mccs,
        metadata = build_dir + "/{build_name}/metadata.tsv",
    output:
        trees = expand("{build_dir}/{{build_name}}/{segment}/tree.nwk",
                        segment=config['segments'][:2], build_dir=[build_dir]),
        node_data = expand("{build_dir}/{{build_name}}/{segment}/branch-lengths.json",
                            segment=config['segments'][:2], build_dir=[build_dir]),
    conda: "../envs/nextstrain.yaml"
    benchmark:
        "benchmarks/refine_{build_name}.txt"
    log:
        "logs/refine_{build_name}.txt"
    resources:
        mem_mb=16000
    shell:
        """
        python scripts/make-branch-length-json.py --timetree {input.treetime_arg}/timetree_1.nexus \
                --divtree {input.treetime_arg}/divergence_tree_1.nexus \
                --molecular-clock {input.treetime_arg}/molecular_clock_1.txt \
                --dates {input.treetime_arg}/dates_1.tsv --mccs {input.mccs} \
                --output-tree {output.trees[0]} --output-node-data {output.node_data[0]}

        python scripts/make-branch-length-json.py --timetree {input.treetime_arg}/timetree_2.nexus \
                --divtree {input.treetime_arg}/divergence_tree_2.nexus \
                --molecular-clock {input.treetime_arg}/molecular_clock_2.txt \
                --dates {input.treetime_arg}/dates_2.tsv  --mccs {input.mccs}  \
                --output-tree {output.trees[1]} --output-node-data {output.node_data[1]}
        """

