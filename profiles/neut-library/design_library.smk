rule download_mlr_json:
    output:
        mlr="builds/{build_name}/ha/mlr/model.json",
    params:
        model_url=lambda wildcards: f"https://data.nextstrain.org/files/workflows/forecasts-flu/gisaid/aa_haplotype/{config['builds'][wildcards.build_name]['lineage']}/region/mlr/MLR_results.json",
    shell:
        r"""
        curl --compressed -o {output.mlr:q} {params.model_url:q}
        """

rule parse_frequencies_and_ga_from_mlr_json:
    input:
        mlr="builds/{build_name}/ha/mlr/model.json",
    output:
        fitnesses="builds/{build_name}/ha/mlr/fitnesses.tsv",
    shell:
        r"""
        python scripts/parse-json.py \
            --input {input.mlr} \
            --outga {output.fitnesses}
        """

rule merge_metadata_and_nextclade_for_library_design:
    input:
        metadata="builds/{build_name}/metadata.tsv",
        nextclade="data/{build_name}/ha/nextclade.tsv",
    output:
        metadata="builds/{build_name}/metadata_with_nextclade.tsv",
    shell:
        r"""
        augur merge \
            --metadata metadata={input.metadata} \
                       nextclade={input.nextclade} \
            --metadata-id-columns metadata=strain \
                                  nextclade=seqName \
            --output-metadata {output.metadata}
        """

rule annotate_metadata_with_library_haplotypes:
    input:
        metadata="builds/{build_name}/metadata_with_nextclade.tsv",
        distance_maps=lambda wildcards: config["builds"][wildcards.build_name]["distance_maps"],
        recurrent_substitutions_map=lambda wildcards: config["builds"][wildcards.build_name]["recurrent_substitutions_map"],
        translations_dir=directory(build_dir + "/{build_name}/ha/translations"),
    output:
        library_metadata="builds/{build_name}/library_metadata.tsv",
    params:
        genes=GENES["ha"],
        alignments=[f"builds/{{build_name}}/ha/translations/{gene}.fasta" for gene in GENES["ha"]],
    shell:
        r"""
        python scripts/prepare_library_haplotypes.py \
            --metadata {input.metadata} \
            --distance-maps {input.distance_maps:q} \
            --recurrent-substitutions-map {input.recurrent_substitutions_map:q} \
            --alignments {params.alignments:q} \
            --gene-names {params.genes:q} \
            --output {output.library_metadata}
        """

rule design_library:
    input:
        metadata="builds/{build_name}/library_metadata.tsv",
        fitnesses="builds/{build_name}/ha/mlr/fitnesses.tsv",
    output:
        library_design="builds/{build_name}/library_design.tsv",
    params:
        mutation_weights=lambda wildcards: " ".join([
            f"{mutation_name}={mutation_weight}"
            for mutation_name, mutation_weight in config["builds"][wildcards.build_name]["mutation_weights"].items()
        ]),
    shell:
        r"""
        python scripts/summarize_haplotypes_for_library_design.py \
            --metadata {input.metadata} \
            --fitnesses {input.fitnesses} \
            --mutation-weights {params.mutation_weights} \
            --output {output.library_design}
        """

rule design_libraries:
    input:
        expand("builds/{build_name}/library_design.tsv", build_name=list(config["builds"].keys()))
