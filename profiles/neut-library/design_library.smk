rule merge_metadata_and_nextclade_for_library_design:
    input:
        metadata="builds/{build_name}/metadata.tsv",
        nextclade=lambda wildcards: f"data/{config['builds'][wildcards.build_name].get('lineage')}/ha/nextclade.tsv",
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
        lcr_weights=lambda wildcards: f"config/{config['builds'][wildcards.build_name]['lineage']}/ha/lcr_weights.tsv",
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
            --lcr-weights {input.lcr_weights:q} \
            --alignments {params.alignments:q} \
            --gene-names {params.genes:q} \
            --output {output.library_metadata}
        """

rule design_library:
    input:
        metadata="builds/{build_name}/library_metadata.tsv",
        emerging_haplotype_fitnesses="builds/{build_name}/ha/mlr/emerging_haplotype/fitnesses.tsv",
        aa_haplotype_fitnesses="builds/{build_name}/ha/mlr/aa_haplotype/fitnesses.tsv",
    output:
        library_design="builds/{build_name}/library_design.tsv",
    params:
        mutation_weights=lambda wildcards: " ".join([
            f"{mutation_name}={mutation_weight}"
            for mutation_name, mutation_weight in config["builds"][wildcards.build_name]["mutation_weights"].items()
        ]),
        min_ha1_coverage=config["min_ha1_coverage"],
    shell:
        r"""
        python scripts/summarize_haplotypes_for_library_design.py \
            --metadata {input.metadata} \
            --aa-haplotype-fitnesses {input.aa_haplotype_fitnesses} \
            --emerging-haplotype-fitnesses {input.emerging_haplotype_fitnesses} \
            --mutation-weights {params.mutation_weights} \
            --min-ha1-coverage {params.min_ha1_coverage} \
            --output {output.library_design}
        """

rule design_libraries:
    input:
        expand("builds/{build_name}/library_design.tsv", build_name=list(config["builds"].keys()))
