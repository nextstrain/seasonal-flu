import os
configfile: os.path.join(workflow.basedir, "strain-name-updates.yaml")
include: "../rules/remote_files.smk"

# print(f"{config=}")

rule all:
    input:
        files = [f"results/strain-name-updates/{d}/{f}" for d in config['strain-lists'] for f in config['strain-lists'][d]]


rule curated_strains:
    input:
        "results/{dataset}/metadata.tsv",
    output:
        "results/strain-name-updates/{dataset}/curated-strains.txt"
    run:
        from augur.io import read_metadata
        m = read_metadata(input[0])
        with open(output[0], 'w') as fh:
            print("\n".join(m.index.tolist()), file=fh)

rule fauna_strains:
    input:
        # "data/{dataset}/fauna-metadata.tsv",
        "../../avian-flu/ingest/fauna/data/metadata_combined.tsv",
    output:
        "results/strain-name-updates/{dataset}/fauna-strains.txt"
    run:
        from augur.io import read_metadata
        m = read_metadata(input[0])
        with open(output[0], 'w') as fh:
            print("\n".join(m.index.tolist()), file=fh)

# rule compute_strain_maps_via_epi_isl_lookup:
#     """
#     Computes maps of (old) fauna strain names to (new) curated strain names, where possible.
#     """
#     input:
#         fauna="data/{dataset}/fauna-metadata.tsv",
#         curated="data/curated_gisaid.ndjson.zst",
#     output:
#         tsv = "results/strain-name-updates/{dataset}/strain-name-map.tsv"
#     wildcard_constraints:
#         dataset= "h3n2|h1n1pdm|vic|yam"
#     shell:
#         """
#         ./devel/match-strain-names-via-epi-isl.py \
#             --fauna {input.fauna} \
#             --curated {input.curated} \
#             --output {output.tsv}
#         """
 

rule match_strains:
    input:
        fauna="results/strain-name-updates/{dataset}/fauna-strains.txt",
        curated="results/strain-name-updates/{dataset}/curated-strains.txt",
        strains=lambda w: path_or_url(config['strain-lists'][w.dataset][w.file_type]),
        # strain_map_tsv="results/strain-name-updates/{dataset}/strain-name-map.tsv",
    output:
        strains="results/strain-name-updates/{dataset}/{file_type}",
    shell:
        """
        ./devel/update-strains.py \
            --fauna-strains {input.fauna:q} \
            --curated-strains {input.curated:q} \
            --query-strains {input.strains:q} \
            > {output.strains:q}
        """