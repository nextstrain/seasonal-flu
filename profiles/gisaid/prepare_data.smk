rule prepare_metadata:
    input:
        metadata="data/{lineage}/metadata.xls",
    output:
        metadata="data/{lineage}/metadata.tsv",
    params:
        old_fields=",".join(config["metadata_fields"]),
        new_fields=",".join(config["renamed_metadata_fields"]),
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/xls2csv.py --xls {input.metadata} --output /dev/stdout \
            | csvtk cut -f {params.old_fields} \
            | csvtk rename -f {params.old_fields} -n {params.new_fields} \
            | csvtk sep -f location --na "N/A" --names region,country,divison,location --merge --num-cols 4 --sep " / " \
            | csvtk replace -f strain -p " " -r "" \
            | csvtk sort -k strain,accession:r \
            | csvtk uniq -T -f strain > {output.metadata}
        """

rule prepare_sequences:
    input:
        sequences="data/{lineage}/raw_sequences_{segment}.fasta",
    output:
        sequences="data/{lineage}/{segment}.fasta",
    conda: "../../workflow/envs/nextstrain.yaml"
    shell:
        """
        seqkit replace -p " " -r "" {input.sequences} \
            | seqkit replace -p "\|" -r " " \
            | seqkit rename > {output.sequences}
        """
