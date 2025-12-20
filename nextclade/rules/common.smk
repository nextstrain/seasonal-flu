import datetime


wildcard_constraints:
    flu_type="[A-Za-z0-9]+",
    lineage=r"h3n2|h1n1pdm|vic|yam",
    segment = r"pb2|pb1|pa|ha|np|na|mp|ns",
    reference="[^/]+",


def genes(w):
    return {
        'ha': ["SigPep", "HA1", "HA2"],
        'na': ["NA"],
        'pb2': ["PB2"],
        'pb1': ["PB1"],
        'pa': ["PA"],
        'np': ["NP"],
        'mp': ["M1", "M2"],
        'ns': ["NEP", "NS1"]
    }.get(w.segment, [])


rule download_sequences:
    output:
        sequences="data/{lineage}/{segment}/sequences.fasta",
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/{segment}/sequences.fasta.xz",
    shell:
        """
        aws s3 cp {params.s3_path} - | xz -c -d > {output.sequences}
        """

rule download_metadata:
    output:
        metadata="data/{lineage}/{segment}/metadata-raw.tsv",
    params:
        s3_path="s3://nextstrain-data-private/files/workflows/seasonal-flu/{lineage}/metadata.tsv.xz",
    shell:
        """
        aws s3 cp {params.s3_path} - | xz -c -d > {output.metadata}
        """

rule get_nextclade_dataset:
    output:
        "nextclade/{lineage}_{segment}/reference.fasta",
    threads: 1
    params:
        nextclade_dataset_name=lambda w: f"flu_{w.lineage}_{w.segment}" if w.segment in ['ha', 'na'] else f"nextstrain/flu/{w.lineage}/{w.segment}"
    shell:
        """
        nextclade3 dataset get -n {params.nextclade_dataset_name} --output-dir nextclade/{wildcards.lineage}_{wildcards.segment}
        """

rule run_nextclade:
    input:
        sequences="data/{lineage}/{segment}/sequences.fasta",
        reference="nextclade/{lineage}_{segment}/reference.fasta",
    output:
        "data/{lineage}/{segment}/nextclade.tsv",
    threads: workflow.cores
    shell:
        """
        nextclade3 run -j {threads} -D nextclade/{wildcards.lineage}_{wildcards.segment} \
                  {input.sequences} --quiet --output-tsv {output}
        """


def get_clade_columns(w):
    return ",".join(["seqName", "qc.overallStatus"] + {
        'h3n2_ha':["clade", "subclade"],
        'h1n1pdm_ha':["clade", "subclade"],
        'vic_ha':["clade", "subclade"]
    }.get(f"{w.lineage}_{w.segment}", ["clade"]))


rule combined_with_metadata:
    input:
        nextclade="data/{lineage}/{segment}/nextclade.tsv",
        metadata="data/{lineage}/{segment}/metadata-raw.tsv"
    output:
        metadata="data/{lineage}/{segment}/metadata.tsv"
    params:
        nextclade_columns=get_clade_columns
    shell:
        """
        tsv-select -H -f {params.nextclade_columns} {input.nextclade} \
            | csvtk join -t --fields "strain;seqName" {input.metadata} /dev/stdin > {output.metadata}
        """


rule download_changelog_dataset:
    message:
        "Downloading previous dataset changelog for {wildcards.lineage} from {params.source} -> {output}"
    output:
        changelog = "data/{lineage}/{segment}/{reference}/dataset-changelog.md"
    params:
        source=lambda w: f"{config['dataset_repo']}/{w.lineage}/{w.segment}/{w.reference}/CHANGELOG.md",
    shell:
        """
        curl {params.source} > {output.changelog}
        """

rule sample_reference_strains:
    input:
        sequences="data/{lineage}/{segment}/sequences.fasta",
        include_strains="../config/{lineage}/reference_strains.txt",
    output:
        sampled_sequences="build/{lineage}/{segment}/reference_strains.fasta",
    shell:
        """
        seqkit grep -f {input.include_strains} {input.sequences} > {output.sampled_sequences}
        """


rule subsample:
    input:
        sequences="data/{lineage}/{segment}/sequences.fasta",
        metadata="data/{lineage}/{segment}/metadata.tsv",
        include_strains="../config/{lineage}/reference_strains.txt",
        nextclade_include="dataset_config/{lineage}/includes.txt",
        exclude="../config/{lineage}/outliers.txt",
    output:
        sampled_sequences="build/{lineage}/{segment}/{reference}/subsample_tmp.fasta",
        sampled_strains="build/{lineage}/{segment}/{reference}/subsample_tmp.txt",
    params:
        filter_arguments=lambda w: config["builds"][w.lineage][w.segment]["refs"][
            w.reference
        ]["filter"],
        reference_EPI_ISL=lambda w: config["builds"][w.lineage][w.segment]["refs"][
            w.reference
        ].get("reference_EPI_ISL", "none"),
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --include {input.include_strains} {input.nextclade_include} \
            --include-where gisaid_epi_isl={params.reference_EPI_ISL} \
            --exclude-where qc.overallStatus='bad' \
            {params.filter_arguments} \
            --output-sequences {output.sampled_sequences} \
            --output-strains {output.sampled_strains}
        """

rule subsample_harddate:
    input:
        sequences=rules.subsample.output.sampled_sequences,
        metadata="data/{lineage}/{segment}/metadata.tsv"
    output:
        sampled_sequences="build/{lineage}/{segment}/{reference}/subsample.fasta",
        sampled_strains="build/{lineage}/{segment}/{reference}/subsample.txt",
    params:
        hardmin=lambda w: config["builds"][w.lineage][w.segment]["refs"][
            w.reference
        ].get("hardmin_date", "1900"),
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-date {params.hardmin} \
            --output-sequences {output.sampled_sequences} \
            --output-strains {output.sampled_strains}
        """

rule align:
    input:
        sequences="build/{lineage}/{segment}/{reference}/subsample.fasta",
        reference_strains="build/{lineage}/{segment}/reference_strains.fasta",
        annotation="dataset_config/{lineage}/{segment}/{reference}/annotation.gff",
        reference="dataset_config/{lineage}/{segment}/{reference}/reference.fasta",
    output:
        alignment="build/{lineage}/{segment}/{reference}/align.aligned.fasta"
    params:
        outdir=lambda w: f"build/{w.lineage}/{w.segment}/{w.reference}/aligned.gene.{{cds}}.fasta",
        nextclade_bin = "nextclade3"
    threads: 3
    shell:
        """
        seqkit rmdup {input.sequences} {input.reference_strains} | \
        {params.nextclade_bin} run \
            --jobs={threads} \
            --input-ref {input.reference} \
            --input-annotation {input.annotation} \
            --gap-alignment-side right \
            --output-translations {params.outdir} \
            --output-fasta {output.alignment} \
            - \
            2>&1
        """


rule virus_specific_jsons:
    input:
        auspice_config= "config/auspice_config.json",
        pathogen = "config/pathogen.json",
        additional_pathogen="dataset_config/{lineage}/{segment}/{reference}/pathogen.json",
    output:
        pathogen = "build/{lineage}/{segment}/{reference}/pathogen.json",
        auspice = "build/{lineage}/{segment}/{reference}/auspice_config.json",
    params:
        reference_name = lambda w: config["builds"][w.lineage][w.segment]['refs'][w.reference]['reference_strain']
    shell:
        """
        python3 scripts/merge_jsons.py --lineage {wildcards.lineage} \
            --reference {wildcards.reference} \
            --reference-name {params.reference_name} \
            --segment {wildcards.segment} \
            --pathogen-jsons {input.pathogen} {input.additional_pathogen} \
            --auspice-config {input.auspice_config} \
            --output-pathogen {output.pathogen} \
            --output-auspice {output.auspice}
        """


rule generate_sample_sequences:
    input:
        sequences="data/{lineage}/{segment}/sequences.fasta",
        metadata="data/{lineage}/{segment}/metadata.tsv"
    output:
        sequences="build/{lineage}/{segment}/{reference}/sample_sequences.fasta",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-date 2020 --group-by year --subsample-max-sequences 50  \
            --exclude-ambiguous-dates-by year \
            --exclude-where 'country!=Usa' 'submitting_lab!=Centers For Disease Control And Prevention' \
            --probabilistic-sampling \
            --output-sequences {output.sequences}
        """

