# Repository for using DMS data from Welsh et al. to estimate phenotypes of circulating H3N2 influenza strains

This repository has all the code used to estimate phenotypes of circulating H3N2 influenza strains from Welsh et al. in the form of a [Nextstrain build](https://docs.nextstrain.org/en/latest/reference/glossary.html#term-build) (forked from the [core seasonal-flu workflow](https://github.com/nextstrain/seasonal-flu)). 

## Summary of strains used in the analysis
We analyzed a set of 1,478 H3N2 influenza strains obtained from a publicly available Nextstrain phylogenetic tree that samples strains from the past 12 years (https://nextstrain.org/flu/seasonal/h3n2/ha/12y).
This tree is regularly updated.
We analyzed the strains from the tree that was available on October 19th, 2023.
We obtained the names of these strains by downloading the `ACKNOWLEDGEMENTS (TSV)` file accessible via the `DOWNLOAD DATA` button.
The contents of this file are in `profiles/dmsa-phenotype/sequences/nextstrain_flu_seasonal_h3n2_ha_12y_acknowledgements_and_accession_numbers.tsv`, along with an added column (`accession_ha`) that gives the GISAID accession number of each HA sequence.
For each strain, we obtained the associated HA gene nucleotide sequence and metedata from the GISAID database (https://gisaid.org/), though we did not include these sequences and data in the repository as they cannot be publicly shared.

## Code and input data for estimating escape scores in a Nextstrain build
<!-- We forked this repository from https://github.com/nextstrain/seasonal-flu, which is a Nextstrain build for seasonal-influenza viruses, where all steps are combined into a single `snakemake` pipeline.
We added a step to this pipeline that uses DMS data to estimate escape scores for the input viruses.
The entire pipeline can be run with the following command, assuming all input files are present (the input sequences and metadata cannot be publicly shared, but can be obtained from GISAID; see above): -->



## Running the workflow

For those unfamiliar with Nextflow, we reccomend checking out the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/) before reading further.

In short, this multi-build workflow has been extended with additional functionality for estimating the phenotype (e.g. antibody escape scores) of the viral sequences within a tree when provided with Deep Mutational Scanning (DMS) data. The snakemake configuration for the workflow is unchanged, except you may specify an arbitrary number of DMS experiments under the `dmsa-phenotype-collections:` key. Each collection is defined by a nested list of key-value pairs defining the name and relavent parameters for a set of related experiments. Once these have been defined, the individual builds workflows (i.e. values under the `build` keys) in the config must then specify the names of the experimental collections you wish to be included in the [Auspice JSON's](https://docs.nextstrain.org/en/latest/reference/glossary.html#term-JSONs) resulting from the completed workflow.

<!-- describe the environment -->

<!-- describe the config -->

```Command-line argument used to run build```

<!-- Below, we describe the code and input data in more detail.
The standard `seasonal-flu` Nextstrain workflow includes the basic steps of reading in HA gene nucleotide sequences and metadata, aligning the sequences (both at the nucleotide and amino-acid level), and generating a phylogenetic tree annotated with strain-specific metadata.

We added `snakemake` rules that perform the following steps:
* one rule uses DMS data from a single serum sample to estimate an escape score for each of the input strains. Our strategy for estimating escape scores is described in more detail in the paper, and takes as input the amino-acid level multiple-sequence alignment from above. This code outputs CSV and JSON files with estimated escape scores.
* the JSON files are fed into the pipeline when generating the auspice output file to allow scores to be visualized on the tree.
* do we add any other rules? -->

<!-- All code and input data that we added are in the directory `profiles/dmsa-phenotype/`, including the following directories and files: -->



* `dmsa-pred/`: a sub-repository with Python scripts used to estimate phenotypes from DMS data.
    * `dmsa_rules.smk`
* `configs/`:
    * `.yaml`: defines all input parameters, including the name of the reference strain used in the DMS experiment, which is set to be `A/Hong Kong/45/2019`.
* `sequences/`:
    * `12y_sequences.fasta.xz`: a FASTA file with HA gene nucleotide sequences for each strain in the analysis. We cannot make this file public because the sequences are from GISAID, but we do provide GISAID accession numbers that can be used to obtain these sequences (see above).
    * `12y_metadata.tsv.xz`: a TSV file with metadata for each strain in the analysis. As above, we cannot make these data public, but they can be obtained using the above GISAID accession numbers.
* `antibody_escape/filtered_data/`: for each serum sample, this directory contains a CSV file reporting the DMS-measured effects of mutations on escape from that serum. These files have been filtered to only include entries for mutations that were seen in both of the replicate DMS libraries and in at least three distinct variants (averaged across libraries) in the DMS data for that serum sample.
* `allow_unmeasured_aa_subs_at_these_sites.txt`: a file listing the set of sites that are used when determining which strains to mask in the analysis. If a strain has one or more mutations that were not measured by DMS or that were not observed in at least two DMS libraries in an average of at least three variants per library, the strain is masked, unless the mutation in question occurs at a site listed in this file. We deemed that these sites are unlikely to affect escape by virtue of not being close to other sites that were found to affect escape.

## Files with estimated escape scores from the above pipeline
* location of estimated escape scores

## Analysis of estimated escape scores
* `analysis_code.ipynb`: a Jupyter notebook that analyzes the escape scores computed above. This notebook creates Figure 6B from the paper and performs the statistical analysis to test whether the relationship between escape score and sampling date is different for different cohorts.

To do:
* add link to paper
* can we make a public version of the file with estimated escape scores that lacks all mutation info?