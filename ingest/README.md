# Ingest

This workflow ingests manually downloaded data from GISAID and outputs curated
metadata and sequences that can be used as input for the phylogenetic workflow.

If you have another data source or private data that needs to be formatted for
the phylogenetic workflow, then you can use a similar workflow to curate your
own data.

## Workflow Usage

The workflow can be run from the top level pathogen repo directory:
```
nextstrain build ingest
```

Alternatively, the workflow can also be run from within the ingest directory:
```
cd ingest
nextstrain build .
```

This produces the default outputs of the ingest workflow:

- metadata      = results/metadata.tsv
- sequences     = results/<segment>/sequences.fasta


## Defaults

The defaults directory contains all of the default configurations for the ingest workflow.

[defaults/config.yaml](defaults/config.yaml) contains all of the default configuration parameters
used for the ingest workflow. Use Snakemake's `--configfile`/`--config`
options to override these default values.

## Snakefile and rules

The rules directory contains separate Snakefiles (`*.smk`) as modules of the core ingest workflow.
The modules of the workflow are in separate files to keep the main ingest [Snakefile](Snakefile) succinct and organized.

The `workdir` is hardcoded to be the ingest directory so all filepaths for
inputs/outputs should be relative to the ingest directory.

Modules are all [included](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes)
in the main Snakefile in the order that they are expected to run.
