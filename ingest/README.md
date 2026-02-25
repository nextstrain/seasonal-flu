# Ingest

The base ingest workflow ingests manually downloaded data from GISAID and outputs curated
metadata and sequences that can be used as input for the phylogenetic workflow.

If you have another data source or private data that needs to be formatted for
the phylogenetic workflow, then you can use a similar workflow to curate your
own data.

> There are additional ingest workflows which Nextstrain uses to manage seasonal-flu
> builds. See the READMEs in `build-configs/{manual-upload,nextstrain-automation}`
> for more details.

## Workflow Usage
Manually download the metadata and sequences from GISAID.
Save the XLS metadata file as `ingest/data/YYYY-MM-DD-N-metadata.xls`.
Make sure the sequences includes all segments and the FASTA header is set to
```
DNA Accession no. | Submitting lab  | Originating lab
```
Save the FASTA file as `ingest/data/YYYY-MM-DD-N-sequences.fasta`.

The workflow can be run from the top level pathogen repo directory:
```
nextstrain build ingest
```

This produces the default outputs of the ingest workflow:

- metadata      = results/<lineage>/metadata.tsv
- sequences     = results/<lineage>/<segment>.fasta


If your downloaded data only includes a subset of lineages, e.g. only include h3n2,
then you can disable lineages by providing a custom config

```yaml
filtering:
  h1n1pdm: ~
  vic: ~
  yam: ~
  b: ~
  h1n1: ~
  h2n2: ~
  avian-flu: ~
```

```
nextstrain build ingest --configfile custom-config.yaml
```


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
