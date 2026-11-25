# Ingest open

Ingest open is the workflow for ingesting public seasonal-flu data from GenSpectrum.
This workflow only includes lineages B/Victoria, H1N1pdm, and H3N2.

## Workflow Usage

### With `nextstrain run`

If you haven't set up the seasonal-flu pathogen, then set it up with:

    nextstrain setup seasonal-flu

Otherwise, make sure you have the latest set up with:

    nextstrain update seasonal-flu

Run the ingest workflow with:

    nextstrain run seasonal-flu ingest-open <analysis-directory>

> [!TIP]
> If you would like to use a custom config to override defaults, add the
> config YAML to your analysis directory as `<analysis-directory>/config.yaml`.

Your `<analysis-directory>` will contain the workflow's intermediate files
and the final outputs:

- `results/{lineage}/metadata.tsv`
- `results/{lineage}/{segment}.fasta`


### With `nextstrain build`

If you prefer to run the workflow with full Snakemake CLI options, then you can
run the workflow with:

      nextstrain build ingest-open

> [!TIP]
> If you would like to use a custom config to override defaults, add the
> config YAML to `ingest-open/config.yaml` or point to a different location of
> the file that must be within the seasonal-flu repository with `--configfile`.

## Defaults

The defaults directory contains all of the default configurations for the ingest workflow.

[defaults/config.yaml](defaults/config.yaml) contains all of the default configuration parameters
used for the ingest workflow.

## Snakefile and rules

The rules directory contains separate Snakefiles (`*.smk`) as modules of the core ingest workflow.
The modules of the workflow are in separate files to keep the main ingest [Snakefile](Snakefile) succinct and organized.

Use the config parameter `custom_rules` to include additional rules for the workflow.
