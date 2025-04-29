# Nextstrain automation

> [!NOTE]
> External users can ignore this directory!
> This build config/customization is tailored for the internal Nextstrain team
> to extend the core ingest workflow for automated workflows.

## Run the workflow

### Via GitHub Actions

Go to the [ingest workflow](https://github.com/nextstrain/seasonal-flu/actions/workflows/ingest.yaml)
and trigger manually via the "Run workflow" dropdown.

You can also trigger the workflow via the [GitHub CLI](https://cli.github.com/):
```
gh workflow run ingest.yaml --repo nextstrain/seasonal-flu
```

### Via Nextstrain CLI

Provide the additional config file to the Snakemake options in order to
include the custom rules from [upload.smk](upload.smk) in the workflow.
Specify the `upload_all` target in order to run the additional upload rules.

The upload rules will require AWS credentials for a user that has permissions
to upload to the Nextstrain data bucket.

The customized workflow can be run from the top level pathogen repo directory with:
```
nextstrain build \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    ingest \
        upload_all \
        --configfile build-configs/nextstrain-automation/config.yaml
```
