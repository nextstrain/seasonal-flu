# Nextstrain manual upload

> [!NOTE]
> External users can ignore this directory!
> This build is tailored for the internal Nextstrain team to
> run manually to upload GISAID files to our private AWS S3 bucket.


## Run the workflow

This workflow is expected to be run manually after downloading files from GISAID.
The GISAID files are expected to be saved as

- `<YYYY-MM-DD-N>-metadata.xls`
- `<YYYY-MM-DD-N>-sequences.fasta`

`<YYYY-MM-DD>` is the date the files were downloaded from GISAID.
`<N>` is the number of the download since GISAID limits the number of records per download.

For example, if you had to split the data between two downloads on 2025-04-11,
then save the files as
- `2025-04-11-1-metadata.xls`
- `2025-04-11-1-sequences.fasta`
- `2025-04-11-2-metadata.xls`
- `2025-04-11-2-sequences.fasta`

The directory in which you save the GISAID files depends on which command you are
using to run the workflow.

### With `nextstrain run`

When running with `nextstrain run`, you can save the GISAID files in any
arbitrary analysis directory. However, you must also create a `config.yaml` within
the analysis directory to specify the `gisaid_pairs` to upload.

Continuing the example above, your analysis directory should look like
```
<analysis-dir>
├── config.yaml
└── data
    ├── 2025-04-11-1-metadata.xls
    ├── 2025-04-11-1-sequences.fasta
    ├── 2025-04-11-2-metadata.xls
    └── 2025-04-11-2-sequences.fasta
```

With the `config.yaml` specifying the `gisaid_pairs` you want to upload

```yaml
gisaid_pairs:
 - 2025-04-11-1
 - 2025-04-11-2
```

Make sure you have the latest seasonal flu pathogen setup.

```
nextstrain update seasonal-flu@master
```

Then run the workflow
```
nextstrain run \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    seasonal-flu@master \
    ingest/build-configs/manual-upload \
    <analysis-directory>
```

### With `nextstrain build`

When running with `nextstrain build` the files must be saved _within_ the
seasonal-flu repo.

Save the downloaded GISAID metadata and sequences as:
- `ingest/build-configs/manual-upload/data/<YYYY-MM-DD-N>-metadata.xls`
- `ingest/build-configs/manual-upload/data/<YYYY-MM-DD-N>-sequences.fasta`

The workflow can be run from the top level pathogen repo directory with:
```
nextstrain build \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    ingest/build-configs/manual-upload \
        --config gisaid_pairs=["2025-04-11-1", "2025-04-11-2"]
```

### Required environment variables

You need to have AWS credentials with permissions to upload to the private
AWS S3 bucket `nextstrain-data-private`

- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
