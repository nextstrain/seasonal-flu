# shared

Shared internal tooling for pathogen workflows.  Used by our individual
pathogen repos which produce Nextstrain builds.  Expected to be vendored by
each pathogen repo using `git subrepo`.

Some tools may only live here temporarily before finding a permanent home in
`augur curate` or Nextstrain CLI.  Others may happily live out their days here.

## Vendoring

Nextstrain maintained pathogen repos will use [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to vendor shared scripts.
(See discussion on this decision in https://github.com/nextstrain/shared/issues/3)

For a list of Nextstrain repos that are currently using this method, use [this
GitHub code search](https://github.com/search?type=code&q=org%3Anextstrain+subrepo+%22remote+%3D+https%3A%2F%2Fgithub.com%2Fnextstrain%2Fingest%22).

If you don't already have `git subrepo` installed, follow the [git subrepo installation instructions](https://github.com/ingydotnet/git-subrepo#installation).
Then add the latest shared scripts to the pathogen repo by running:

```
git subrepo clone https://github.com/nextstrain/shared shared/vendored
```

Any future updates of sahred scripts can be pulled in with:

```
git subrepo pull shared/vendored
```

If you run into merge conflicts and would like to pull in a fresh copy of the
latest shared scripts, pull with the `--force` flag:

```
git subrepo pull shared/vendored --force
```

> **Warning**
> Beware of rebasing/dropping the parent commit of a `git subrepo` update

`git subrepo` relies on metadata in the `shared/vendored/.gitrepo` file,
which includes the hash for the parent commit in the pathogen repos.
If this hash no longer exists in the commit history, there will be errors when
running future `git subrepo pull` commands.

If you run into an error similar to the following:
```
$ git subrepo pull shared/vendored
git-subrepo: Command failed: 'git branch subrepo/shared/vendored '.
fatal: not a valid object name: ''
```
Check the parent commit hash in the `shared/vendored/.gitrepo` file and make
sure the commit exists in the commit history. Update to the appropriate parent
commit hash if needed.

## History

Much of this tooling originated in
[ncov-ingest](https://github.com/nextstrain/ncov-ingest) and was passaged thru
[mpox's ingest/](https://github.com/nextstrain/mpox/tree/@/ingest/). It
subsequently proliferated from [mpox][] to other pathogen repos ([rsv][],
[zika][], [dengue][], [hepatitisB][], [forecasts-ncov][]) primarily thru
copying.  To [counter that
proliferation](https://bedfordlab.slack.com/archives/C7SDVPBLZ/p1688577879947079),
this repo was made.

[mpox]: https://github.com/nextstrain/mpox
[rsv]: https://github.com/nextstrain/rsv
[zika]: https://github.com/nextstrain/zika/pull/24
[dengue]: https://github.com/nextstrain/dengue/pull/10
[hepatitisB]: https://github.com/nextstrain/hepatitisB
[forecasts-ncov]: https://github.com/nextstrain/forecasts-ncov

## Elsewhere

The creation of this repo, in both the abstract and concrete, and the general
approach to "ingest" has been discussed in various internal places, including:

- https://github.com/nextstrain/private/issues/59
- @joverlee521's [workflows document](https://docs.google.com/document/d/1rLWPvEuj0Ayc8MR0O1lfRJZfj9av53xU38f20g8nU_E/edit#heading=h.4g0d3mjvb89i)
- [5 July 2023 Slack thread](https://bedfordlab.slack.com/archives/C7SDVPBLZ/p1688577879947079)
- [6 July 2023 team meeting](https://docs.google.com/document/d/1FPfx-ON5RdqL2wyvODhkrCcjgOVX3nlXgBwCPhIEsco/edit)
- _…many others_

## Scripts

Scripts for supporting workflow automation that don’t really belong in any of our existing tools.

- [assign-colors](scripts/assign-colors) - Generate colors.tsv for augur export based on ordering, color schemes, and what exists in the metadata. Used in the phylogenetic or nextclade workflows.
- [notify-on-diff](scripts/notify-on-diff) - Send Slack message with diff of a local file and an S3 object
- [notify-on-job-fail](scripts/notify-on-job-fail) - Send Slack message with details about failed workflow job on GitHub Actions and/or AWS Batch
- [notify-on-job-start](scripts/notify-on-job-start) - Send Slack message with details about workflow job on GitHub Actions and/or AWS Batch
- [notify-on-record-change](scripts/notify-on-recod-change) - Send Slack message with details about line count changes for a file compared to an S3 object's metadata `recordcount`.
  If the S3 object's metadata does not have `recordcount`, then will attempt to download S3 object to count lines locally, which only supports `xz` compressed S3 objects.
- [notify-slack](scripts/notify-slack) - Send message or file to Slack
- [s3-object-exists](scripts/s3-object-exists) - Used to prevent 404 errors during S3 file comparisons in the notify-* scripts
- [trigger](scripts/trigger) - Triggers downstream GitHub Actions via the GitHub API using repository_dispatch events.
- [trigger-on-new-data](scripts/trigger-on-new-data) - Triggers downstream GitHub Actions if the provided `upload-to-s3` outputs do not contain the `identical_file_message`
  A hacky way to ensure that we only trigger downstream phylogenetic builds if the S3 objects have been updated.


NCBI interaction scripts that are useful for fetching public metadata and sequences.

- [fetch-from-ncbi-entrez](scripts/fetch-from-ncbi-entrez) - Fetch metadata and nucleotide sequences from [NCBI Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25501/) and output to a GenBank file.
  Useful for pathogens with metadata and annotations in custom fields that are not part of the standard [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) outputs.

Historically, some pathogen repos used the undocumented NCBI Virus API through [fetch-from-ncbi-virus](https://github.com/nextstrain/shared/blob/c97df238518171c2b1574bec0349a55855d1e7a7/fetch-from-ncbi-virus) to fetch data. However we've opted to drop the NCBI Virus scripts due to https://github.com/nextstrain/shared/issues/18.

Potential Nextstrain CLI scripts

- [sha256sum](scripts/sha256sum) - Used to check if files are identical in upload-to-s3 and download-from-s3 scripts.
- [cloudfront-invalidate](scripts/cloudfront-invalidate) - CloudFront invalidation is already supported in the [nextstrain remote command for S3 files](https://github.com/nextstrain/cli/blob/a5dda9c0579ece7acbd8e2c32a4bbe95df7c0bce/nextstrain/cli/remote/s3.py#L104).
  This exists as a separate script to support CloudFront invalidation when using the upload-to-s3 script.
- [upload-to-s3](scripts/upload-to-s3) - Upload file to AWS S3 bucket with compression based on file extension in S3 URL.
  Skips upload if the local file's hash is identical to the S3 object's metadata `sha256sum`.
  Adds the following user defined metadata to uploaded S3 object:
    - `sha256sum` - hash of the file generated by [sha256sum](sha256sum)
    - `recordcount` - the line count of the file
- [download-from-s3](scripts/download-from-s3) - Download file from AWS S3 bucket with decompression based on file extension in S3 URL.
  Skips download if the local file already exists and has a hash identical to the S3 object's metadata `sha256sum`.

## Snakemake

Snakemake workflow functions that are shared across many pathogen workflows that don’t really belong in any of our existing tools.

- [config.smk](snakemake/config.smk) - Shared functions for parsing workflow configs.
- [remote_files.smk](snakemake/remote_files.smk) - Exposes the `path_or_url` function which will use Snakemake's storage plugins to download/upload files to remote providers as needed.


## Software requirements

Some scripts may require Bash ≥4. If you are running these scripts on macOS, the builtin Bash (`/bin/bash`) does not meet this requirement. You can install [Homebrew's Bash](https://formulae.brew.sh/formula/bash) which is more up to date.

## Testing

Most scripts are untested within this repo, relying on "testing in production". That is the only practical testing option for some scripts such as the ones interacting with S3 and Slack.

## Working on this repo

This repo is configured to use [pre-commit](https://pre-commit.com),
to help automatically catch common coding errors and syntax issues
with changes before they are committed to the repo.

If you will be writing new code or otherwise working within this repo,
please do the following to get started:

1. [install `pre-commit`](https://pre-commit.com/#install) by running
   either `python -m pip install pre-commit` or `brew install
   pre-commit`, depending on your preferred package management
   solution
2. install the local git hooks by running `pre-commit install` from
   the root of the repo
3. when problems are detected, correct them in your local working tree
   before committing them.

Note that these pre-commit checks are also run in a GitHub Action when
changes are pushed to GitHub, so correcting issues locally will
prevent extra cycles of correction.
