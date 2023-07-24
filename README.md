# nextstrain.org/flu

[![Build Status](https://github.com/nextstrain/seasonal-flu/actions/workflows/ci.yaml/badge.svg?branch=master)](https://github.com/nextstrain/seasonal-flu/actions/workflows/ci.yaml)

This is the [Nextstrain](https://nextstrain.org) build for seasonal influenza viruses,
available online at [nextstrain.org/flu](https://nextstrain.org/flu).

The build encompasses fetching data, preparing it for analysis, doing quality control,
performing analyses, and saving the results in a format suitable for visualization (with
[auspice][]).  This involves running components of Nextstrain such as [fauna][] and
[augur][].

All influenza virus specific steps and functionality for the Nextstrain pipeline should be
housed in this repository.

This build is more complicated than other standard nextstrain build because all four
currently circulating seasonal influenza lineages (A/H3N2, A/H1N1pdm, B/Vic and B/Yam)
are analyzed using the same Snakefile with appropriate wildcards. In addition, we run
analyses of both the HA and NA segments of the influenza virus genome and analyze datasets
that span different time intervals (eg 2, 3, 6 years). Furthermore, the Nextstrain analysis
of influenza virus evolution also uses antigenic and serological data from different
WHO collaborating centers.

The different builds for the general public and the different WHO collaborating centers
are configured via separate config files. The Nextstrain build configs
([upload](profiles/upload.yaml), [nextstrain-public](profiles/nextstrain-public.yaml), [private.nextflu.org](profiles/private.nextflu.org.yaml))
are used for our semi-automated builds through our [GitHub Action workflows](.github/workflows/).

## Example build

You can run an example build using the example data provided in this repository.

First follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html)
for Nextstrain's suite of software tools.

Then run the example build via:

```
nextstrain build .  --configfile profiles/ci/builds.yaml
```

When the build has finished running, view the output Auspice trees via:

```
nextstrain view auspice/
```

## History

 - Prior to March 31, 2023, we selected strains for each build using a custom Python script called [select_strains.py](https://github.com/nextstrain/seasonal-flu/blob/64b5204d23c0b95e4b06f943e4efb8db005759c0/scripts/select_strains.py). With the merge of [the refactored workflow](https://github.com/nextstrain/seasonal-flu/pull/76), we have since used a configuration file to define the `augur filter` query logic we want for strain selection per build.

[Nextstrain]: https://nextstrain.org
[fauna]: https://github.com/nextstrain/fauna
[augur]: https://github.com/nextstrain/augur
[auspice]: https://github.com/nextstrain/auspice
[snakemake cli]: https://snakemake.readthedocs.io/en/stable/executable.html#all-options
[nextstrain-cli]: https://github.com/nextstrain/cli
[nextstrain-cli README]: https://github.com/nextstrain/cli/blob/master/README.md
