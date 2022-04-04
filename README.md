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
that span different time intervals (eg 2, 3, 6 years).

Furthermore, the Nextstrain analysis of influenza virus evolution also uses antigenic and
serological data from different WHO collaborating centers. These antigenic data come in
four flavors depending on the assay that passage history of the antigens. The influenza
virus output files have the wildcard set

`{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}`

that currently use the following values:

* center: [`who`, `cdc`, `crick`, `niid`, `vidrl`]
* lineage: [`h3n2`, `h1n1pdm`, `vic`, `yam`]
* segment: [`ha`, `na`]
* resolution: [`6m`, `2y`, `3y`, `6y`, `12y`]
* assay: [`hi`, `fra`]
* passage: [`cell`, `egg`]

Intermediate files follow this wildcard ordering, but may omit irrelevant wildcards, eg
`filtered_h3n2_ha.fasta`.

To manage both builds for the general public and the different WHO collaborating centers,
the Snakefiles are split into a `Snakefile_base` that contains the rules for the core
analysis and the files, alongside:

* `Snakefile` for the standard "live" build housed at
  [nextstrain.org/flu](https://nextstrain.org/flu)
* `Snakefile_WHO` for the WHO CC builds
* `Snakefile_report` to generate figures and additional analysis for the biannual reports
to the WHO

The latter Snakefiles import the rules specified in `Snakefile_base`, define additional
rules, and specify the build targets.

### fauna / RethinkDB credentials

This build starts by pulling sequences from our live [fauna][] database (a RethinkDB
instance). This requires environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY` to be
set.

If you don't have access to our database, you can run the build using the example data
provided in this repository. Before running the build, copy the example sequences into the
`data/` directory like so:

```
mkdir data/
cp example_data/* data/
```

Then run the the build via:

```
nextstrain build . targets/flu_seasonal_h3n2_ha_12y
```

[Nextstrain]: https://nextstrain.org
[fauna]: https://github.com/nextstrain/fauna
[augur]: https://github.com/nextstrain/augur
[auspice]: https://github.com/nextstrain/auspice
[snakemake cli]: https://snakemake.readthedocs.io/en/stable/executable.html#all-options
[nextstrain-cli]: https://github.com/nextstrain/cli
[nextstrain-cli README]: https://github.com/nextstrain/cli/blob/master/README.md
