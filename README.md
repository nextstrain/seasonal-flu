# nextstrain.org/flu

This is the [Nextstrain](https://nextstrain.org) build for seasonal influenza viruses, available online at [nextstrain.org/flu](https://nextstrain.org/flu).

The build encompasses fetching data, preparing it for analysis, doing quality
control, performing analyses, and saving the results in a format suitable for
visualization (with [auspice][]).  This involves running components of
Nextstrain such as [fauna][] and [augur][].

All influenza virus specific steps and functionality for the Nextstrain pipeline should be
housed in this repository.

This build is more complicated than other standard nextstrain build because all four currently circulating seasonal influenza lineages (A(H3N2), A(H1N1pdn), B(vic) and B(yam)) are analyzed using the same snakefile with appropriate wildcards.
In addition, we run analysis of the HA and NA segment of the influenza virus genome and analyze data sets that span different time intervals (e.g. 2, 6, 12 years).

Furthermore, the nextstrain analysis of influenza virus evolution also uses antigenic and serological data from different WHO collaborating centers.
These antigenic data come in four flavors depending on the assay that passage history of the antigens.
The influenza virus output files have the wildcard set

`{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}`

that currently use the following values:

 * center: [`who`, `cdc`, `vidrl`, `crick`, `niid`]
 * lineage: [`h3n2`, `h1n1pdm`, `vic`, `yam`]
 * segment: [`ha`, `na`]
 * resolution: [`12y`, `6y`, `3y`, `2y`]
 * assay: [`hi`, `cell`]
 * passage: [`cell`, `egg`]


To manage both builds for the general public and the different WHO collaborating centers, the Snakefiles are split into a `Snakefile_base` that contains the rules for the core analysis and the files

 * `Snakefile` for the standard build
 * `Snakefile_WHO` for the WHO CC builds
 * `Snakefile_reports` to generate figures and additional analysis for the biannual reports to the WHO

The latter snakefiles import the rules specified in `Snakefile_base`, define additional rules, and specify the build targets.


### fauna / RethinkDB credentials

This build starts by pulling sequences from our live [fauna][] database (a RethinkDB instance). This
requires environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY` to be set.


[Nextstrain]: https://nextstrain.org
[fauna]: https://github.com/nextstrain/fauna
[augur]: https://github.com/nextstrain/augur
[auspice]: https://github.com/nextstrain/auspice
[snakemake cli]: https://snakemake.readthedocs.io/en/stable/executable.html#all-options
[nextstrain-cli]: https://github.com/nextstrain/cli
[nextstrain-cli README]: https://github.com/nextstrain/cli/blob/master/README.md
