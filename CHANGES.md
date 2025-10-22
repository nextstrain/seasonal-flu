This is the changelog for the Nextstrain seasonal-flu workflow.
All notable changes in a release will be documented in this file.

This changelog is intended for _humans_ and follows many of the principles from [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

Changes for this project _do not_ currently follow the [Semantic Versioning rules](https://semver.org/spec/v2.0.0.html).
Instead, changes appear below grouped by the date they were added to the workflow.
The "__NEXT__" heading below describes changes in the unreleased development source code and as such may not be routinely kept up to date.

# 22 October 2025

- Set default clade and branch label to "subclade" for all public Nextstrain builds. See commit [6fb8478](https://github.com/nextstrain/seasonal-flu/commit/6fb8478234f65a576210ef6222899751d04a811c).
- Nextclade: Migrate unmaintained clade labels to a new `legacy-clade` column, report subclades in the `clade` column, and retain the original `subclade` column for backward compatibility. See [#262](https://github.com/nextstrain/seasonal-flu/pull/262/) for details.

# 30 September 2025

- Add new vaccine strains for H1N1pdm and H3N2 based on [the WHO recommendations for the Southern Hemisphere 2026](https://www.who.int/news/item/26-09-2025-recommendations-announced-for-influenza-vaccine-composition-for-the-2026-southern-hemisphere-influenza-season). See [#256](https://github.com/nextstrain/seasonal-flu/pull/256) for details.

# 17 September 2025

- Allow titer collections to define their own list of `genes` which get used for fitting titer models. For example, setting `genes: ["HA1"]` in a titer collection limits the titer substitution model to use only HA1 substitutions.

# 15 September 2025

- Add `force_run_nextclade` top-level configuration parameter to force the workflow to run Nextclade from scratch even when files already exist on S3.
- Add `s3_path` top-level configuration parameter to override the default path for builds.

# 18 August 2025

- Explicitly root divergence tree with reference node before pruning reference. See [#246](https://github.com/nextstrain/seasonal-flu/pull/246) for details.

# 11 August 2025

- Parameterize the local clock filter cutoff as a top-level build configuration option. This change allows users to effectively disable the threshold-based filtering on the z-score from the local clock filter by setting the cutoff to a high value or make the threshold more stringent. See [#245](https://github.com/nextstrain/seasonal-flu/pull/245) for details.

# 8 August 2025

- Parameterize frequencies options by adding a new build-level config parameter section `frequencies` which allows users to override parameters including `narrow_bandwidth`, `wide_bandwidth`, `proportion_wide`, `pivot_interval`, and `pivot_interval_units`. These parameters map to the corresponding command line arguments for `augur frequencies`. See [#243](https://github.com/nextstrain/seasonal-flu/pull/243) for details.

# 1 August 2025

- Add optional configuration parameter, `nextclade_server`, to specify a Nextclade dataset server to download Nextclade datasets. This parameter allows users to run Nextclade from GitHub branches where new clades are being defined or other URLs with custom datasets. For example, to run the GISAID quickstart workflow with the July 2025 proposed subclades, run the following command: `nextstrain build . --configfile profiles/gisaid/builds.yaml -np --config nextclade_server="https://raw.githubusercontent.com/nextstrain/nextclade_data/refs/heads/flu-update-2025-07/data_output"`.
- Update emerging haplotype definitions to match new clade definitions (e.g., replaces H3N2 HA haplotype J.2:158K-189R with the new clade J.2.3). See [#241](https://github.com/nextstrain/seasonal-flu/pull/241) for details.

# 31 July 2025

 - Fix parsing of GISAID metadata and sequences for the quickstart guide including removing newlines that appear within Excel cells of the metadata. See [#242](https://github.com/nextstrain/seasonal-flu/pull/242) for details.
 - Generate a Nextclade annotations TSV as part of the standard phylogenetic workflow and merge those annotations with the complete input metadata prior to subsampling. With this change, users can now define subsampling filters on Nextclade columns and refer to these columns in other parts of the workflow (e.g., in Auspice config JSONs as colorings). For example, to keep only sequences with a Nextclade overall quality status of "good", users can define a filter query like ```--query "\`qc.overallStatus\` == 'good'"``` (note the use of backticks to escape the QC column name with a `.` in it and the use of backward slashes to escape those backticks in the shell). We also now provide the Nextclade QC status as a coloring for HA trees. See [#240](https://github.com/nextstrain/seasonal-flu/pull/240) for details.

# 16 May 2025

 - Added the ability to annotate haplotype labels for viruses based on the specific combination of nucleotide or amino acid substitutions that appear in their sequences. [See the pull request for more details](https://github.com/nextstrain/seasonal-flu/pull/221).
