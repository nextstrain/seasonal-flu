This is the changelog for the Nextstrain seasonal-flu workflow.
All notable changes in a release will be documented in this file.

This changelog is intended for _humans_ and follows many of the principles from [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

Changes for this project _do not_ currently follow the [Semantic Versioning rules](https://semver.org/spec/v2.0.0.html).
Instead, changes appear below grouped by the date they were added to the workflow.
The "__NEXT__" heading below describes changes in the unreleased development source code and as such may not be routinely kept up to date.

# 31 July 2025

 - Fix parsing of GISAID metadata and sequences for the quickstart guide including removing newlines that appear within Excel cells of the metadata. See [#242](https://github.com/nextstrain/seasonal-flu/pull/242) for details.
 - Generate a Nextclade annotations TSV as part of the standard phylogenetic workflow and merge those annotations with the complete input metadata prior to subsampling. With this change, users can now define subsampling filters on Nextclade columns and refer to these columns in other parts of the workflow (e.g., in Auspice config JSONs as colorings). For example, to keep only sequences with a Nextclade overall quality status of "good", users can define a filter query like ```--query "\`qc.overallStatus\` == 'good'"``` (note the use of backticks to escape the QC column name with a `.` in it and the use of backward slashes to escape those backticks in the shell). We also now provide the Nextclade QC status as a coloring for HA trees. See [#240](https://github.com/nextstrain/seasonal-flu/pull/240) for details.

# 16 May 2025

 - Added the ability to annotate haplotype labels for viruses based on the specific combination of nucleotide or amino acid substitutions that appear in their sequences. [See the pull request for more details](https://github.com/nextstrain/seasonal-flu/pull/221).
