# Monthly reports on seasonal influenza evolution

## Build and deploy trees

[Run the nextflu-private builds with the GitHub Action](https://github.com/nextstrain/seasonal-flu/actions/workflows/run-nextflu-private-builds.yaml).
The workflow automatically deploys dated Auspice JSONs to [the nextflu-private group](https://nextstrain.org/groups/nextflu-private/) (e.g., https://nextstrain.org/groups/nextflu-private/flu/seasonal/2024-01-08/h3n2/2y/ha).
View the GitHub Action summary for details about how to download the build artifacts from AWS Batch.

## Prepare a report

Start with an existing narrative Markdown file as a template for the new report.

``` bash
mkdir -p narratives/
nextstrain remote download \
  https://nextstrain.org/groups/nextflu-private/narratives/nextstrain-cdc/2022-10-03 \
  narratives/
```

Edit the narrative file in your favorite editor, using URLs for the builds you just uploaded to the group as the content.
As you edit the file, upload it to the group and verify the content in the browser as you go.

``` bash
nextstrain remote upload \
  https://nextstrain.org/groups/nextflu-private/ \
  narratives/nextstrain-cdc_2022-10-03.md
```

The nextflu-private workflow automatically generates several outputs used for reporting including:

 - Plots of the number of genomes available in GISAID for all lineages (e.g., `figures/total-sample-count-by-lineage.png`)
 - Markdown tables of counts of samples collected in the last month per build and clade (e.g., `builds/h1n1pdm_2y_titers/counts_of_recent_tips_by_clade.md`)
 - Markdown tables of current frequency, change in frequency in the last 4 weeks, and titer reference strains for derived haplotypes (e.g., `builds/h1n1pdm_2y_titers/ha/haplotype_summary/cell_hi.md` per titer collection)

Copy these files (or their contents) into the corresponding sections of the narrative Markdown.

## Update the group overview

[The nextflu-private group's overview](https://nextstrain.org/groups/nextflu-private/) lists recent reports and trees by date with the most recent at the top of the page.
[Open the group's settings page](https://nextstrain.org/groups/nextflu-private/settings) and update the list of builds to include the builds you just uploaded.
