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

## Update the group overview

[The nextflu-private group's overview](https://nextstrain.org/groups/nextflu-private/) lists recent reports and trees by date with the most recent at the top of the page.
[Open the group's settings page](https://nextstrain.org/groups/nextflu-private/settings) and update the list of builds to include the builds you just uploaded.

## Plot counts per lineage for reports

Plot the number of genomes available in GISAID for all lineages in a single plot and then plot number per region within each lineage.
Plots use Altair which is not part of the Nextstrain Docker or Conda base environments, so we use a custom Conda environment with Snakemake.

``` bash
snakemake \
  figures/total-sample-count-by-lineage.png \
  --use-conda \
  --conda-frontend mamba \
  -j 1 \
  -p \
  --configfile profiles/nextflu-private.yaml
```
