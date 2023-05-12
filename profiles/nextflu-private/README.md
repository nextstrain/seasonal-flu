# Monthly reports on seasonal influenza evolution

## Build trees

Clone this repository and change into its directory.

``` bash
git clone https://github.com/nextstrain/seasonal-flu.git seasonal-flu-aws
cd seasonal-flu-aws
```

Run the build in dry-run mode, to make sure you see all of the rules you expect to see.

``` bash
nextstrain build . -n -p --configfile profiles/nextflu-private.yaml
```

Run the build on AWS Batch.
This step assumes your AWS credentials live in the ambient environment or in the `~/.aws/credentials` file.
Adjust CPUs and memory as needed.
This command sets the max number of CPUs used for tree building to an odd number and total CPUs requested to an even number, so we always have at least one CPU free for other steps in the workflow.

``` bash
nextstrain build \
  --aws-batch \
  --cpus 32 \
  --memory 48Gib \
  . \
    -p \
    --configfile profiles/nextflu-private.yaml \
    --set-threads tree=7
```

Detach the job from the current terminal with Ctrl-Z and save the AWS Batch job id to reattach later.
After reattaching, wait for files to download, and then explore the trees locally with auspice.us or Auspice.

Alternately, run the jobs on the Hutch's SLURM cluster with the following command.

```bash
snakemake \
  -j 10 \
  --configfile profiles/nextflu-private.yaml \
  --profile profiles/hutch
```

If trees look reasonable, rename them to match the current date (or the date corresponding to when the data were updated).

``` bash
export DATE=2022-10-03
cd auspice/

for file in {h1n1pdm,h3n2,vic}*
do
  mv ${file} "flu_seasonal_${DATE}_${file}"
done
```

Upload the files to the Nextstrain group.

``` bash
nextstrain remote upload \
  https://nextstrain.org/groups/nextflu-private/ \
  *.json
```

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
  --use-conda \
  --conda-frontend mamba \
  -j 1 \
  -p \
  --configfile profiles/nextflu-private.yaml \
  figures/total-sample-count-by-lineage.png
```
