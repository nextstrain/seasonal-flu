# Monthly reports on seasonal influenza evolution

## Build trees

Clone this repository and change into its directory.

``` bash
git clone https://github.com/nextstrain/seasonal-flu.git seasonal-flu-aws
cd seasonal-flu-aws

# Check out this development branch.
git checkout refactor-workflow
```

Setup your environment variables for fauna.

``` bash
. ~/environment_rethink.sh
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
Download the group overview markdown file.
Your Nextstrain Group authorization will grant you permission to access the group overview file through the REST API with curl.

``` bash
curl https://nextstrain.org/groups/nextflu-private/settings/overview \
    --header "`nextstrain authorization`" \
    --header 'Content-Type: text/markdown' \
    -o group-overview.md
```

Edit this file to include links to the latest report and builds by date.
Upload the updated overview to the group.

``` bash
curl https://nextstrain.org/groups/nextflu-private/settings/overview \
    --header "`nextstrain authorization`" \
    --header 'Content-Type: text/markdown' \
    --upload-file group-overview.md
```

In the future, we will replace these curl commands with `nextstrain remote` commands that use the group customization endpoints through nextstrain.org.
