# Strain name updates

This workflow / code is to facilitate our move from fauna-ingested data to our new curation pipeline.
Specifically, strain names have changed and we need to update them where they are being used as the canonical IDs.
Once we have sunsetted fauna this file & associated code should be removed.

### Files needed

The (seasonal-flu) outputs of the new ingest workflow, including intermediate (temp) files:
```sh
snakemake --cores 2 -pf --notemp --forceall \
  --configfile build-configs/nextstrain-automation/config.yaml --config gisaid_pairs='["gisaid_cache"]' \
  results/h3n2/metadata.tsv results/vic/metadata.tsv results/h1n1pdm/metadata.tsv results/yam/metadata.tsv
```

The previous fauna outputs:
```sh
for SUBTYPE in h3n2 h1n1pdm vic yam; do
  aws s3 cp s3://nextstrain-data-private/files/workflows/seasonal-flu/${SUBTYPE}/metadata.tsv.xz - | xz -c -d > data/${SUBTYPE}/fauna-metadata.tsv
done
```

### How to run

```sh
snakemake --cores 5 --snakefile devel/strain-name-updates.snakefile -pf 
```

> Note that the workflow produces files in `results` and doesn't overwrite the originals. You should do this manually.
