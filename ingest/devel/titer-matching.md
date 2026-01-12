# Cross reference titer strain names against metadata strain names


## Files needed:

Titers TSVs (gzipped) from S3 and store them in `data/titers/<subtype>/`

```sh
mkdir -p data/titers
cd data/titers
aws s3 cp --exclude "*" --include "*.tsv.gz" --recursive --dryrun s3://nextstrain-data-private/files/workflows/seasonal-flu/ .
```


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

## Run:

```sh
./devel/titer-matching.py \
  --titers data/titers \
  --fauna data \
  --curated results \
  --output results/titer-matching.png \
  --log results/titer-matching.tsv \
  && open results/titer-matching.png
```
