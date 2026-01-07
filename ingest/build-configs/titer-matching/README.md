


1. Get the titers TSVs (gzipped) from S3 and store them in `data/titers/<subtype>/`
```
aws s3 cp --exclude "*" --include "*.tsv.gz" --recursive --dryrun s3://nextstrain-data-private/files/workflows/seasonal-flu/ .
```

2. Get the fauna metadata (will get some extra that we don't need but that's ok)
```
mkdir data/fauna-metadata
cd data/fauna-metadata/

aws s3 cp --exclude "*" --include "metadata.tsv.xz" --recursive --dryrun s3://nextstrain-data-private/files/workflows/seasonal-flu/ .

rm -rf trials
```

3. Curate (ingest) data as per our new snakemake pipeline
this gets files like `results/h3n2/metadata.tsv`

4. Run script
```
./scripts/titer-strain-matching.py --titers data/titers/ --fauna data/fauna-metadata/ --curated results/ --output results/titer-matching.png && open results/titer-matching.png
```

