# Strain name matching

With our move away from fauna and to a curated all-influenza ingest pipeline we have introduced a number of strain name changes. This immediately presents a problem as we have myriad lists of hardcoded strain names, such as outliers-to-drop and force-include-lists. This workflow (which is rather ad-hoc!) attempts to match up the old (i.e. fauna) strain names with their updated strain names.

We use a combination of fuzzy-matching and, where possible a hardcoded map of fauna strain names to new strain names.

## How to run

```
cd ingest
snakemake --cores 4 --snakefile build-configs/strain-name-matching/Snakefile -pf -n
```

## What files you'll need

* The original fauna metadata files e.g. (from the base directory) `aws s3 cp s3://nextstrain-data-private/files/workflows/seasonal-flu/h1n1pdm/metadata.tsv.xz - | xz -c -d > data/h1n1pdm/metadata.tsv`. Note that this file won't represent fauna data forever!

* A number of files from the normal ingest pipeline. `snakemake --cores 2 --config gisaid_pairs='["gisaid_cache"]' -pf data/curated_gisaid.ndjson data/avian-flu/curated_gisaid.ndjson results/h3n2/metadata.tsv results/h1n1pdm/metadata.tsv results/vic/metadata.tsv results/yam/metadata.tsv`


## Hardcoded strain maps

For seasonal-flu datasets we can query the EPI_ISLs of the fauna data against the curated data to create lookups. 
This table reports how many of the fauna strain names have matches in our data. For avian-flu it's a little tricker so we leverage the existing diff-avian-flu script to create lookups.


| Dataset  | Updated | Missing | Unchanged |
| -------- | ------- | -- | -- |
| H1N1pdm  | 1,729    | 1,659 | 149,143
| H3N2  | 2,125    | 2,674 | 177,009
| vic  | 1,338    | 286 | 66,283
| yam  | 367    | 10 | 21,946
| avian-flu | 34,426 | 1,531 | 26,754



## TITERS _WIP!



```
cat ../builds/seasonal-flu_h3n2_3y/all_titers.tsv | cut -f 1 | sort | uniq > results/strain-name-matching/h3n2/titers_virus_strain.txt
cat ../builds/seasonal-flu_h3n2_3y/all_titers.tsv | cut -f 2 | sort | uniq > results/strain-name-matching/h3n2/titers_serum_strain.txt

./scripts/strain-name-fuzzer.py --curated-strains results/strain-name-matching/h3n2/metadata.txt --query-strains results/strain-name-matching/h3n2/titers_serum_strain.txt --strain-map results/strain-name-matching/h3n2/strain-name-map.tsv > /dev/null

./scripts/strain-name-fuzzer.py --curated-strains results/strain-name-matching/h3n2/metadata.txt --query-strains results/strain-name-matching/h3n2/titers_virus_strain.txt --strain-map results/strain-name-matching/h3n2/strain-name-map.tsv > /dev/null
```


```
cd data/titers/
aws s3 cp --exclude "*" --include "*.tsv.gz" --recursive s3://nextstrain-data-private/files/workflows/seasonal-flu/ .
```

```
for x in data/titers/*/*.tsv.gz; do
    echo ${x}
    gzcat ${x} | csvtk -t cut -f virus_strain | sort | uniq > ${x%.tsv.gz}_virus-strain.txt
    gzcat ${x} | csvtk -t cut -f serum_strain | sort | uniq > ${x%.tsv.gz}_serum-strain.txt
done
```

Following is too slow, so....
```
for x in data/titers/h3n2/*.txt; do
    echo ${x}
    ./scripts/strain-name-fuzzer.py \
        --curated-strains results/strain-name-matching/h3n2/metadata.txt \
        --query-strains ${x} \
        --strain-map results/strain-name-matching/h3n2/strain-name-map.tsv \
        > /dev/null
done
```

cat data/titers/h3n2/*.txt | sort | uniq > data/titers/h3n2.all-strains.txt

./scripts/strain-name-fuzzer.py \
    --curated-strains results/strain-name-matching/h3n2/metadata.txt \
    --query-strains data/titers/h3n2.all-strains.txt \
    --strain-map results/strain-name-matching/h3n2/strain-name-map.tsv \
    > /dev/null