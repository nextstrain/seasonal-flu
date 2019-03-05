# augur builds for seattle flu study

This is the (in development) [augur][] build for understanding influenza dynamics in Seattle.
It is based off [this nextstrain build](https://github.com/nextstrain/seasonal-flu).

## Reassortment

#### Aim
Given phylogenies for each segement, we aim to identify sets of tips which have not undergone reassortment.
This allows the concatenation of the segments for each set of (non-reassorting) tips, allowing the construction of a more informative phylogeny.


#### Method
In development


## Running
This build starts by pulling sequences from our live [fauna][] database (a RethinkDB instance).
This requires environment variables `FAUNA_PATH`, `RETHINK_HOST` and `RETHINK_AUTH_KEY` to be set.

```bash
snakemake
```



[Nextstrain]: https://nextstrain.org
[fauna]: https://github.com/nextstrain/fauna
[augur]: https://github.com/nextstrain/augur
[auspice]: https://github.com/nextstrain/auspice

