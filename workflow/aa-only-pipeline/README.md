# Augur protein analysis development using seasonal-flu

This workflow was used to test and develop [Augur PR 1958](https://github.com/nextstrain/augur/pull/1958)

Provision some starting data via the canonical seasonal-flu workflow:

```sh
# working directory: seasonal-flu
snakemake --cores 4 --configfile profiles/nextstrain-public.yaml \
    -pf data/h3n2/ha.fasta data/h3n2/metadata.tsv data/h3n2/ha/nextclade.tsv.xz
```

Downsample this alignment for speed / dev purposes:

```sh
mkdir -p workflow/aa-only-pipeline/{nuc,HA1,HA2}
augur filter --metadata data/h3n2/metadata.tsv --sequences data/h3n2/ha.fasta \
    --min-date 2023-01-01 --group-by month --subsample-max-sequences 500 \
    --output-sequences workflow/aa-only-pipeline/results/nuc/sequences.fasta \
    --output-metadata workflow/aa-only-pipeline/results/metadata.tsv
```

## References

I used  `CY163680` (it gets pruned out of trees due to the clock filter). These are part of the repo in `./config`

## Snakemake workflow targets

`snakemake --cores 1 --snakefile Snakefile-AA-only -pf auspice/<TARGET>.json`

1. `auspice/nuc.json` Normal nucleotide workflow, inferring mutations in SigPep, HA1, HA2 based on the nucleotide tree
2. `auspice/HA1.json`, `auspice/HA2.json`. Generate tree based on AA seqs for that peptide. Infer mutations in said peptide across tree.
3. `auspice/HA12.json`. Concatenate HA1 & HA2 AA alignments, generate protein tree, infer mutations on HA1 & HA2 independently across tree.



## Development notes

  1. We have nucleotide data, but we want to infer the tree (including branch lengths) from protein(s),
    then infer ancestral DNA seqs/muts on tree (as well as AA), and view both AA & DNA in Auspice.
    
      * The only hard part of this is using proteins for refine. Maybe you also want to concatenate them, but that can be done via a separate python script and is beyond the scope of this PR
      * Reconstructing proteins is already possible (independent from nuc reconstruction) via `augur ancestral`
      * Auspice should just work - we have nucleotide data after all.
  
  2. We don't have nucleotide data. We have a single protein.
      
      * Augur tree/refine as per (1).
      * Modify `augur ancestral` to do an AA reconstruction without a nuc reconstruction
      * Export just the AA muts & annotation structure for the gene in `augur export v2` (already done?)
      * Auspice needs to work without the `nuc` annotation, but this should be easy-ish?
  
  3. We don't have nucleotide data. We have multiple proteins.
  
      * concatenate the proteins first if you want to for tree reconstruction resolution
      * Augur tree/refine as per (1)
      * `augur ancestral` as per (1) - ensure we allow multiple proteins
      * export as per (1)
      * Auspice needs a little more effort to work here, but it's _possible_


