# Influenza A(H3N2) HA dataset based on reference "A/Wisconsin/67/2005"

| Key                  | Value                |
| -------------------- | -------------------- |
| authors                | [Richard Neher](https://neherlab.org), [Nextstrain](https://nextstrain.org)                         |
| name                 | Influenza A(H3N2) HA                      |
| reference            | A/Wisconsin/67/2005                      |
| dataset path         | flu/h3n2/ha/CY163680                     |
| reference accession  | CY163680   |
| clade definitions  |  [github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/](https://github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/)  |

## Scope of this dataset
This dataset uses an older reference sequence (A/Wisconsin/67/2005) and recent sequences will differ at a large number of positions from this reference.
For the analysis of currently circulating viruses, the dataset using A/Darwin/6/2021 as reference might be more appropriate.

## Features
This dataset supports

 * Assignment to clades and subclades based on the nomenclature defined in [github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/](https://github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/)
 * Identification of glycosilation motifs
 * Sequence QC
 * Phylogenetic placement

## Clades of seasonal influenza viruses

The WHO Collaborating centers define "clades" as genetic groups of viruses with signature mutations to facilitate discussion of circulating diversity of the viruses.
Clade demarcation do not always coincide with significantly different antigenic properties of the viruses.
Clade names are structured as _Number-Letter_ binomials (with exceptions) separated by periods as in `3C.2a1b.2a.2a.1a`. These sometimes get shortened by omission of leading binomials like `2a.1`.

In addition to these clades, "subclades" are defined to break down diversity at higher resolution and allow following the spread of different viral groups.
These follow a Pango-like nomenclature consisting of a letter followed by a numbers separated by periods as in `G.1.3.1`.
The leading letter is an alias of a previous name.
Details of the nomenclature system can be found at [github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/](https://github.com/influenza-clade-nomenclature/seasonal_A-H3N2_HA/).

## What is Nextclade dataset

Read more about Nextclade datasets in Nextclade documentation: https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html
