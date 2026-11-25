# Tidy tree: status and plan

Working notes for putting the lineage guided tree ("tidy tree") on a stable footing.
Written 2026-11-24, updated 2026-11-25. The `founders` branches of the nomenclature
repos are pushed; this repo's `WIP/tidy-tree` is not.

Tidy tree builds a phylogeny by building a subtree per lineage of a hierarchical
nomenclature and stitching the subtrees together along a guide tree, using the founder
sequence of each lineage as the attachment point. This gives more stable trees than a
single global tree, especially when there is little historical context in the sample.
The long term goal is for this to become an augur command, and for the clade annotation
to come from the tree structure rather than from signature SNPs in `clades.tsv`.

## Where the pieces live

| what | where |
|---|---|
| tidy tree script | `scripts/tidy_tree.py` (this repo, branch `WIP/tidy-tree`) |
| phylo workflow rules | `workflow/snakemake_rules/core.smk` |
| nextclade dataset workflow | `nextclade/Snakefile`, `nextclade/rules/common.smk` |
| founders used by the phylo builds | `.auto-generated/founders.fasta` of the nomenclature repos, downloaded to `nomenclature/{lineage}/{segment}/founders.fasta` |
| founders used by the nextclade datasets | `nextclade/config/{lineage}/{segment}/{ref}/founder_sequences*.fasta` |
| nomenclature repos | `../clades/seasonal_*`, branch `founders` |
| nomenclature scripts | `../clades/helper-scripts`, branch `founders` |
| full sequence data used for testing | `../seasonal-flu-tidy-tree/nextclade/data` |

## Done

### Upstream: founders in the nomenclature repos (branch `founders` in all 7 repos)

Each subclade yaml now carries a `founder` block giving the nucleotide differences to
the founder of its parent, so any founder can be reconstructed from the reference
without access to sequence data:

```yaml
founder:
  method: representative-consensus      # or ancestral-inference
  provenance: "consensus of A/WesternAustralia/2534415989/2025, A/Georgia/38/2026"
  mutations:
  - {position: 300, state: 'A'}
```

Positions are 1-based in the coordinates of `config/reference.fasta`, the same frame as
the `nuc` locus of `defining_mutations`. Storing differences to the parent keeps the
entries short, is independent of the reference a downstream build aligns to, and makes
updates readable in a diff.

**181 founders across the six repos**, all validating: 52 H3N2 HA, 31 H1N1pdm HA, 30
H1N1pdm NA, 23 H3N2 NA, 26 B/Vic HA, 19 B/Vic NA. They were seeded from the ancestral
sequences inferred by the nextclade dataset workflow and reproduce those files exactly,
apart from one deliberate repair (G.1.3, see below).

Scripts in `../clades/helper-scripts`:

- `founder_utils.py` — shared yaml/reference/annotation handling, founder assembly,
  indel left alignment, repair of founders missing their own defining mutations
- `assemble_founders.py` — build and validate founders, write `.auto-generated/founders.fasta`
  and `.auto-generated/representatives.tsv`; `--strict` fails CI
- `compute_founders.py` — derive a founder from representative sequences
- `seed_founders.py` — one-off backfill from existing full length founder sequences

Each repo has `.github/workflows/validate-founders.yml`: it checks out the sister
`helper-scripts` repo, runs `assemble_founders.py --strict`, and verifies the
`.auto-generated` files are up to date. Public data only, no secrets.

### This repo: tidy tree phase 0

Two commits on `WIP/tidy-tree`:

- `64191585` tidy-tree: account for every sequence instead of dropping some silently
- `bdf51121` phylo: choose the tree method per build instead of always using tidy tree

**Silent sequence loss.** Sequences were discarded without notice when their lineage
call was not a node of the guide tree (newly designated lineages, finer resolution than
the guide tree, anything above `--root-lineage`) and when a lineage of the guide tree
had no founder, which with `--ignore-missing-founders` removed the lineage together with
its sequences and all of its descendants. Now an unknown lineage is placed in the
nearest enclosing lineage or the root (`--unknown-lineage-policy`), a lineage without a
founder is folded into its parent keeping its sequences and re-attaching its children,
every sequence is reported in exactly one category, and `--max-dropped-fraction`
(default 0.05) fails the build. On a test of 24 sequences spanning a founderless lineage
and a lineage absent from the guide tree, the tree went from retaining 8 to retaining 24.

**Per build tree method.** `uses_tidy_tree(build_name, segment)` in `core.smk` requires
both a guide tree and founders. Builds without a nomenclature (internal segments, yam,
h1n1) go back through `augur tree` and `prune_reference`, which had been commented out,
and keep `--include-reference` in their alignment. Builds with one get the founders in
their alignment and go through `tidytree`, which writes `tree_tidy.nwk` so the two rules
are never ambiguous; `prune_outliers` takes whichever tree the build produced.
`tidy_tree.enabled: false` turns it off globally. Before this, a dry run of
`profiles/allflu` failed with `KeyError: 'yam'`; it now resolves to 6 `tidytree` jobs
and 2 `tree` plus 2 `prune_reference` jobs.

Also fixed: the stray `temp_tree.nwk` written into the working directory, and an
unbound variable crash when a lineage ends up with no sequences.

### This repo: tidy tree phase 1, phylo builds use the upstream founders

The clades, subclades, guide tree and founders all come from one nomenclature repo per
lineage and segment (`nomenclature_repo_by_lineage_and_segment` and `nomenclature_url`
in the `Snakefile`), downloaded into `nomenclature/{lineage}/{segment}/`, which is not
version controlled. Lineages without a repo (yam) copy their frozen definitions from
`config/` into the same place, so all profiles point at `nomenclature/`. The branch is
configurable via `nomenclature: {branch: ...}` and defaults to `founders` until that is
merged, then it should become `main`. Snakemake records the download URL, so changing
the branch re-downloads.
`download_founders` strips the alignment gaps and `align` realigns the founders to the
build reference with `--gap-alignment-side right`. `founder_refs` is gone. The download
rules use `curl -fsSL`, so a missing file fails instead of writing a 404 page.

`tidy_tree.py` now recognises founders against the full guide tree. A founder of a
lineage outside the root lineage used to be counted as a sequence without a lineage
call, i.e. dropped; it is now reported as unused.

A root lineage without a founder, such as `unassigned`, is rooted on an outgroup:
`--explicit-root` if given, otherwise the oldest exactly dated sequence of that lineage
from `--metadata`.

### This repo: root lineage and subclade labels from the tree

The root lineage is chosen automatically (`--root-lineage auto`, the default) as the
most specific lineage of the guide tree containing the lineages of all sequences. The
per-resolution table in the `Snakefile` is gone; a root can still be set with
`tidy_tree.root_lineage.{lineage}.{segment}.{resolution}`, and sequences outside it are
then excluded and reported, not counted as dropped. Previously anything outside the root
lineage was built into the root lineage's subtree as if it belonged to it. In practice
every build contains old reference strains forced in by the subsampling (H3N2 back to
1968), so the root is `unassigned`, rooted on the oldest of them, and the old diversity
sits at the base of every tree.

Subclade labels for tidy tree builds come from `scripts/assign_subclades_from_tree.py`
instead of `augur clades`. `tidy_tree.py --output-lineages` writes the lineage of each
sequence as a node of the guide tree (call, else enclosing lineage, else root). For each
lineage, from the root of the guide tree down, the label goes on the MRCA of its tips
(including sublineages) and the subtree below is assigned to it. If the MRCA subtree also
holds tips of other lineages, the lineage has no node of its own and its largest pure
subtrees are labelled instead. This happens when the stem of a lineage has zero length
and `prune_short_branches` in `sanitize_trees` collapses it; plain MRCA labelling let
H1N1pdm NA `B.1` take over most of the tree. `augur clades` still runs alongside
(`subclades_from_mutations`) and the disagreements are in
`logs/subclades_from_tree_*.txt`.

Tested on the 12y `profiles/nextstrain` builds for h3n2, h1n1pdm and vic, HA and NA,
subsampled to 300 sequences: all six run to auspice JSON with every sequence placed.
Tree labels agree with nextclade for every tip (largely by construction), `augur clades`
for 89% (H3N2 NA) to 100%. Every lineage gets a single branch label except four with a
collapsed stem: H1N1pdm NA `B.1` (4 branches), `C.5`, `D.2` and B/Vic NA `B.2` (2 each).

## To do

### Phase 1 leftovers

- The nextclade dataset workflow still uses the per-reference copies in
  `nextclade/config/*/`, including the per-CDS files it needs for `founderMuts`. Move it
  to the upstream founders too, which retires `nextclade/update_founders.sh` and the
  hand renaming of `founder_sequences_new*.fasta`. The `extract_clade_founders` rule
  should become a step that proposes an updated `founder` block upstream, ideally as a
  PR, instead of writing files that are moved by hand.

### Phase 2 — founder lifecycle

A newly designated lineage has no founder until a tree has been built with it, which is
the chicken and egg problem. Two sources:

- **inferred** (`method: ancestral-inference`): the current loop, tidy tree then
  `augur ancestral` then the sequence at the lineage root node
- **representative consensus** (`method: representative-consensus`): computed from the
  representatives listed in the yaml, available as soon as the lineage is named

Where both could be computed they agree exactly (see the cross check below), so the
representative route is sound for new designations. Requires at least two
representatives that **cover the position and sample the basal diversity** of the
lineage.

### Phase 3 — refactor `scripts/tidy_tree.py`

Split into a library and a CLI wrapper. Replace the recursion with loops so the
`sys.setrecursionlimit(10000)` workaround can go. Reuse `augur tree`'s IQ-TREE runner
and name escaping rather than the local `run_iqtree`. Build independent subtrees in
parallel. Add small synthetic test cases. (The pre-nomenclature context is now handled
by rooting at the guide tree root on its oldest sequence, which replaces the planned
`__context__` pseudo lineage.)

### Phase 4 — clade labels from the tree structure (done, in transition)

Remaining: drop `subclades_from_mutations` for tidy tree builds once the comparison has
been watched for a while. HA clade labels (`clades.json`) still come from `augur clades`
with `clades.tsv` and are copied to other segments by `import_clades`.

### Phase 5 — augur

Either `augur tree --method iqtree` with guide options (`--lineage-guide-tree`,
`--lineages` plus a column, `--founders`, `--root-lineage auto`) or a separate
`augur tidy-tree`, with the phase 4 labels as an extra `--output-node-data`.

## Open decisions

1. Should the generated `founders.fasta` be committed or built in CI? (The `founder`
   block itself will be accepted upstream.)
2. Should founder tips stay in the nextclade reference trees? Useful as placement
   targets, noise in auspice.
3. Should synthesised founders for brand new lineages be generated upstream at
   designation time or here at build time?
4. Should `left_align_gaps` stay in `compute_founders.py` now that the nextclade flag is
   documented? It is currently a safety net that warns and corrects; the alternative is
   to trust the flag and delete it (two line removal).

## Gotchas worth not rediscovering

**Gap alignment side.** A deletion in a repeat has several equally parsimonious
placements. nextclade defaults to `--gap-alignment-side left`, which places it further
3' than the convention the founders follow. Always align with
`--gap-alignment-side right`, which is what the seasonal-flu workflows already use and
which puts the B/Vic HA1 162-164 deletion at nucleotides 562-570 as the community
expects. Getting this wrong made all six B/Vic HA subclades look as if they differed
from their founders at twelve positions. `assemble_founders.py --strict` now enforces
the convention and all 181 founders satisfy it.

**Representatives must be basal, not just numerous.** Unanimity among representatives
that all sit on one branch within the lineage puts that branch's mutations into the
founder. H1N1pdm NA D.2.1 had three representatives that agreed on four mutations the
founder should not have. The fix is a representative that branches early, not more
representatives.

**A representative must actually be a member.** A/Albania/317901/2025 carries NA:46V and
NA:95N but NA:343E, so it lacks D.2.1's defining constellation and sits on the branch
leading to the lineage rather than within it. Using it pulled the defining NA:343K out
of the founder; `compute_founders.py` refused to write and the defining mutation check
is what caught it. It was replaced by A/Tokyo/EIS14-800/2026, after which the computed
founder matches the inferred one exactly.

**`--strict` alone is not enough.** The defining mutation check only covers defining
mutations, so removing any other founder mutation passes it. The CI step that diffs the
regenerated `.auto-generated` files against the committed ones is the real guard.

**Ancestral inference can place a mutation one node too late.** H3N2 HA G.1.3's inferred
founder lacked HA1:192F although all 18 of its descendants carry it.
`seed_founders.py --enforce-defining-mutations` repairs this using the codon of the
descendants and reports what it changed. Only G.1.3 needed it across all six repos.

**Founders are in the alignment but not in the metadata.** When accounting for sequences
in `tidy_tree.py` they must be excluded explicitly, otherwise they count as dropped; at
52 of 352 that alone would trip `--max-dropped-fraction`.

## Cross check: representative consensus vs ancestral inference

Computed for every subclade with at least two available representatives, after left
aligning indels. **15 of 23 reproduce the inferred founder exactly.**

| repo | exact | remaining |
|---|---|---|
| H3N2 HA | 6/6 | — |
| H1N1pdm HA | 5/7 | D.2 (2 reps disagree), D.3.1 (terminal position only) |
| H1N1pdm NA | 2/2 | — |
| H3N2 NA | 1/2 | C.2 (last base of the reference) |
| B/Vic HA | 1/6 | five subclades with 2 representatives each |

Everything still outstanding is a subclade with only two representatives that may not
sample basal diversity, the pattern described above. B/Vic NA has no representatives in
any of its 19 yamls and could not be checked at all; curating those is the obvious next
piece of nomenclature work.

Reproduce with:

```
cd ../clades/<repo>
nextclade run --gap-alignment-side right --input-ref config/reference.fasta \
    --input-annotation config/annotation.gff \
    --output-fasta /tmp/reps.aligned.fasta /tmp/reps.fasta
python3 ../helper-scripts/compute_founders.py --input-dir subclades \
    --reference config/reference.fasta --annotation config/annotation.gff \
    --sequences /tmp/reps.aligned.fasta --dry-run
```

## Corrections to earlier analysis

Two findings from the first review of this branch were wrong and should not be chased
again:

- The phylo workflow does **not** use the wrong lineage column. The `clade` column of
  the nextclade output holds subclade names in every lineage, identical to `subclade`
  where both exist, so `--lineage-column clade` is correct.
- `data/{lineage}/{segment}/nextclade.tsv.xz` **is** produced by the workflow, by
  `run_nextclade` in `workflow/snakemake_rules/select_strains.smk`, with
  `ruleorder: download_nextclade > run_nextclade` preferring the copy on S3. There is no
  need to take lineage calls from the `align` step; the existing route covers the full
  dataset rather than just the build's subsample.
