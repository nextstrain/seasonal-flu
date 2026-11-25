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
| founders used by the nextclade datasets | the same, downloaded to `nextclade/data/{lineage}/{segment}/founders.fasta` |
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

### This repo: branches are collapsed once, in tidy tree

The 3y builds of `profiles/nextstrain-public.yaml` (full size, ~3000 sequences per tree,
glycosylation, LBI, titer models and embeddings switched off) split large lineages into
many labelled pieces, e.g. vic NA `B.8` into 9 and H1N1pdm NA `C.5.3` into 4. The tidy
tree itself was fine: when a founder sits at zero distance from the node where tips of
the parent lineage branch off, `prune_short_branches` in `sanitize_trees` collapsed the
stem, and `refine --stochastic-resolve` then resolved the polytomy at random and buried
those parent lineage tips 5-10 nodes deep inside the lineage.

Now branches are only ever removed in `tidy_tree.py`, with the lineages in view:

- it drops the founder tips (unless `--keep-founders`), collapses internal branches
  shorter than `--collapse-threshold` substitutions (0.1, the cutoff TreeTime uses) and
  merges the single child nodes that grafting leaves where it removed a founder
- subtrees with fewer than 3 sequences, built without IQ-TREE, get the distance of each
  sequence to the founder as branch length. They had none, which counted as zero and
  collapsed the stems of lineages represented only by founders, e.g. H3N2 HA `J` into
  `G.1.3` in the 6m build
- `flag_outliers.py` and `sanitize_trees.py` get `--keep-short-branches`
- `refine` uses `--keep-polytomies --branch-length-inference input` for tidy tree
  builds. TreeTime otherwise reoptimizes the branch lengths and then prunes those below
  0.1 substitutions, which collapsed branches of up to 6 substitutions, including
  lineage nodes. TreeTime ladderizes the tree only on setup and after resolving
  polytomies, so these collapses had also left the trees unladderized. The time trees
  barely change: on H3N2 3y HA, `J` and `K` move by 0.01-0.02 years, the root of the
  old references by one year

Founder tips used to reach `sanitize_trees` and survive it when an HA and an NA lineage
share a name; only `refine` removed them, for lack of metadata.

`assign_subclades_from_tree.py` does not label a lineage without tips of its own (only
its founder represents it), and labels a lineage on the polytomy of its MRCA also when
that holds subtrees made only of tips of ancestral lineages, which keep their own label,
unless an ancestral lineage is labelled on the same node.

Results on the 3y, 2y and 6m builds of `nextstrain-public` (18 trees): every sequence
placed, no internal node lost after tidy tree other than by pruning tips, all trees
ladderized, every lineage a tip is assigned to has a branch label. `tidytree` takes 1-2
minutes per tree when run alone. The automatic root follows the reference window: A,
A.2 or A.1.1 for 3y, B.1.2 or A.2.2 for H3N2 2y and 6m, A.3 and A.1 for vic 6m; H1N1pdm
HA is rooted at `unassigned` in all windows because some 2018-2020 references predate
its subclades. Every lineage has a single label except:

- H1N1pdm NA `C.5` (3y) and H3N2 NA `B.4` (all three windows): the founder of a sibling
  lineage descends from theirs, so the sibling (`C.3`, `B.3`) branches off at their
  polytomy. `C.3` carries all three founder mutations of `C.5` (607G, 1049A, 1415T),
  `B.3` is the `B.4` founder plus 989A, and `C.5.3` reverts 607 and 1049 relative to
  `C.5`. The NA founders need review upstream; this is not a tree problem.
- H1N1pdm NA `B.1`, two old reference strains in two places.

Disagreements with `augur clades` come from the signature mutations: up to 89 tips in
H3N2 6m NA, mostly `B.4.2` here and `B.2` there.

Dry runs of every profile resolve, except for reasons unrelated to tidy tree: `gisaid`
needs conda (fine without it), `upload` lacks pb2 references (as on `master`), and
`nextclade` and `nextflu-public` have since been deleted on `master`. `allflu` and
`nextstrain` mix tidy tree with `augur tree` builds (yam, internal segments).

### Merge of `master` (2026-11-25)

`master` had moved on by about 250 commits. Relevant for tidy tree: nextclade no longer
runs in the phylo workflow, ingest adds its calls to the metadata, so `tidytree` now
reads `subclade_nextclade_{segment}` from the build's `metadata.tsv` (same lineages for
every sequence of vic 6m NA as from `nextclade.tsv.xz`). `master`'s URL tables for
`clades.tsv` and `subclades.tsv` are replaced by the nomenclature repos of this branch;
the new `h1n1` lineage, like yam, has its definitions in `config/` and is copied into
`nomenclature/`. All profiles point at `nomenclature/`. Every profile dry runs except
`full-trees`, which defines `download_mlr_json` twice on `master` too, and `upload`
(no pb2 references, as before).

## Next steps (as of 2026-11-25)

Everything up to here is committed on `WIP/tidy-tree`. The reduced 12y builds
(300 sequences each) of `profiles/nextstrain` and the full 3y, 2y and 6m builds of
`nextstrain-public` have been run before the merge of `master`, and every profile has
been dry run after it. In order:

1. **Review the NA founders** of H1N1pdm `C.5`/`C.5.3`/`C.3` and H3N2 `B.3`/`B.4` in the
   nomenclature repos (see above), and check the other NA repos for sibling founders
   that descend from one another.
2. **Run a build with embeddings switched on.** Glycosylation, LBI, the titer models,
   MLR fitness and LCR sums ran on the 3y builds after the merge of `master`; embeddings
   have not yet run on a tidy tree.
3. **Merge `founders` into `main`** in the nomenclature repos, then change the default of
   `nomenclature_branch` in the `Snakefile` to `main`.
4. Phase 3 and phase 2 below. For phase 3, also declare the scripts as
   inputs of their rules so that editing them reruns the rules, and turn the toy cases
   used so far (nested lineages, a single tip lineage, a misplaced tip, a collapsed stem)
   into tests of `assign_subclades_from_tree.py`.

### How to test

Data needed in `data/{h3n2,h1n1pdm,vic}/`: `metadata.tsv` as written by ingest, which
holds the nextclade calls as `subclade_nextclade_{ha,na}`, `ha.fasta`, `na.fasta` and the
titer tables of the builds. The nomenclature files are downloaded by the workflow into
`nomenclature/`.

Reduce the builds with an extra config file, e.g. `test.yaml` outside the repo (or
`config_local.yaml`, which is gitignored). Set `resolution` and
`--subsample-max-sequences` as needed:

```yaml
array-builds:
  nextstrain:
    patterns:
      lineage: [h3n2, h1n1pdm, vic]
      center: [cdc]
      passage: [cell]
      assay: [hi, fra]
      resolution: [12]
    subsamples:
      global:
        filters: "--group-by region year --subsample-max-sequences 300 --min-date {min_date} --include {include} --exclude {exclude} --exclude-where 'ha!=True' 'na!=True'"
        priorities: "titers"
```

The builds with titer data are `h3n2_cdc_cell_fra`, `h1n1pdm_cdc_cell_hi` and
`vic_cdc_cell_hi`. Targets go **before** `--configfile`, which takes a list:

```
T=$(for b in h3n2_cdc_cell_fra h1n1pdm_cdc_cell_hi vic_cdc_cell_hi; do for s in ha na; do
      echo auspice/${b}_12y_${s}.json; done; done)
snakemake $T --profile profiles/nextstrain \
    --configfile profiles/nextstrain/builds.yaml test.yaml
```

- Snakemake does not rerun a rule when only its script changed. After editing
  `tidy_tree.py` or `assign_subclades_from_tree.py` add
  `--forcerun tidytree subclades_from_tree`.
- `--config 'nomenclature={"branch": "main"}'` switches the nomenclature branch; the
  files are fetched again because the URL changed.

What to look at:

- `logs/tidytree_<build>_<segment>.txt`: the chosen root lineage and outgroup, then the
  sequence assignment: `N of M sequences placed, 0 dropped, 0 excluded` is the expected
  outcome with the automatic root
- `logs/subclades_from_tree_<build>_<segment>.txt`: lineages without a node of their own,
  and with `--config 'tidy_tree={"compare_subclades": true}'` the tip labels that
  disagree with `augur clades`
- the auspice JSONs: each lineage should carry its `Subclade` branch label once, except
  the lineages listed as without a node of their own

## To do

### Phase 1 leftovers (done)

The nextclade dataset workflow uses the upstream founders, see above.

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

### Phase 4 — clade labels from the tree structure (done)

Subclades of tidy tree builds come from the tree. `subclades_from_mutations` (`augur
clades` with `subclades.tsv`) no longer runs for them; `tidy_tree.compare_subclades:
true` brings it back and reports the disagreements in `logs/subclades_from_tree_*`.

The older clades are no longer labelled, in any build: the `clades` and `import_clades`
rules, `clades.json` and the `clades` build parameter are gone, as are the `Clade`
colorings and filters of the auspice configs, which now default to the subclade (to the
region where a build has no subclades). Recent clades correspond to subclades; the older
ones are only defined by mutations. Derived haplotypes, which are relative to a lineage,
are now computed only for builds with subclades.

To bring the older clades back, colour the tips by a metadata column rather than label
branches. Nextclade has them, but ingest does not pass them on yet: `short-clade`
matches the former labels for all H3N2 and all but 7 old H1N1pdm tips of the 3y builds,
`legacy-clade` for all B/Vic tips. The clade files in `config/yam`, `config/h1n1`,
`config/h3n2/ha/emerging_clades.tsv` and `profiles/nextflu-private/vic/ha/clades.tsv`
are unused now.

### Nextclade dataset workflow (done)

The founders come from the nomenclature repos (`founders` in `nextclade/config/
config_dict.yaml`), like the guide tree and the clade tables, all from the branch in
`nomenclature_branch`. `align` realigns them to the reference of each dataset; the
per-reference copies in `nextclade/config/*/`, `extract_clade_founders`,
`translate_clade_founders` and `scripts/extract_founder_sequences.py` are gone.
`nextclade/update_founders.sh` and the untracked `founder_sequences_new*.fasta` are
obsolete. New founders are made upstream when a lineage is designated.

The tidy tree chooses its root lineage automatically; `root` in the config is the
outgroup where the root lineage has no founder (`unassigned`). The founder tips stay in
the dataset trees (`--keep-founders`), and the tree is a divergence tree only.

The subclade (`clade` and its back-compat copy `subclade`) is labelled from the tree,
legacy and short clades and the proposals with `augur clades`. The lineage calls come
from the released dataset, which does not know lineages designated or renamed since; so
`assign_subclades_from_tree.py --refine-calls` takes the call of `augur clades` with the
current `subclades.tsv` where it is a sublineage of the original call (H3N2 NA: 286 tips
`B.4.2` -> `C`, the former `B.4.2.2`; H1N1pdm HA: 206 tips `D.3.1` -> `E` and its
sublineages). The decisions are logged in `build/*/*/*/clade_subclade.log`.

Compared with the released datasets (April 2026) on 4000-6000 sequences per dataset,
the subclade calls differ only by nomenclature updates (renames, new designations such
as `J.2.4.2`, `K.1`, `D.1.2`) and a few old strains; legacy clades agree for all H1N1pdm
and all but 13 old H3N2 sequences, QC is unchanged. The `b` datasets lack
`config/b/{outliers,reference_strains}.txt` and could not be built.

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
  where both exist. Ingest writes it to the metadata as `subclade_nextclade_{segment}`,
  which is what `tidytree` reads.
- There is no need to take lineage calls from the `align` step. Since the merge of
  `master`, which removed nextclade from the phylo workflow, they come from the
  metadata, as for every other nextclade annotation.
