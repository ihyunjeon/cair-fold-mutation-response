# data/

Nothing in here is committed except this file and the directory placeholders.
Structures, alignments and DMS tables are all regenerated or re-downloaded.

## Layout

| path | what |
|---|---|
| `raw/` | Downloads as received: ProteinGym DMS tables, reference PDB/AlphaFold DB structures. Never edited in place. |
| `processed/` | Anything we generate: mutant FASTA sets, predicted structures, geometry tables. |

## The lab-notebook rule

**This file is the project's memory.** When you discover a row count, a numbering
offset, a trap, or a reason something is the way it is, write it here in the same
commit as the code that found it. A finding that lives only in a notebook cell is
lost the moment someone restarts a kernel.

## Findings

### Noise floor (TASK Step 2) — Julie, 2026-11-13, dev subset `subset_v0`

- `pairs.csv` in `subset_v0` has 31 pairs: 21 `mut_vs_wt` and 10 `wt_vs_wt`.
- ESMFold WT replicates are identical, but `superpose` returns RMSD and
  displacement around 1e-15 Å, not exactly 0. This is floating-point rounding in
  the superposition, not motion. Values below 1e-6 Å are set to 0 in
  `noise_floor.csv`. PDB coordinates have 3 decimals, so no real difference can
  be that small.
- AF2 has a noise floor only for `ubq_human`, from one pair (seed 0 vs seed 1).
  `pdz3_rat` and `ubq_yeast` have a single AF2 WT fold each, so they have no
  floor and no row in `noise_floor.csv`.
- The floor is computed by grouping the Step 1 tables by protein × predictor,
  not by calling `geometry.noise_floor`, because `pairs.csv` already defines
  which WT folds are compared (same model, different seed).
- Open question: which column of `noise_floor.csv` is the floor value for the
  rendering color scale.
