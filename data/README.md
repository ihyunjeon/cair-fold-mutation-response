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

### 2026-11-08: constructs, numbering, first fold runs

**Where the data lives.** Raw downloads, weights and run directories are on the lead's
external SSD at `/Volumes/LaCie/data/cair_fold/` (`raw/`, `weights/`, `subset_v0/`,
`full_v1/`). Members get a zipped copy of a run directory and point the code at it.

**Assays (ProteinGym v1.3).** `RL40A_YEAST_Roscoe_2013` (1195 single mutants) and
`DLG4_RAT_McLaughlin_2012` (1576 single mutants). Both ProteinGym target sequences are
*full-length precursors*, not the folded domain:

| protein | ProteinGym target | folded construct | `offset` for `apply_variant` on the construct |
|---|---|---|---|
| ubiquitin | RL40A_YEAST, 128 aa (ubiquitin-L40 fusion) | residues 1-76 | 0 |
| PSD95 PDZ3 | DLG4_RAT, 724 aa | residues 302-402 (101 aa) | -301 |

DMS positions are used unchanged (UniProt numbering); the construct string is indexed
with the offset above. `assert_variant_set` passes on every variant (subset 4/4, 6/6, 3/3;
full panel 97/97, 101/101, 3/3).

**PDZ3 crystal construct.** 5HDY/5HET/5HF4 carry expression tags at 297-301 (GSPEF) and
403-415 (NSRVDSSGRIVTD). These are *not* DLG4 sequence; we fold only 302-402.

**Yeast vs human ubiquitin.** The DMS is yeast; every ubiquitin crystal pair is human. They
differ at P19S, E24D, A28S. So the DMS arm folds yeast 1-76 (`ubq_yeast`) and the
experimental 6-way arm folds human 1-76 (`ubq_human`). Burial for yeast variants is taken
from human 1UBQ; the three differing positions are on the surface.

**Burial.** Relative SASA (Shrake-Rupley, Tien 2013 max ASA) on 1UBQ chain A and 5HDY
chain A, protein atoms only; buried = RSA < 0.20.

**Fitness strata.** Within each assay: deleterious = bottom quartile of DMS_score,
neutral = at or above the median, the middle band unused. McLaughlin 2012 is a
*peptide-binding* selection, not stability: a "deleterious" PDZ3 variant may be fully
folded and just a poor binder.

**Experimental pairs (RCSB, X-ray, single protein entity).**
ubiquitin 1UBQ vs L67S 5OXI ("C-terminally retracted", so a real backbone change is
expected), K33S 4PIH, K11S 4PIG (surface-Lys controls); PDZ3 5HDY vs G330T 5HET and
H372A 5HF4 (same lab, same construct).

**MSAs (AF2 arm).** WT MSA reused for every mutant, with only the query row substituted
(decided 2026-11-08). ColabFold MMseqs2 server depth: ubq_yeast 21413, ubq_human 21577,
pdz3_rat 12789 sequences.

**Traps found getting folds to run.**
- ESMFold weights are stored fp16 (7.2 GB). Loading them as fp32 needs ~14 GB of RAM and
  thrashed an 18 GB laptop to 28 GB of swap. Load in fp16/bf16 (`--lm-dtype`).
- ESMFold with an fp16 language model on Apple MPS aborts (`Destination NDArray and
  Accumulator NDArray cannot have different datatype`). Don't use MPS.
- ColabFold 1.5.5 needs `jax==0.4.23` + `dm-haiku==0.0.10`. Newer jax removes
  `jax.linear_util` and haiku fails on import.
- Hoffman2 nodes are CentOS 7 / glibc 2.17. Current scipy/torch wheels need glibc 2.28,
  so pip falls back to source builds that fail. `hoffman2/setup_envs.sh` pins versions
  that ship manylinux2014 wheels.
- The public ESM Atlas fold API returned 504s on uncached sequences. Not usable.

**ESMFold is deterministic.** Three repeat folds of each WT on the same machine gave
bit-identical coordinates (max C-alpha difference 0.000 A, 9 pairs, subset_v0). The
"replicate folds" noise floor is therefore exactly 0 for ESMFold and is only meaningful
for AF2 (seed/model spread). The noise-floor milestone needs an ESMFold-appropriate
definition before week 7 is called done: for example, spread over recycle count, or the
AF2-style floor applied to both predictors. **Open decision for the lead.**

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
