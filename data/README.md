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
AF2-style floor applied to both predictors. **Decided 2026-11-16 — see below.**

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

### 2026-11-16: the noise floor is recycle convergence, scored on the core

**The floor definition (closes the open decision above).** The only mutation-free spread
that exists for *both* predictors is **recycle convergence**: refold the same sequence at
increasing recycle counts and measure the C-alpha RMSD between consecutive steps. AF2's
seed floor has no ESMFold counterpart (ESMFold has no seed), and AF2's model-to-model
spread is not a floor at all — five models are five different predictors. `pipeline/
noise_floor.py` computes all six frames and marks which are available.

Quote `rmsd_plateau_max` (the worst consecutive step in the second half of the sweep),
not `rmsd_mean`. The mean averages in the early steps, which are nowhere near converged,
and overstates the floor by roughly an order of magnitude.

**Answering Julie's open question above** (which column is "the" floor for the rendering
colour scale): use `rmsd_plateau_max` from `floors_core.csv`, per construct. Her 1e-15 A
observation is right and the 1e-6 A clamp is the right call — PDB coordinates carry three
decimals, so nothing below ~1.7e-3 A is measurable and anything smaller is superposition
rounding. Her `subset_v0` AF2 row (one seed pair, `ubq_human` only) is a seed floor, which
exists for AF2 but has no ESMFold counterpart, so it cannot be the shared scale; it stays
in the ladder as its own rung.

**The C-terminal tail contaminates every ubiquitin number.** `ubq_yeast` looked
non-convergent — late steps oscillating 0.05-0.20 A instead of decaying. It is not the
fold. At the worst late step (recycle 5->6, global RMSD 0.2012 A) residues 76 and 75 move
1.297 A and 1.124 A while residues 1-71 move 0.0186 A RMS. That is ubiquitin's flexible
C-terminal tail (72-76) swinging about 1 A between recycles *regardless of sequence*.
PDZ3 has no such feature and is unaffected.

So `noise_floor.py` takes `--exclude-tail N`: superposition still uses the whole chain,
only the RMSD scoring is restricted. **Use `--exclude-tail 5` for both ubiquitins, 0 for
PDZ3.** Two tables are kept side by side, `floors.csv` (global) and `floors_core.csv`
(core), so the contamination stays visible rather than being quietly cropped out.

| construct | core floor | global floor | median core mutation response | ratio |
|---|---|---|---|---|
| `ubq_yeast` | 0.0353 A | 0.2012 A | 0.1456 A | 4.1x |
| `pdz3_rat` | 0.0104 A | 0.3372 A | 0.1580 A | 15.2x |
| `ubq_human` | 0.0044 A | 0.1343 A | — | — |

ESMFold's *replicate* floor stays exactly 0.0000 A for all three (3 pairs each), as
before. `ubq_yeast`'s core steps still drift slightly upward after recycle 4
(0.0101 -> 0.0186 -> 0.0248 -> 0.0353), so 0.0353 A is a conservative drift amplitude,
not a settled residual — worth one line on the limitations slide.

**This is not just a floor correction; it changes the headline result.** Dropping the
tail sharpens the ubiquitin DMS correlation from rho = -0.254 (p = 0.012) to
**rho = -0.344 (p = 5.7e-4)** and leaves PDZ3 alone (-0.439 -> -0.432), exactly as
predicted. Across the 97 `ubq_yeast` mutants the median tail RMS is 1.1370 A against a
core median of 0.1456 A — the global RMSD we had been quoting for ubiquitin was mostly
tether, not fold. This is also why the 8 A shell displacement beat global RMSD for
ubiquitin earlier (rho = -0.482): the shell metric was already excluding the tail.

**ESMFold stays on the laptop.** Its folds are all CPU/bf16. `run_esmfold.py` skips any
output that already exists, so running the arm on a GPU node would fold only the gaps, at
a different device and precision, and put a hardware difference inside the replicate
floor. The recycle sweep runs locally too. The cluster runs the AF2 arm only.

**Hoffman2 scheduler.** The cluster moved to Slurm on 2026-11-05 and a Grid Engine GPU
job submitted after that date never dispatches. Do not pass `-p/--partition` or `--qos`:
`pi_kappel` owns nodes, so Slurm routes `-A pi_kappel` to `pi_shared`/`pi_shared24`
(24 h, 4 concurrent GPUs) by itself, and `-p campus` is rejected. See `hoffman2/README.md`.
