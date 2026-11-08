# Disagreement data analysis

**Your role:** build the analysis that answers the project's central question. For the
same point mutation, does ESMFold (no MSA) move the structure more than AlphaFold2 (MSA),
and is that motion real or noise?

Everything you need is in `inputs/`. The shared code is in `src/` at the repo root. Background
and hypotheses are in the repo-root `README.md`; read the HYPOTHESES and REFERENCE FRAMES sections first.

---

## 1. The data you have, and what it is not

`inputs/` is a **development subset** (`subset_v0`), folded on a laptop on 2026-11-08.
Use it to write and debug your code. **Do not report any number from it.** It used
reduced AF2 settings (one model, shallow MSA) to fit on a CPU.

The real panel (`full_v1`: about 100 stratified mutants per protein, all 5 AF2 models, 5
WT seeds) is folding on Hoffman2. It will arrive with **the same files, columns and
layout**, so code that runs on `inputs/` runs on it unchanged. Write nothing that hardcodes
a variant name, a row count or a protein list.

### What was folded

| protein | what it is | why it's here |
|---|---|---|
| `ubq_yeast` | yeast ubiquitin 1-76 | DMS arm (Roscoe 2013 fitness) |
| `pdz3_rat` | rat PSD95 PDZ3, residues 302-402 | DMS arm (McLaughlin 2012) + crystal pairs |
| `ubq_human` | human ubiquitin 1-76 | crystal pairs only (the crystals are human) |

Mutants are sampled on two axes: **fitness** (deleterious = bottom quartile of the DMS
score, neutral = at or above the median) × **burial** (buried = relative SASA < 0.20 on the
crystal WT). Crystallized mutants are tagged `set=experimental_pair`.

Coverage in the subset: ESMFold folded all 16 sequences. AF2 folded the WTs plus I23N,
K48R, A375D, G330T, H372A, L67S, K33S, K11S. `pairs.csv` lists only pairs where both
structures exist.

### Files

```
inputs/tables/structures.csv     one row per predicted structure (34)
inputs/tables/per_residue.csv    one row per residue per structure: pLDDT + C-alpha xyz
inputs/tables/pairs.csv          the comparisons to make (mut_vs_wt, wt_vs_wt)
inputs/tables/variants.csv       mutant metadata: stratum, RSA, DMS score, crystal IDs
inputs/tables/sequences.csv      every folded sequence (+ sha1)
inputs/tables/constructs.csv     the 3 constructs: sequence, construct_start, notes
inputs/preds/esmfold/*.pdb       ESMFold structures (+ .json metadata: pTM, settings)
inputs/preds/colabfold/*.pdb     AF2 structures
inputs/experimental/*.pdb        the 7 crystal structures, trimmed to the construct
inputs/experimental/experimental_structures.csv   which crystal is which
```

`structures.csv:path` is relative to `inputs/`, e.g. `inputs/` + `preds/esmfold/pdz3_rat__WT__rep1.pdb`.

### Key columns

`structures.csv`

| column | meaning |
|---|---|
| `structure_id` | `seq_id\|predictor\|model\|s<seed>\|r<replicate>`; the join key |
| `protein`, `variant` | e.g. `pdz3_rat`, `G330T` (or `WT`) |
| `predictor` | `esmfold` or `alphafold2` |
| `model`, `seed`, `replicate` | AF2 model/seed; ESMFold replicate number |
| `set`, `fitness`, `burial`, `rsa`, `DMS_score` | sampling stratum and DMS data |
| `pos_uniprot` | mutated residue, UniProt numbering (empty for WT) |
| `mean_plddt`, `site_plddt`, `ptm`, `max_pae` | model confidence (pLDDT 0-100) |
| `settings` | exact fold settings; check this before mixing runs |

`per_residue.csv`: `structure_id, resid_construct, resid_uniprot, aa, plddt, ca_x, ca_y, ca_z`.
Coordinates are in each model's own frame, **not superposed**.

`pairs.csv`: `kind` (`mut_vs_wt` / `wt_vs_wt`), `structure_a` (always the reference / WT),
`structure_b`, `predictor`, `model`, `seed` (`-1` = the two differ in seed), `pos_uniprot`.
Mutant-vs-WT pairs are already like-for-like: same predictor, model and seed.

---

## 2. Numbering: the one thing that will silently break your results

There are two residue numberings in play:

- **construct numbering**: 1..L, as the predictor wrote it. **The predicted PDBs use this.**
- **UniProt numbering**: what DMS tokens (`E373V` = residue 373) and the crystal PDBs use.

`resid_uniprot = resid_construct + construct_start - 1`. For ubiquitin `construct_start = 1`,
so the two numberings coincide. For PDZ3 `construct_start = 302`: construct residue 1 is
UniProt 302, and **the mutated site E373V is residue 72 inside the predicted PDB**.

Consequences:

- `geometry.locality_profile(..., site=...)` on predicted PDBs needs the **construct** number:
  `site = pos_uniprot - construct_start + 1`.
- `geometry.superpose` matches residues **by number**. Predicted PDZ3 (1-101) against a
  crystal (302-402) has nothing in common and raises an error, which is the good case.
  Renumber the prediction first (add 301) before comparing it to a crystal.
- In the tables, always join on `resid_uniprot`. Never join on row order.

---

## 3. Your task, in order

Each step produces a table. Save outputs as CSV under `disagreement_data_analysis/outputs/`
(create it; it is git-ignored, so share deliverables with the team explicitly).

**Step 1: displacement per pair.** For every row in `pairs.csv`, superpose B onto A on
C-alpha and record global RMSD, per-residue displacement and per-residue ΔpLDDT (B − A).
`src/geometry.py` already does this: `load_ca`, `superpose`. Output a per-pair summary
(RMSD, max displacement, mean ΔpLDDT, site ΔpLDDT) and a long per-residue table.

**Step 2: the noise floor.** From the `wt_vs_wt` pairs, per protein × predictor: the
distribution of RMSD and per-residue displacement (`geometry.noise_floor`). A mutant
effect counts only if it clears this.
> **Known:** ESMFold is deterministic. Its repeat folds are bit-identical, so its floor
> is exactly 0. Your code must handle a zero floor (no division by 0, no log of 0). The
> lead is deciding the ESMFold-appropriate floor; for now, report both predictors' floors
> as measured and use the AF2 seed floor as the common reference.

**Step 3: locality.** For each mutant: mean displacement inside vs outside an 8 Å shell
around the mutated residue, measured on the WT structure (`geometry.locality_profile`).
A ratio near 1 = motion spread over the chain (noise-like). Well above 1 = a local response.

**Step 4: the disagreement (the headline).** Pair ESMFold and AF2 **by variant**.
For each variant: ESMFold displacement minus AF2 displacement (global RMSD, shell
displacement, site ΔpLDDT), each expressed relative to that predictor's noise floor.
Also measure **WT_ESMFold vs WT_AF2** for each protein: how different the two predictors
are before any mutation.

**Step 5: contact-map change.** C-alpha contacts at 8 Å, WT vs mutant per predictor:
fraction of contacts lost/gained, overall and within the shell. Packing changes can leave
RMSD flat but change contacts.

**Step 6: stratify.** Break steps 1-5 down by fitness × burial. Bootstrap 95% CIs.
**Per protein, never pooled across proteins.** Report n per stratum next to every number.

**Step 7: the experimental reference frames** (`inputs/experimental/`):
- *accuracy floor*: predicted WT vs crystal WT, per predictor (renumber first, see §2).
- *ground-truth deformation*: crystal WT vs crystal mutant. Then compare predicted
  displacement to real displacement for the same mutation.
- Crystal chains can be missing a terminal residue (5OXI and 4PIG chain A end at 75).
  `superpose` matches by number and handles it; report `n_matched`.

**Deliverable for the rendering member.** For each 6-way pair, write
`outputs/displacement_<protein>_<variant>.csv` with columns:

```
object,resid_uniprot,displacement_A
exp_mut,1,<float>
af_mut,1,<float>
esm_mut,1,<float>
...
```

`object` is one of `exp_mut`, `af_mut`, `esm_mut`: that structure's displacement against
its own WT (exp_mut vs exp_wt, af_mut vs af_wt, esm_mut vs esm_wt), in UniProt
numbering. The rendering member colors the 3D figures with this exact format, so don't change it.
Also give them each predictor's noise-floor value (one number per protein × predictor)
in `outputs/noise_floor.csv`; it sets the white point of the color scale.

---

## 4. Setup

```bash
conda create -n cair -c conda-forge python=3.11 numpy pandas biopython scipy matplotlib jupyter
conda activate cair
cd julie            # work from your folder; your notebook is julie/julie.ipynb
jupyter lab julie.ipynb
```

In the notebook (paths are relative to `julie/`):

```python
import os
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"  # before numpy
import sys; sys.path.append("../src")          # repo-root src/
from geometry import load_ca, superpose, locality_profile, noise_floor
import pandas as pd
IN = "disagreement_data_analysis/inputs/"
st = pd.read_csv(IN + "tables/structures.csv").set_index("structure_id")
pairs = pd.read_csv(IN + "tables/pairs.csv")
a = load_ca(IN + st.loc[pairs.structure_a[0], "path"])
```

Pin BLAS to one thread (the `OMP_NUM_THREADS` line). Small-matrix numpy work is otherwise
thousands of times slower on a laptop.

Code that more than one person will use goes in the repo-root `src/`. Analysis goes in your
notebook (`julie/julie.ipynb`). Save outputs to `disagreement_data_analysis/outputs/`
(git-ignored). Findings (row counts, surprises, decisions) go in the repo-root
`data/README.md` in the same commit that found them.

---

## 5. Rules of the road

- **A null result is a result.** If neither predictor clears its noise floor, that's the
  finding. Nobody needs to go looking for a positive.
- Global RMSD is reported for comparability with prior work; it's not the finding.
  Per-residue displacement and locality are.
- pLDDT is the model's confidence, not stability. A pLDDT drop doesn't mean the
  mutation destabilizes. Compare ESMFold and AF2 pLDDT only as within-predictor
  *differences*: they come from differently trained heads.
- McLaughlin 2012 (PDZ3) measures **peptide binding**, not folding. A "deleterious" PDZ3
  variant may be perfectly folded.
- Before claiming a trend from a plot, check it against the numbers. A plot that looks
  like it shows something is not evidence.

## 6. Provenance of these inputs

- **ESMFold:** HF `facebook/esmfold_v1`, CPU, ESM-2 LM in bf16, 4 recycles, WT × 3.
  Per-residue pLDDT is rounded to 1 unit in these PDBs (fixed for full_v1).
- **AF2:** ColabFold 1.5.5 `alphafold2_ptm` model_1, 3 recycles, no templates,
  `--max-msa 512:1024`. The WT MSA is reused for every mutant, with the query row substituted.
- **Crystals:** RCSB, chain A, altloc A, protein atoms only, trimmed to the construct.
  1UBQ (WT), 5OXI L67S, 4PIH K33S, 4PIG K11S (human ubiquitin); 5HDY (WT), 5HET G330T,
  5HF4 H372A (rat PDZ3; 297-301 and 403-415 are expression tags, removed).
- Every mutant sequence was checked in code to carry the intended residue at the
  intended position (`src/mutations.py:assert_variant_set`) before folding.
