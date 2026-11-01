# Does a structure predictor notice a mutation?

**CAIR Collective — Structural Biology / Molecular Biophysics team, 2026–27**
Project lead: Ihyun Jeon · Members: Julie, Yash, Kiran

---

## Research question

AlphaFold2 returns a near-wild-type structure when you mutate a single residue, even
for substitutions known to destabilize the protein. The usual explanation is that the
model leans on a multiple sequence alignment of thousands of homologs, so one changed
position gets outvoted.

**If that explanation is right, a predictor with no alignment should behave differently.**

ESMFold folds from a single sequence through a protein language model, with no MSA at
all. So:

> When you mutate one residue, does a single-sequence structure predictor move where an
> MSA-based one will not — and if it moves, does the motion mean anything?

## Hypotheses, stated so they can fail

| | Statement | What we would see |
|---|---|---|
| **H1** | MSA dominance explains the insensitivity | ESMFold deforms measurably more than AlphaFold2 for the same substitution |
| **H0** | The insensitivity is intrinsic to sequence-based structure prediction | Neither predictor exceeds its own run-to-run noise floor |
| **H2** | Any motion is physically meaningful | Displacement concentrates near the mutated site and scales with measured fitness cost |
| **H2-null** | Motion is numerical noise | Displacement is delocalized and unrelated to fitness |

A null result is a result here — provided the noise floor below is measured.

## Approach

**Panel.** Small, deeply characterized proteins with published deep mutational scanning
data: ubiquitin (76 aa), PSD95 PDZ3 (~83 aa), avGFP (~238 aa), and TEM-1 β-lactamase
(~286 aa) as a stretch. Short chains so folding is seconds, and structures an audience
recognizes.

**Sampling.** ~100 substitutions per protein, stratified on two axes: measured fitness
(neutral vs. strongly deleterious) and burial (buried vs. exposed). Random sampling
would only let us report an average; stratifying is what makes "does the response depend
on structural context?" answerable.

**Pipeline.**

| # | Step | Tool | Owner |
|---|---|---|---|
| 1 | Build mutant sequence set, stratified | Python, ProteinGym | Lead |
| 2 | Assert every mutant carries the intended substitution | `src/mutations.py` | Lead |
| 3 | Fold wild type **with replicates** | ESMFold + ColabFold (GPU) | Lead |
| 4 | Fold all mutants, both predictors | ESMFold + ColabFold (GPU) | Lead |
| 5 | Superpose, per-residue displacement | `src/geometry.py` | Member |
| 6 | Contact maps, local SASA, ΔpLDDT | BioPython | Member |
| 7 | Render overlays colored by displacement | PyMOL | Member |
| 8 | Paired comparison, locality, stratification | scipy | Lead + member |

Only steps 3–4 need a GPU, and those belong to the lead. Everything else runs on a
laptop against structure files of a few hundred KB. **No member's week can be blocked by
a cluster queue** — that is deliberate.

## The two controls that carry the project

**1. The noise floor.** Structure predictors are not deterministic across seeds and
recycles. Before any mutant is interpreted, we fold the *same wild-type sequence*
repeatedly and measure the spread in RMSD, per-residue displacement and pLDDT. A mutant
RMSD of 0.3 Å means nothing if repeated wild-type folds differ by 0.4 Å. We have not
found this floor reported in the papers claiming insensitivity.

**2. The substitution assertion.** Every generated mutant sequence is verified in code to
carry the intended residue at the intended position, with the wild-type residue matching
the DMS table. DMS, UniProt and PDB numbering disagree routinely. This class of error
never crashes — it silently produces a full set of plausible structures answering a
question nobody asked.

## Evaluation

Per-protein, never pooled across proteins. Paired by substitution when comparing
predictors. Bootstrap confidence intervals. Every claimed effect is reported against the
noise floor from control 1.

## Fall timeline (weeks 6–10)

| Week | Milestone |
|---|---|
| 6 | Panel frozen, data verified, sampling designed, environments working |
| 7 | Wild-type folds and replicates → **the noise floor** |
| 8 | ESMFold mutants for two proteins; geometry pipeline end to end |
| 9 | AlphaFold2 arm; paired comparison; locality and burial stratification |
| 10 | Figures, slides, 7-minute talk by all members |

## Winter extension

The causal test: if MSA depth is what makes AlphaFold2 insensitive, then progressively
**subsampling the alignment should restore sensitivity**, scaling with depth. That turns
the project from a comparison of two tools into a statement about why either behaves the
way it does.

## Repository layout

```
├── README.md          this file
├── notebooks/         one per member — see below
├── data/              raw + processed; contents gitignored (see data/README.md)
└── src/               shared scripts and utilities
```

**Notebooks.** Coding members use code cells for code and markdown cells for notes.
Non-coding members use markdown cells throughout, as a digital lab notebook, and must
complete at least one literature review.

**`data/README.md` is the project's memory.** Row counts, numbering offsets, traps and
reasons go there in the same commit as the code that found them.

## Status

Week 6. Nothing measured yet. The premise rests on published results that the team has
so far read only as abstracts — verifying them in full is a week 6 deliverable, not an
assumption.
