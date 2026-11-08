# Protein structure rendering

**Your role:** produce the project's primary result figures: PyMOL renders of wild type
and mutant structures, overlaid and colored by how far each residue moved. In this
project the renders are the result, not decoration added at the end. A good render lets an
audience *see* whether a predictor responded to a mutation, where, and how that compares
to what the real protein does.

Everything you need is in `inputs/`. Background is in the repo-root `README.md`; read the REFERENCE
FRAMES section, which explains why every figure has an experimental panel.

---

## 1. What you have

`inputs/sixway/` holds five **6-way comparisons**. Each is one point mutation that has been
crystallized, shown six ways: {experimental, AlphaFold2, ESMFold} × {wild type, mutant}.

| folder | protein | mutation | crystals (WT vs mutant) | what to expect |
|---|---|---|---|---|
| `ubq_human_L67S` | human ubiquitin | L67S (buried core Leu) | 1UBQ vs 5OXI | the real protein changes: "C-terminally retracted" |
| `ubq_human_K33S` | human ubiquitin | K33S (surface Lys) | 1UBQ vs 4PIH | control: little real change |
| `ubq_human_K11S` | human ubiquitin | K11S (surface Lys) | 1UBQ vs 4PIG | control: little real change |
| `pdz3_rat_G330T` | rat PSD95 PDZ3 | G330T (binding-groove loop) | 5HDY vs 5HET | same lab, same construct |
| `pdz3_rat_H372A` | rat PSD95 PDZ3 | H372A (α2 His, ligand pocket) | 5HDY vs 5HF4 | same lab, same construct |

"What to expect" is from the crystal paper titles and residue positions. **It is not a
measured result**, so don't caption a figure with it until the numbers back it up.

Each folder contains:

```
exp_wt.pdb    exp_mut.pdb     crystal structures (chain A, altloc A, waters/ligands removed)
af_wt.pdb     af_mut.pdb      AlphaFold2 predictions
esm_wt.pdb    esm_mut.pdb     ESMFold predictions
load_sixway.pml               loads + aligns all six, sets colors and representation
manifest.csv                  where each file came from + the residue found at the mutated site
```

All six files share **one residue numbering** (UniProt: ubiquitin 1-76, PDZ3 302-402), so
`resi 67` means the same residue in every object. That's already been checked:
`manifest.csv:site_residue` shows the WT residue in every `*_wt` file and the mutant
residue in every `*_mut` file.

**B-factor column:** crystals hold crystallographic B-factors. Predictions hold **pLDDT**
(model confidence, 0-100). These are different quantities; never put them on one color scale.

**This is development data** (`subset_v0`), folded on a laptop with reduced settings. Use
it to build and tune the rendering pipeline. Final figures will be re-rendered from the
Hoffman2 run (`full_v1`). That run delivers the same folder layout, plus more mutants for
the "strongest / weakest responder" panels.

---

## 2. Setup

```bash
conda create -n cair_pymol -c conda-forge python=3.11 pymol-open-source pandas numpy
conda activate cair_pymol
cd yash/protein_structure_rendering/inputs/sixway/ubq_human_L67S
pymol load_sixway.pml          # must be run from inside the folder (paths are relative)
```

Tested with PyMOL open-source 3.1.0. All five scenes load with no errors.

What `load_sixway.pml` does:
- loads the six objects and superposes the other five onto `exp_wt` using **every C-alpha with
  `cycles=0`** (no outlier rejection). PyMOL's default `align` throws away the residues that
  moved most, which hides exactly the motion we're trying to show. Keep `cycles=0`.
- cartoon for all, sticks for the mutated side chain (`select site`).
- colors: experimental grey (WT `grey70`, mutant `grey30`), AF2 blue (`lightblue` / `tv_blue`),
  ESMFold orange (`lightorange` / `orange`). WT light, mutant dark.
- grid slots: 1 = experimental, 2 = AF2, 3 = ESMFold, each with WT and mutant overlaid.
  `set grid_mode, 1` shows the three panels side by side.
- **It does not ray-trace or save images.** That's your job.

---

## 3. Your task, in order

**Step 1: a consistent house style.** Write one shared PyMOL script (e.g. `style.pml` or a
Python function) that every figure uses: white background, cartoon settings, the color
scheme above, lighting, and the site stick representation. Then every figure in the talk
looks like it came from the same lab.
Useful settings: `bg_color white`, `set ray_opaque_background, 1`, `set cartoon_fancy_helices, 1`,
`set ray_trace_mode, 1`, `set antialias, 2`, `set orthoscopic, 1`, `set ray_shadows, 0`.

**Step 2: a fixed camera per protein.** Orient once per protein so the mutated site is
visible and not hidden behind the chain. Save the view with `get_view` and reuse it
with `set_view` for every panel and every mutant of that protein. The same protein must
always appear from the same angle, or the audience can't compare panels.

**Step 3: overlay renders (no data needed).** For each of the five folders, render the
3-panel figure (experimental | AF2 | ESMFold, each WT+mutant overlaid) with `grid_mode 1`,
then `ray` + `png`. Make a batch script that loops over all `inputs/sixway/*/` folders
and writes to `outputs/`, so the full run re-renders with one command.
Size: slides are 16:9. `ray 2400,1350` then `png <file>.png, dpi=300` works well.

**Step 4: displacement coloring (the primary figure).** Color each mutant by how far each
residue moved from its own wild type. The analysis member delivers this as
`outputs/displacement_<protein>_<variant>.csv`:

```
object,resid_uniprot,displacement_A
exp_mut,67,<float>
af_mut,67,<float>
esm_mut,67,<float>
```

Load it into the B-factor column and color by it:

```python
import pandas as pd
d = pd.read_csv("displacement_ubq_human_L67S.csv")
for obj, g in d.groupby("object"):
    vals = dict(zip(g.resid_uniprot, g.displacement_A))
    cmd.alter(f"{obj}", "b = vals.get(int(resi), -1)", space={"vals": vals})
    cmd.spectrum("b", "white_red", f"{obj} and b > -1", minimum=FLOOR, maximum=VMAX)
```

Rules for the color scale:
- **the same `minimum`/`maximum` on every panel of a figure** (experimental, AF2, ESMFold).
  Otherwise a 0.2 Å wiggle and a 2 Å shift look equally red.
- `minimum` = the **noise floor** (from the analysis member's `outputs/noise_floor.csv`).
  Anything at or below what the predictor does when nothing changes stays white.
- residues missing from a structure (e.g. 5OXI and 4PIG chain A end at residue 75) get
  `b = -1`. Color them grey, not white, so "missing" never reads as "didn't move".
- include a color bar with units (Å) in the figure or on the slide.

*Until the CSVs exist*, build and test this code on the pLDDT already in the B-factor
column of `af_mut`/`esm_mut`, then swap in the real file. Don't compute displacement
yourself in PyMOL for the figures. One source of numbers for the whole project.

**Step 5: the responder panel (from full_v1).** For each protein, the strongest and
weakest responders (picked by the analysis member from the full run), rendered with the
same house style, camera and color scale. These are the slides the talk is built around.

---

## 4. Deliverables and timeline

| week | deliverable |
|---|---|
| 8 | house style + fixed cameras + batch overlay renders of the 5 dev pairs |
| 9 | displacement-colored 6-way figures; responder panel from full_v1 |
| 10 | final figures at slide resolution, color bars, consistent captions |

Save scripts in `protein_structure_rendering/` (commit them). Save images in
`protein_structure_rendering/outputs/` (PNG; git-ignored unless the team decides otherwise). Write down any choice that changes how a
figure reads (camera, color range, which residues are shown) in your notebook
`yash/yash.ipynb`.

## 5. Things that make a render misleading (avoid all of them)

- default `align` with outlier rejection (hides the motion). Always `cycles=0`.
- a different color range per panel.
- coloring pLDDT and calling it displacement, or coloring crystal B-factors next to pLDDT.
- a camera that hides the mutated site, or changes between WT and mutant panels.
- captions that claim a result the analysis hasn't measured yet.
