"""Stage the 6-way comparison (experimental / AF2 / ESMFold  x  WT / mutant) for PyMOL.

For every experimental_pair in variants.csv, writes <run>/sixway/<protein>_<variant>/:

  exp_wt.pdb   exp_mut.pdb    crystal chain, protein atoms only, altloc A (or blank)
                              only, trimmed to the folded construct; B = crystallographic B
  af_wt.pdb    af_mut.pdb     AF2 model_1_ptm seed 0 (the like-for-like pair)
  esm_wt.pdb   esm_mut.pdb    ESMFold rep1
                              predictions renumbered to the crystal/UniProt numbering,
                              chain A, B = pLDDT (0-100)
  load_sixway.pml             loads + superposes all six onto exp_wt on C-alpha,
                              sets colors/representation. Does NOT ray-trace or save PNGs.
  manifest.csv                provenance of every file

Every staged structure is checked to carry the expected residue at the mutated site.
Usage:  python prep_sixway.py --run-dir RUN_DIR
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.SeqUtils import seq1

RAW = Path("/Volumes/LaCie/data/cair_fold/raw/pdb")


class CrystalSelect(Select):
    def __init__(self, chain, lo, hi):
        self.chain, self.lo, self.hi = chain, lo, hi

    def accept_chain(self, ch):
        return ch.id == self.chain

    def accept_residue(self, r):
        return r.id[0] == " " and self.lo <= r.id[1] <= self.hi

    def accept_atom(self, a):
        return a.get_altloc() in (" ", "A")


def write_crystal(pdb_id, chain, lo, hi, out):
    s = PDBParser(QUIET=True).get_structure(pdb_id, RAW / f"{pdb_id}.pdb")
    for a in s.get_atoms():                      # drop altloc labels after selecting A
        if a.get_altloc() == "A":
            a.set_altloc(" ")
    io = PDBIO(); io.set_structure(s); io.save(str(out), CrystalSelect(chain, lo, hi))


def write_prediction(src, offset, out):
    s = PDBParser(QUIET=True).get_structure("p", src)
    ch = next(s[0].get_chains())
    for r in sorted(ch, key=lambda r: -r.id[1]):  # renumber high->low to avoid id clashes
        r.id = (" ", r.id[1] + offset, " ")
    ch.id = "A"
    io = PDBIO(); io.set_structure(s); io.save(str(out))


def residue_at(path, resid):
    s = PDBParser(QUIET=True).get_structure("x", path)
    for r in next(s[0].get_chains()):
        if r.id[1] == resid:
            return seq1(r.get_resname())
    return None


PML = """\
# 6-way comparison: {title}
# exp = crystal (grey/black), af = AlphaFold2 (blue), esm = ESMFold (orange). WT light, mutant dark.
# Superposition: every object onto exp_wt, C-alpha only, cycles=0 (no outlier rejection,
# so real motion is not hidden). RMSDs print to the PyMOL console.
# Nothing here ray-traces or writes images.
# Run from inside this folder:  cd <this folder> && pymol load_sixway.pml
reinitialize
load exp_wt.pdb
load exp_mut.pdb
load af_wt.pdb
load af_mut.pdb
load esm_wt.pdb
load esm_mut.pdb
remove hydro
{aligns}
hide everything
show cartoon
set cartoon_transparency, 0
color grey70, exp_wt
color grey30, exp_mut
color lightblue, af_wt
color tv_blue, af_mut
color lightorange, esm_wt
color orange, esm_mut
select site, resi {site}
show sticks, site and not name N+C+O
util.cnc site
set grid_mode, 1
set grid_slot, 1, exp_wt
set grid_slot, 1, exp_mut
set grid_slot, 2, af_wt
set grid_slot, 2, af_mut
set grid_slot, 3, esm_wt
set grid_slot, 3, esm_mut
set grid_mode, 0
orient exp_wt
deselect
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    run = Path(ap.parse_args().run_dir)
    st = pd.read_csv(run / "tables" / "structures.csv")
    var = pd.read_csv(run / "variants.csv")
    cons = pd.read_csv(run / "constructs.csv").set_index("protein")

    for _, v in var[var.set == "experimental_pair"].iterrows():
        start = int(cons.loc[v.protein, "construct_start"])
        end = start + int(cons.loc[v.protein, "length"]) - 1
        d = run / "sixway" / f"{v.protein}_{v.variant}"
        d.mkdir(parents=True, exist_ok=True)
        man = []

        for tag, pdb_id in [("exp_wt", v.pdb_wt), ("exp_mut", v.pdb_mut)]:
            write_crystal(pdb_id, "A", start, end, d / f"{tag}.pdb")
            man.append(dict(file=f"{tag}.pdb", source=f"RCSB {pdb_id} chain A altloc A, resi {start}-{end}",
                            bfactor="crystallographic B"))

        def pick(variant, predictor):
            sel = st[(st.protein == v.protein) & (st.variant == variant) & (st.predictor == predictor)]
            sel = sel[(sel.model == "model_1_ptm") & (sel.seed == 0)] if predictor == "alphafold2" else sel[sel.replicate == 1]
            if len(sel) != 1:
                raise ValueError(f"{v.protein} {variant} {predictor}: expected 1 structure, got {len(sel)}")
            return sel.iloc[0]

        try:
            picks = {tag: pick(variant, predictor) for tag, variant, predictor in
                     [("af_wt", "WT", "alphafold2"), ("af_mut", v.variant, "alphafold2"),
                      ("esm_wt", "WT", "esmfold"), ("esm_mut", v.variant, "esmfold")]}
        except ValueError as exc:
            print(f"SKIP {d.name}: {exc}")
            continue
        for tag, s in picks.items():
            write_prediction(run / s.path, start - 1, d / f"{tag}.pdb")
            man.append(dict(file=f"{tag}.pdb", source=s.structure_id, bfactor="pLDDT 0-100",
                            mean_plddt=s.mean_plddt))

        # the gate: WT residue in every *_wt file, mutant residue in every *_mut file
        for m in man:
            want = v.variant[-1] if "_mut" in m["file"] else v.variant[0]
            got = residue_at(d / m["file"], int(v.pos))
            m["site_residue"] = got
            if got != want:
                raise ValueError(f"{d.name}/{m['file']}: residue {v.pos} is {got}, expected {want}")

        aligns = "\n".join(f"align {o} and name CA, exp_wt and name CA, cycles=0, object=aln_{o}"
                           for o in ["exp_mut", "af_wt", "af_mut", "esm_wt", "esm_mut"])
        (d / "load_sixway.pml").write_text(PML.format(title=f"{v.protein} {v.variant} ({v.pdb_wt} vs {v.pdb_mut})",
                                                      here=d, aligns=aligns + "\ndelete aln_*", site=int(v.pos)))
        pd.DataFrame(man).to_csv(d / "manifest.csv", index=False)
        print(f"staged {d}")


if __name__ == "__main__":
    main()
