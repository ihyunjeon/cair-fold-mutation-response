"""Per-residue displacement tables for the 3D renders, one CSV per six-way case.

The rendering arm colours each structure by how far every residue moved from its OWN
wild type -- crystal mutant against crystal WT, AF2 mutant against AF2 WT, ESMFold
mutant against ESMFold WT. Computing that in PyMOL would give the project two sources
for the same number, so it is computed here, once, and shipped as data.

Writes  <out>/displacement_<protein>_<variant>.csv   object,resid_uniprot,displacement_A
        <out>/noise_floor.csv                        protein,predictor,floor_A,...

Residue numbering is UniProt throughout: the six-way PDBs are already renumbered onto
the crystal's numbering, so predictions and crystals share a scale. A residue present in
the mutant but absent from its WT (5OXI and 4PIG chain A stop at 75) is simply not
emitted -- the renderer paints anything missing grey, which is the honest colour for it.

Usage:  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -I pipeline/make_render_inputs.py \
            --sixway RUN/sixway --out notebooks/yash/protein_structure_rendering/inputs/displacement
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from geometry import load_ca, superpose  # noqa: E402

# tag -> (wild-type file, mutant file). The object names are what the PyMOL session
# calls each structure; they are part of the contract with the rendering arm.
OBJECTS = {"exp_mut": ("exp_wt.pdb", "exp_mut.pdb"),
           "af_mut": ("af_wt.pdb", "af_mut.pdb"),
           "esm_mut": ("esm_wt.pdb", "esm_mut.pdb")}


def case_rows(case_dir: Path) -> tuple[list[tuple[str, int, float]], list[str]]:
    rows: list[tuple[str, int, float]] = []
    notes: list[str] = []
    for obj, (wt_f, mut_f) in OBJECTS.items():
        wt_p, mut_p = case_dir / wt_f, case_dir / mut_f
        if not (wt_p.exists() and mut_p.exists()):
            notes.append(f"{obj}: missing structure, skipped")
            continue
        wt, mut = load_ca(str(wt_p)), load_ca(str(mut_p))
        d = superpose(wt, mut)
        if d.n_matched < len(mut):
            notes.append(f"{obj}: {d.n_matched}/{len(mut)} residues matched")
        rows += [(obj, int(r), round(float(v), 4))
                 for r, v in zip(d.resids, d.per_residue)]
    return rows, notes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sixway", required=True, help="run_dir/sixway")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    six, out = Path(args.sixway), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for case_dir in sorted(p for p in six.iterdir() if p.is_dir()):
        rows, notes = case_rows(case_dir)
        if not rows:
            print(f"{case_dir.name}: nothing to write")
            continue
        dest = out / f"displacement_{case_dir.name}.csv"
        with dest.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["object", "resid_uniprot", "displacement_A"])
            w.writerows(rows)
        per_obj = {o: sum(1 for r in rows if r[0] == o) for o, _ in OBJECTS.items()}
        print(f"{dest.name}: {len(rows)} rows  " +
              "  ".join(f"{k}={v}" for k, v in per_obj.items() if v) +
              ("  [" + "; ".join(notes) + "]" if notes else ""))


if __name__ == "__main__":
    main()
