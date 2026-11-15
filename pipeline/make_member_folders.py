"""Copy a fold run into the two member folders under notebooks/.

  notebooks/julie/disagreement_data_analysis/inputs/   tables + every predicted PDB + trimmed crystals
  notebooks/yash/protein_structure_rendering/inputs/  the 6-way folders (PDBs + load_sixway.pml)

The TASK.md files are written by hand and are not touched. Rerun this with
--run-dir pointing at full_v1 once the Hoffman2 run is back; layouts are identical.

Usage:  python make_member_folders.py --run-dir /Volumes/LaCie/data/cair_fold/subset_v0
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    run = Path(ap.parse_args().run_dir)

    # ---- member 1: disagreement analysis
    m1 = REPO / "notebooks" / "julie" / "disagreement_data_analysis" / "inputs"
    if m1.exists():
        shutil.rmtree(m1)
    (m1 / "tables").mkdir(parents=True)
    for f in ["structures.csv", "per_residue.csv", "pairs.csv"]:
        shutil.copy(run / "tables" / f, m1 / "tables" / f)
    for f in ["variants.csv", "sequences.csv", "constructs.csv"]:
        shutil.copy(run / f, m1 / "tables" / f)
    st = pd.read_csv(run / "tables" / "structures.csv")
    for p in st.path:                       # keep paths identical to structures.csv:path
        dst = m1 / p
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(run / p, dst)
        js = (run / p).with_suffix(".json")
        if js.exists():
            shutil.copy(js, dst.with_suffix(".json"))
    # crystals, trimmed to the folded construct and on UniProt numbering (made by prep_sixway)
    exp = m1 / "experimental"
    exp.mkdir()
    rows = []
    for d in sorted((run / "sixway").iterdir()):
        man = pd.read_csv(d / "manifest.csv")
        for tag in ["exp_wt", "exp_mut"]:
            src = man.loc[man.file == f"{tag}.pdb", "source"].iloc[0]   # "RCSB 1UBQ chain A ..."
            pdb_id = src.split()[1]
            if not (exp / f"{pdb_id}.pdb").exists():
                shutil.copy(d / f"{tag}.pdb", exp / f"{pdb_id}.pdb")
                rows.append(dict(pdb_id=pdb_id, protein=d.name.rsplit("_", 1)[0],
                                 role="WT" if tag == "exp_wt" else d.name.rsplit("_", 1)[1],
                                 file=f"experimental/{pdb_id}.pdb", source=src))
    pd.DataFrame(rows).to_csv(exp / "experimental_structures.csv", index=False)

    # ---- member 2: rendering
    m2 = REPO / "notebooks" / "yash" / "protein_structure_rendering" / "inputs"
    if m2.exists():
        shutil.rmtree(m2)
    shutil.copytree(run / "sixway", m2 / "sixway")

    for root in [m1, m2]:
        n = sum(1 for _ in root.rglob("*") if _.is_file())
        size = sum(f.stat().st_size for f in root.rglob("*") if f.is_file()) / 1e6
        print(f"{root.relative_to(REPO)}: {n} files, {size:.1f} MB")


if __name__ == "__main__":
    main()
