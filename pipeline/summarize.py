"""Collect every ESMFold and ColabFold prediction in a run directory into flat tables.

This only collects; it computes no comparisons. Outputs in <run>/tables/:

  structures.csv   one row per predicted structure (predictor, model, seed, replicate,
                   global pLDDT, pTM, pLDDT at the mutated site, path) joined with the
                   variant metadata (set, fitness/burial stratum, RSA, DMS score)
  per_residue.csv  one row per residue per structure: construct + UniProt numbering,
                   amino acid, pLDDT, C-alpha x/y/z (no PDB parsing needed downstream)
  pairs.csv        manifest of like-for-like comparisons to make:
                     kind=mut_vs_wt   mutant vs WT, same predictor/model/seed/rep
                     kind=wt_vs_wt    WT replicate vs WT replicate (noise floor)
                   Pairs are listed, not computed.

Usage:  python summarize.py --run-dir RUN_DIR
"""
from __future__ import annotations

import argparse
import itertools
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

CF_RE = re.compile(r"^(?P<seq_id>.+?)_unrelaxed_(?:rank_(?P<rank>\d+)_)?alphafold2_ptm_model_(?P<model>\d)_seed_(?P<seed>\d+)\.pdb$")


def read_ca(path: Path):
    s = PDBParser(QUIET=True).get_structure("x", path)
    rows = []
    for res in next(s.get_models()).get_list()[0]:
        if "CA" in res:
            ca = res["CA"]
            rows.append((res.id[1], seq1(res.get_resname()), ca.get_bfactor(), *ca.get_coord()))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    run = Path(ap.parse_args().run_dir)
    tdir = run / "tables"
    tdir.mkdir(exist_ok=True)

    seqs = pd.read_csv(run / "sequences.csv").set_index("seq_id")
    variants = pd.read_csv(run / "variants.csv")
    structs = []

    for js in sorted((run / "preds" / "esmfold").glob("*.json")):
        m = json.loads(js.read_text())
        structs.append(dict(seq_id=m["seq_id"], predictor="esmfold", model="esmfold_v1",
                            seed=0, replicate=m["rep"], rank=np.nan, ptm=m["ptm"],
                            path=str(js.with_suffix(".pdb").relative_to(run)),
                            runtime_s=m["runtime_s"], device=m["device"],
                            settings=f"backend={m.get('backend', 'hf_transformers')};lm_dtype={m.get('esm_dtype')};"
                                     f"recycles={m.get('num_recycles')}"))

    cf = run / "preds" / "colabfold"
    cfg = json.loads((cf / "config.json").read_text()) if (cf / "config.json").exists() else {}
    cf_settings = (f"colabfold={cfg.get('version', '?')};recycles={cfg.get('num_recycles')};"
                   f"max_msa={cfg.get('max_msa')};templates={cfg.get('use_templates')};msa=wt_reused")
    for pdb in sorted(cf.glob("*_unrelaxed_*.pdb")):
        g = CF_RE.match(pdb.name)
        if not g:
            continue
        score = cf / pdb.name.replace("_unrelaxed_", "_scores_").replace(".pdb", ".json")
        sj = json.loads(score.read_text()) if score.exists() else {}
        structs.append(dict(seq_id=g["seq_id"], predictor="alphafold2", model=f"model_{g['model']}_ptm",
                            seed=int(g["seed"]), replicate=1, rank=int(g["rank"]) if g["rank"] else np.nan, ptm=sj.get("ptm", np.nan),
                            max_pae=sj.get("max_pae", np.nan),
                            path=str(pdb.relative_to(run)), runtime_s=np.nan, device="cpu",
                            settings=cf_settings))

    st = pd.DataFrame(structs).reindex(columns=[
        "seq_id", "predictor", "model", "seed", "replicate", "rank", "ptm", "max_pae", "path",
        "runtime_s", "device", "settings"])
    st["structure_id"] = (st.seq_id + "|" + st.predictor + "|" + st.model + "|s" + st.seed.astype(str)
                          + "|r" + st.replicate.astype(str))
    st = st.join(seqs[["protein", "variant", "construct_start", "length", "seq_sha1", "sequence"]], on="seq_id")

    per_res = []
    for _, s in st.iterrows():
        rows = read_ca(run / s.path)
        pred_seq = "".join(r[1] for r in rows)
        if pred_seq != s.sequence:
            raise ValueError(f"{s.structure_id}: predicted sequence != input sequence")
        for resid, aa, b, x, y, z in rows:
            per_res.append((s.structure_id, resid, resid + s.construct_start - 1, aa, round(b, 2),
                            round(x, 3), round(y, 3), round(z, 3)))
    pr = pd.DataFrame(per_res, columns=["structure_id", "resid_construct", "resid_uniprot", "aa",
                                        "plddt", "ca_x", "ca_y", "ca_z"])
    # ESMFold PDBs round pLDDT to 1 unit; use the full-precision values from the JSON when present
    for _, s in st[st.predictor == "esmfold"].iterrows():
        meta = json.loads((run / s.path).with_suffix(".json").read_text())
        if "plddt_ca" in meta:
            m = pr.structure_id == s.structure_id
            pr.loc[m, "plddt"] = meta["plddt_ca"]

    g = pr.groupby("structure_id").plddt
    st["mean_plddt"] = st.structure_id.map(g.mean()).round(2)
    st["n_res"] = st.structure_id.map(g.size())
    st = st.merge(variants, on=["protein", "variant"], how="left")
    st["site_plddt"] = [
        pr.loc[(pr.structure_id == sid) & (pr.resid_uniprot == pos), "plddt"].squeeze() if pd.notna(pos) else np.nan
        for sid, pos in zip(st.structure_id, st.pos)
    ]
    st["site_plddt"] = pd.to_numeric(st.site_plddt, errors="coerce")

    # like-for-like pairing manifest
    key = ["protein", "predictor", "model", "seed", "replicate"]
    wt = st[st.variant == "WT"]
    mut = st[st.variant != "WT"]
    pairs = []
    for _, m in mut.iterrows():
        w = wt[(wt[key] == m[key]).all(axis=1)]
        if len(w) == 1:
            pairs.append(dict(kind="mut_vs_wt", protein=m.protein, predictor=m.predictor, model=m.model,
                              seed=m.seed, variant=m.variant, pos_uniprot=m.pos,
                              structure_a=w.structure_id.iloc[0], structure_b=m.structure_id))
    for (prot, pred), grp in wt.groupby(["protein", "predictor"]):
        for a, b in itertools.combinations(grp.itertuples(), 2):
            pairs.append(dict(kind="wt_vs_wt", protein=prot, predictor=pred,
                              model=a.model if a.model == b.model else f"{a.model}~{b.model}",
                              seed=a.seed if a.seed == b.seed else -1, variant="WT", pos_uniprot=np.nan,
                              structure_a=a.structure_id, structure_b=b.structure_id,
                              same_model=a.model == b.model))
    pairs = pd.DataFrame(pairs)

    cols = ["structure_id", "protein", "variant", "predictor", "model", "seed", "replicate", "rank",
            "set", "fitness", "burial", "rsa", "DMS_score", "DMS_score_bin", "pos", "mean_plddt",
            "site_plddt", "ptm", "max_pae", "n_res", "seq_sha1", "pdb_wt", "pdb_mut", "path",
            "runtime_s", "device", "settings"]
    st[cols].rename(columns={"pos": "pos_uniprot"}).to_csv(tdir / "structures.csv", index=False)
    pr.to_csv(tdir / "per_residue.csv", index=False)
    pairs.to_csv(tdir / "pairs.csv", index=False)
    print(st.groupby(["protein", "predictor", "variant"]).size().unstack(0).fillna(0).astype(int).to_string())
    print(f"structures={len(st)} per_residue_rows={len(pr)} pairs={len(pairs)} "
          f"({(pairs.kind == 'mut_vs_wt').sum()} mut_vs_wt, {(pairs.kind == 'wt_vs_wt').sum()} wt_vs_wt)")


if __name__ == "__main__":
    main()
