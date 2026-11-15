"""Build the mutant sequence set for a fold run (subset or full panel).

For each protein:
  1. take the folded construct from the ProteinGym target sequence,
  2. compute per-residue burial (relative SASA) on the experimental WT structure,
  3. stratify DMS single mutants on fitness x burial and sample per stratum,
  4. add the experimentally crystallized mutants (for the 6-way comparison),
  5. run the substitution assertion (src/mutations.py) over everything,
  6. write variants.csv and one FASTA per sequence.

Usage:
  python build_subset.py --out OUT_DIR --per-stratum 2 --seed 0
  (per-stratum 25 gives the ~100-per-protein full panel)
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from mutations import apply_variant, assert_variant_set  # noqa: E402

RAW = Path("/Volumes/LaCie/data/cair_fold/raw")

# Tien et al. 2013 theoretical max ASA (A^2)
MAX_ASA = dict(A=129, R=274, N=195, D=193, C=167, E=223, Q=225, G=104, H=224, I=197,
               L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)
BURIED_RSA = 0.20

# construct = target_seq[start-1:end], numbered start..end (UniProt / DMS numbering).
# offset passed to apply_variant() against the construct string = 1 - start.
PROTEINS = {
    "ubq_yeast": dict(dms="RL40A_YEAST_Roscoe_2013", start=1, end=76,
                      burial_pdb="1UBQ", burial_chain="A",
                      note="yeast ubiquitin 1-76 (DMS organism); burial from human 1UBQ, "
                           "differs at P19S, E24D, A28S (all surface)"),
    "pdz3_rat": dict(dms="DLG4_RAT_McLaughlin_2012", start=302, end=402,
                     burial_pdb="5HDY", burial_chain="A",
                     note="rat PSD95 PDZ3 302-402 = native part of the 5HDY construct "
                          "(297-301 GSPEF and 403-415 are expression tags)"),
}

# Constructs that exist only for the experimental comparison (no DMS of this exact seq).
HUMAN_UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
EXPERIMENTAL = [
    # protein, variant, pdb_wt, pdb_mut, chain, note
    ("ubq_human", "L67S", "1UBQ", "5OXI", "A", "C-terminally retracted mutant; real backbone change expected"),
    ("ubq_human", "K33S", "1UBQ", "4PIH", "A", "surface Lys; near-zero change expected (control)"),
    ("ubq_human", "K11S", "1UBQ", "4PIG", "A", "surface Lys; near-zero change expected (control)"),
    ("pdz3_rat", "G330T", "5HDY", "5HET", "A", "binding-groove beta2-beta3 loop; same lab/construct as WT"),
    ("pdz3_rat", "H372A", "5HDY", "5HF4", "A", "alpha2 His372 (ligand pocket); same lab/construct as WT"),
]


def rsa_table(pdb_id: str, chain: str) -> pd.DataFrame:
    s = PDBParser(QUIET=True).get_structure(pdb_id, RAW / "pdb" / f"{pdb_id}.pdb")
    model = s[0]
    # strip everything but the protein chain of interest
    for ch in list(model):
        if ch.id != chain:
            model.detach_child(ch.id)
    ch = model[chain]
    for res in list(ch):
        if res.id[0] != " ":
            ch.detach_child(res.id)
    ShrakeRupley().compute(model, level="R")
    rows = []
    for res in ch:
        aa = seq1(res.get_resname())
        rows.append(dict(resid=res.id[1], pdb_aa=aa, sasa=res.sasa,
                         rsa=min(res.sasa / MAX_ASA[aa], 1.5) if aa in MAX_ASA else np.nan))
    return pd.DataFrame(rows)


def seq_hash(s: str) -> str:
    return hashlib.sha1(s.encode()).hexdigest()[:10]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-stratum", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    out = Path(args.out)
    (out / "fasta").mkdir(parents=True, exist_ok=True)

    ref = pd.read_csv(RAW / "proteingym" / "DMS_substitutions.csv").set_index("DMS_id")
    rng = np.random.default_rng(args.seed)
    rows, constructs = [], {}

    for prot, cfg in PROTEINS.items():
        target = ref.loc[cfg["dms"], "target_seq"]
        construct = target[cfg["start"] - 1: cfg["end"]]
        offset = 1 - cfg["start"]
        constructs[prot] = dict(seq=construct, start=cfg["start"], note=cfg["note"])

        dms = pd.read_csv(RAW / "proteingym" / "DMS_ProteinGym_substitutions" / f"{cfg['dms']}.csv")
        dms = dms[~dms.mutant.str.contains(":")].copy()
        dms["pos"] = dms.mutant.str[1:-1].astype(int)
        dms = dms[(dms.pos >= cfg["start"]) & (dms.pos <= cfg["end"])].copy()

        rsa = rsa_table(cfg["burial_pdb"], cfg["burial_chain"]).set_index("resid")
        dms["rsa"] = dms.pos.map(rsa.rsa)
        dms = dms.dropna(subset=["rsa"])
        dms["burial"] = np.where(dms.rsa < BURIED_RSA, "buried", "exposed")

        # fitness strata from the DMS score distribution of THIS assay
        lo, hi = dms.DMS_score.quantile([0.25, 0.5]).values
        dms["fitness"] = np.select([dms.DMS_score <= lo, dms.DMS_score >= hi],
                                   ["deleterious", "neutral"], default="intermediate")

        for (fit, bur), grp in dms[dms.fitness != "intermediate"].groupby(["fitness", "burial"]):
            # one variant per position so a stratum is not two hits at one site
            grp = grp.sample(frac=1, random_state=int(rng.integers(1 << 31))).drop_duplicates("pos")
            for _, r in grp.head(args.per_stratum).iterrows():
                rows.append(dict(protein=prot, variant=r.mutant, pos=r.pos, set="dms_stratified",
                                 fitness=fit, burial=bur, rsa=round(r.rsa, 3),
                                 DMS_score=r.DMS_score, DMS_score_bin=r.DMS_score_bin,
                                 dms_id=cfg["dms"], pdb_wt=cfg["burial_pdb"], pdb_mut=""))
        dms_scores = dms.set_index("mutant")

        for p, v, pw, pm, ch, note in EXPERIMENTAL:
            if p != prot:
                continue
            pos = int(v[1:-1])
            r = dms_scores.loc[v] if v in dms_scores.index else None
            rows.append(dict(protein=prot, variant=v, pos=pos, set="experimental_pair",
                             fitness=(r.fitness if r is not None else ""),
                             burial=("buried" if rsa.rsa.get(pos, 1) < BURIED_RSA else "exposed"),
                             rsa=round(float(rsa.rsa.get(pos, np.nan)), 3),
                             DMS_score=(r.DMS_score if r is not None else np.nan),
                             DMS_score_bin=(r.DMS_score_bin if r is not None else np.nan),
                             dms_id=cfg["dms"], pdb_wt=pw, pdb_mut=pm, note=note))

    # human ubiquitin construct for the experimental pairs (crystals are human)
    constructs["ubq_human"] = dict(seq=HUMAN_UBQ, start=1, note="human ubiquitin 1-76, matches 1UBQ/5OXI/4PIH/4PIG")
    rsa_h = rsa_table("1UBQ", "A").set_index("resid")
    for p, v, pw, pm, ch, note in EXPERIMENTAL:
        if p == "ubq_human":
            pos = int(v[1:-1])
            rows.append(dict(protein=p, variant=v, pos=pos, set="experimental_pair", fitness="",
                             burial=("buried" if rsa_h.rsa[pos] < BURIED_RSA else "exposed"),
                             rsa=round(float(rsa_h.rsa[pos]), 3), DMS_score=np.nan,
                             DMS_score_bin=np.nan, dms_id="", pdb_wt=pw, pdb_mut=pm, note=note))

    df = pd.DataFrame(rows).drop_duplicates(["protein", "variant"], keep="last")

    # ---- the gate: every variant must carry the intended WT residue at the intended spot
    for prot, c in constructs.items():
        toks = df[df.protein == prot].variant.tolist()
        summ = assert_variant_set(c["seq"], toks, offset=1 - c["start"])
        print(f"[assert] {prot}: {summ['passed']}/{summ['total']} passed")
        if summ["failed"]:
            raise SystemExit(f"substitution assertion failed for {prot}: {summ['failures']}")

    seq_rows = []
    for prot, c in constructs.items():
        seq_rows.append(dict(seq_id=f"{prot}__WT", protein=prot, variant="WT",
                             sequence=c["seq"], construct_start=c["start"]))
    for _, r in df.iterrows():
        c = constructs[r.protein]
        mseq = apply_variant(c["seq"], r.variant, offset=1 - c["start"])
        seq_rows.append(dict(seq_id=f"{r.protein}__{r.variant}", protein=r.protein, variant=r.variant,
                             sequence=mseq, construct_start=c["start"]))
    seqs = pd.DataFrame(seq_rows)
    seqs["seq_sha1"] = seqs.sequence.map(seq_hash)
    seqs["length"] = seqs.sequence.str.len()

    df.to_csv(out / "variants.csv", index=False)
    seqs.to_csv(out / "sequences.csv", index=False)
    pd.DataFrame([dict(protein=k, construct_start=v["start"], length=len(v["seq"]),
                       offset_vs_dms=1 - v["start"], note=v["note"], sequence=v["seq"])
                  for k, v in constructs.items()]).to_csv(out / "constructs.csv", index=False)
    for _, s in seqs.iterrows():
        (out / "fasta" / f"{s.seq_id}.fasta").write_text(f">{s.seq_id}\n{s.sequence}\n")
    print(df.groupby(["protein", "set", "fitness", "burial"]).size().to_string())
    print(f"{len(seqs)} sequences -> {out}")


if __name__ == "__main__":
    main()
