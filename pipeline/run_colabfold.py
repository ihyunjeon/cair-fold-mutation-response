"""AlphaFold2 (ColabFold) arm.

MSA policy (decided 2026-11-08): **reuse the WT MSA**. One MMseqs2 server query per
construct; each mutant gets the WT alignment with only the query row substituted. This
is standard practice and is the condition the AF2-insensitivity papers critique.

Model policy:
  WT sequences  -> --wt-models models (default all 5) x --wt-seeds seeds (noise floor + model spread)
  mutants       -> --mut-models models (default 1 = model_1), seed 0
Comparisons in the analysis must pair like with like (mutant model_k seed_s vs WT model_k seed_s).

No templates, no Amber relax, 3 recycles.

Step 1  python run_colabfold.py --run-dir R --step msa      (needs internet, 1 query/construct)
Step 2  python run_colabfold.py --run-dir R --step fold     (offline; laptop CPU or GPU node)
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd

COMMON = ["--num-recycle", "3", "--model-type", "alphafold2_ptm", "--random-seed", "0"]


def sh(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def msa_step(run: Path, colabfold_batch: str) -> None:
    seqs = pd.read_csv(run / "sequences.csv")
    wt = seqs[seqs.variant == "WT"]
    qdir = run / "colabfold" / "msa_query"
    qdir.mkdir(parents=True, exist_ok=True)
    csv = qdir / "wt_queries.csv"
    wt[["seq_id", "sequence"]].rename(columns={"seq_id": "id"}).to_csv(csv, index=False)
    sh([colabfold_batch, "--msa-only", str(csv), str(run / "colabfold" / "msa_wt")])


def build_mutant_a3m(run: Path) -> Path:
    """Copy each WT a3m with the query row (and its header) replaced by the mutant."""
    seqs = pd.read_csv(run / "sequences.csv")
    src = run / "colabfold" / "msa_wt"
    dst = run / "colabfold" / "a3m"
    dst.mkdir(parents=True, exist_ok=True)
    for _, s in seqs.iterrows():
        wt_id = f"{s.protein}__WT"
        lines = (src / f"{wt_id}.a3m").read_text().splitlines()
        # first '>' line is the query header, the next line its sequence
        qi = next(i for i, ln in enumerate(lines) if ln.startswith(">"))
        wt_seq = lines[qi + 1]
        if len(wt_seq) != len(s.sequence):
            raise ValueError(f"{s.seq_id}: query length {len(wt_seq)} != {len(s.sequence)}")
        ndiff = sum(a != b for a, b in zip(wt_seq, s.sequence))
        if ndiff != (0 if s.variant == "WT" else 1):
            raise ValueError(f"{s.seq_id}: {ndiff} differences between WT query row and mutant")
        lines[qi] = f">{s.seq_id}"
        lines[qi + 1] = s.sequence
        (dst / f"{s.seq_id}.a3m").write_text("\n".join(lines) + "\n")
    return dst


def fold_step(run: Path, colabfold_batch: str, data: str, wt_seeds: int, mut_models: int,
              wt_models: int = 5, max_msa: str | None = None, shard: str = "0/1",
              groups: str = "wt,mut") -> None:
    a3m = build_mutant_a3m(run)
    seqs = pd.read_csv(run / "sequences.csv")
    k, n = map(int, shard.split("/"))          # run every n-th sequence starting at k
    seqs = seqs.iloc[k::n]
    for group, sel, extra in [
        ("wt", seqs.variant == "WT", ["--num-models", str(wt_models), "--num-seeds", str(wt_seeds)]),
        ("mut", seqs.variant != "WT", ["--num-models", str(mut_models), "--num-seeds", "1"]),
    ]:
        if not sel.any() or group not in groups.split(","):
            continue
        gdir = run / "colabfold" / f"in_{group}_shard{k}of{n}"
        gdir.mkdir(parents=True, exist_ok=True)
        for sid in seqs[sel].seq_id:
            link = gdir / f"{sid}.a3m"
            if not link.exists():
                link.symlink_to(a3m / f"{sid}.a3m")
        if max_msa:  # DEV ONLY: shallower MSA stacks for CPU speed; never for the real run
            extra = [*extra, "--max-msa", max_msa]
        sh([colabfold_batch, *COMMON, *extra, "--data", data,
            str(gdir), str(run / "preds" / "colabfold")])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--step", choices=["msa", "fold"], required=True)
    ap.add_argument("--colabfold-batch", default="colabfold_batch")
    ap.add_argument("--data", default="/Volumes/LaCie/data/cair_fold/weights/colabfold")
    ap.add_argument("--wt-seeds", type=int, default=2)
    ap.add_argument("--mut-models", type=int, default=1)
    ap.add_argument("--wt-models", type=int, default=5)
    ap.add_argument("--max-msa", default=None, help="e.g. 512:1024 -- dev subset on CPU only")
    ap.add_argument("--groups", default="wt,mut", help="which groups to fold: wt, mut, or both")
    ap.add_argument("--shard", default="0/1", help="k/n: fold every n-th sequence from k (parallel CPU jobs)")
    a = ap.parse_args()
    run = Path(a.run_dir)
    if a.step == "msa":
        msa_step(run, a.colabfold_batch)
    else:
        fold_step(run, a.colabfold_batch, a.data, a.wt_seeds, a.mut_models, a.wt_models, a.max_msa, a.shard, a.groups)


if __name__ == "__main__":
    main()
