"""Every reference frame that does not involve a mutation, measured from one run dir.

A displacement means nothing without these. This script produces the rows of the
"reference-frame ladder": each is a wild-type-vs-wild-type comparison, so any motion it
reports is motion the predictor produced while the sequence stayed the same.

Frames, smallest expected to largest:

  esm_replicate     repeated ESMFold folds of one sequence. ESMFold has no seed and no
                    inference dropout, so this is expected to be exactly 0. Recorded,
                    not assumed.
  esm_recycle       consecutive recycle counts (--recycle-sweep). How much the structure
                    is still moving as it converges.
  af2_recycle       the same quantity for AF2 (--save-recycles). THIS is the floor that
                    compares across predictors: both have recycles, only AF2 has a seed.
  af2_seed          AF2 at fixed model, different --random-seed. AF2's native floor:
                    the seed sets the MSA cluster draw and masking.
  af2_model         AF2 at fixed seed, different model weights (1-5). NOT a floor -- a
                    separate, larger reference line, and the honest upper bound on
                    "AF2 moved and nothing about the sequence changed".
  af2_vs_esm        the two predictors on the identical WT sequence. Inter-predictor
                    disagreement, which bounds how much any cross-predictor claim can mean.

Floors are reported PER CONSTRUCT and never pooled: the seed draws from the MSA stack,
and the constructs have different alignment depths, so their floors are different numbers.

Usage:
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -I pipeline/noise_floor.py --run-dir RUN
Writes RUN/floors.csv and prints the table.
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from geometry import CAModel, load_ca, superpose  # noqa: E402

AF2_TAG = "alphafold2_ptm"


EXCLUDE_TAIL = 0        # set from --exclude-tail


def pair_stats(models: list[CAModel], pairs: list[tuple[int, int]]) -> dict | None:
    """Superpose the given index pairs and summarise the spread.

    With EXCLUDE_TAIL > 0 the last N residues are dropped before the RMSD is taken --
    superposition still uses the whole chain, only the scoring is restricted. Ubiquitin's
    C-terminal tail (72-76) swings ~1 A between recycles whatever the sequence is, which
    inflates a global-RMSD floor by roughly 6x and dilutes the mutation signal it is meant
    to be compared against. Dropping it measures the fold rather than the tether.
    """
    if not pairs:
        return None
    rmsds, maxima = [], []
    for i, j in pairs:
        d = superpose(models[i], models[j])
        if EXCLUDE_TAIL > 0 and len(d.resids):
            keep = d.resids < d.resids.max() - (EXCLUDE_TAIL - 1)
            if keep.sum() < 3:
                raise ValueError(f"--exclude-tail {EXCLUDE_TAIL} leaves <3 residues")
            per = d.per_residue[keep]
            rmsds.append(float(np.sqrt((per ** 2).mean())))
            maxima.append(float(per.max()))
            continue
        rmsds.append(d.rmsd)
        maxima.append(d.max_displacement)
    r = np.asarray(rmsds, dtype=float)
    return {
        "n_structures": len(models),
        "n_pairs": len(r),
        "rmsd_mean": float(r.mean()),
        "rmsd_p95": float(np.percentile(r, 95)),
        "rmsd_max": float(r.max()),
        "max_per_residue": float(np.max(maxima)),
    }


def all_pairs(models: list[CAModel]) -> dict | None:
    n = len(models)
    return pair_stats(models, [(i, j) for i in range(n) for j in range(i + 1, n)])


def consecutive_pairs(models: list[CAModel]) -> dict | None:
    """For an ordered series (recycles): only neighbouring steps.

    All-pairs would mix "recycle 0 vs recycle 3" (a convergence trajectory) into the
    same number as "recycle 2 vs recycle 3" (the residual wobble at convergence). Only
    the latter is a floor.

    `rmsd_mean` here still averages the whole trajectory, early steps included, so it
    OVERSTATES the floor. `rmsd_final_step` -- the last consecutive step -- is the number
    to quote. Both are reported so the convergence shape stays visible: if the final step
    is not the smallest, the prediction is still wobbling rather than converged, and the
    final step is an amplitude rather than a residual.
    """
    pairs = [(i, i + 1) for i in range(len(models) - 1)]
    stats = pair_stats(models, pairs)
    if stats is None:
        return None
    steps = [pair_stats(models, [pr])["rmsd_mean"] for pr in pairs]
    stats["rmsd_final_step"] = steps[-1]
    # the conservative floor in the converged regime: worst step in the second half
    stats["rmsd_plateau_max"] = max(steps[len(steps) // 2:])
    # strictly decreasing steps = converging; otherwise still wobbling, and the final
    # step is an oscillation amplitude rather than a converged residual
    stats["converging"] = all(b < a for a, b in zip(steps, steps[1:]))
    stats["steps"] = ";".join(f"{v:.4f}" for v in steps)
    return stats


def load_all(paths: list[Path]) -> list[CAModel]:
    return [load_ca(str(p)) for p in paths]


def esm_paths(esm: Path, protein: str, kind: str) -> list[Path]:
    """kind 'rep' -> replicate folds; kind 'rec' -> recycle sweep, ordered by count."""
    pat = re.compile(rf"^{re.escape(protein)}__WT__{kind}(\d+)\.pdb$")
    hits = [(int(m.group(1)), p) for p in sorted(esm.glob(f"{protein}__WT__{kind}*.pdb"))
            if (m := pat.match(p.name))]
    return [p for _, p in sorted(hits)]


def af2_final(cf: Path, protein: str, model: int | None = None,
              seed: int | None = None) -> list[tuple[int, int, Path]]:
    """Final (converged) AF2 WT structures as (model, seed, path), recycle files excluded."""
    pat = re.compile(rf"^{re.escape(protein)}__WT_unrelaxed_{AF2_TAG}"
                     rf"_model_(\d+)_seed_(\d+)\.pdb$")
    out = []
    for p in sorted(cf.glob(f"{protein}__WT_unrelaxed_{AF2_TAG}_model_*_seed_*.pdb")):
        m = pat.match(p.name)          # excludes the ".r<N>.pdb" recycle files
        if not m:
            continue
        mi, si = int(m.group(1)), int(m.group(2))
        if (model is None or mi == model) and (seed is None or si == seed):
            out.append((mi, si, p))
    return out


def af2_recycles(cf: Path, protein: str, model: int = 1, seed: int = 0) -> list[Path]:
    """--save-recycles output for one model/seed, ordered by recycle index."""
    stem = f"{protein}__WT_unrelaxed_{AF2_TAG}_model_{model}_seed_{seed:03d}"
    pat = re.compile(rf"^{re.escape(stem)}\.r(\d+)\.pdb$")
    hits = [(int(m.group(1)), p) for p in sorted(cf.glob(f"{stem}.r*.pdb"))
            if (m := pat.match(p.name))]
    return [p for _, p in sorted(hits)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--out", default=None, help="default: RUN/floors.csv")
    ap.add_argument("--exclude-tail", type=int, default=0,
                    help="drop the last N residues from the RMSD (superposition still uses "
                         "the whole chain). Use 5 for ubiquitin: its C-terminal tail swings "
                         "~1 A between recycles regardless of sequence.")
    a = ap.parse_args()

    global EXCLUDE_TAIL
    EXCLUDE_TAIL = a.exclude_tail
    if EXCLUDE_TAIL:
        print(f"scoring on all but the last {EXCLUDE_TAIL} residues\n")
    run = Path(a.run_dir)
    esm = run / "preds" / "esmfold"
    cf = run / "preds" / "colabfold"
    proteins = list(pd.read_csv(run / "constructs.csv")["protein"])

    rows: list[dict] = []

    def add(frame: str, protein: str, stats: dict | None, note: str = "") -> None:
        if stats is None:
            print(f"  [skip] {frame:14s} {protein:12s} -- not enough structures")
            return
        rows.append({"frame": frame, "protein": protein, **stats, "note": note})

    for p in proteins:
        print(f"--- {p} ---")

        # ESMFold replicate floor
        reps = esm_paths(esm, p, "rep")
        if len(reps) >= 2:
            add("esm_replicate", p, all_pairs(load_all(reps)),
                f"{len(reps)} reps, ESMFold has no seed")
        else:
            print(f"  [skip] esm_replicate  {p:12s} -- {len(reps)} replicate(s)")

        # ESMFold recycle-convergence floor
        recs = esm_paths(esm, p, "rec")
        if len(recs) >= 2:
            add("esm_recycle", p, consecutive_pairs(load_all(recs)),
                "consecutive recycle counts: " + ",".join(
                    re.search(r"__rec(\d+)\.pdb$", f.name).group(1) for f in recs))
        else:
            print(f"  [skip] esm_recycle    {p:12s} -- run --recycle-sweep first")

        # AF2 recycle-convergence floor (model_1 seed_000)
        ar = af2_recycles(cf, p)
        if len(ar) >= 2:
            add("af2_recycle", p, consecutive_pairs(load_all(ar)),
                "model_1 seed_000, consecutive recycles: " + ",".join(
                    re.search(r"\.r(\d+)\.pdb$", f.name).group(1) for f in ar))
        else:
            print(f"  [skip] af2_recycle    {p:12s} -- run --save-recycles first")

        finals = af2_final(cf, p)
        by_model: dict[int, list[Path]] = {}
        by_seed: dict[int, list[Path]] = {}
        for mi, si, path in finals:
            by_model.setdefault(mi, []).append(path)
            by_seed.setdefault(si, []).append(path)

        # AF2 seed floor: pool pairs WITHIN each model, never across models
        seed_rows = [all_pairs(load_all(paths)) for paths in
                     (v for v in by_model.values() if len(v) >= 2)]
        if seed_rows:
            tot = sum(s["n_pairs"] for s in seed_rows)
            add("af2_seed", p, {
                "n_structures": sum(s["n_structures"] for s in seed_rows),
                "n_pairs": tot,
                # pair-count-weighted, so a model with more seeds counts proportionally
                "rmsd_mean": sum(s["rmsd_mean"] * s["n_pairs"] for s in seed_rows) / tot,
                "rmsd_p95": max(s["rmsd_p95"] for s in seed_rows),
                "rmsd_max": max(s["rmsd_max"] for s in seed_rows),
                "max_per_residue": max(s["max_per_residue"] for s in seed_rows),
            }, f"within-model across seeds, {len(seed_rows)} model(s)")
        else:
            print(f"  [skip] af2_seed       {p:12s} -- need >=2 seeds at one model")

        # AF2 model spread: within each seed, across models. NOT a floor.
        model_rows = [all_pairs(load_all(paths)) for paths in
                      (v for v in by_seed.values() if len(v) >= 2)]
        if model_rows:
            tot = sum(s["n_pairs"] for s in model_rows)
            add("af2_model", p, {
                "n_structures": sum(s["n_structures"] for s in model_rows),
                "n_pairs": tot,
                "rmsd_mean": sum(s["rmsd_mean"] * s["n_pairs"] for s in model_rows) / tot,
                "rmsd_p95": max(s["rmsd_p95"] for s in model_rows),
                "rmsd_max": max(s["rmsd_max"] for s in model_rows),
                "max_per_residue": max(s["max_per_residue"] for s in model_rows),
            }, f"within-seed across models, {len(model_rows)} seed(s) -- model spread, not a floor")
        else:
            print(f"  [skip] af2_model      {p:12s} -- need >=2 models at one seed")

        # Inter-predictor, WT only
        ref = af2_final(cf, p, model=1, seed=0)
        if ref and reps:
            d = superpose(load_ca(str(reps[0])), load_ca(str(ref[0][2])))
            add("af2_vs_esm", p, {
                "n_structures": 2, "n_pairs": 1,
                "rmsd_mean": d.rmsd, "rmsd_p95": d.rmsd, "rmsd_max": d.rmsd,
                "max_per_residue": d.max_displacement,
            }, "ESMFold rep1 vs AF2 model_1 seed_000, same WT sequence")

    if not rows:
        print("\nNothing measured. Fold the WT replicates/seeds first.")
        return

    df = pd.DataFrame(rows)
    order = ["esm_replicate", "esm_recycle", "af2_recycle", "af2_seed",
             "af2_model", "af2_vs_esm"]
    df["_o"] = df.frame.map({f: i for i, f in enumerate(order)})
    df = df.sort_values(["_o", "protein"]).drop(columns="_o")

    out = Path(a.out) if a.out else run / "floors.csv"
    df.to_csv(out, index=False)

    print("\n" + "=" * 104)
    print("REFERENCE FRAMES WITHOUT A MUTATION  (every row is wild type vs wild type)")
    print("=" * 104)
    print(f"{'frame':14s} {'protein':12s} {'n':>3s} {'pairs':>5s} "
          f"{'mean RMSD':>10s} {'final step':>10s} {'p95':>8s} {'max':>8s} {'max/res':>8s}  note")
    print("-" * 118)
    has_final = "rmsd_final_step" in df.columns
    for _, r in df.iterrows():
        fs = r.rmsd_final_step if has_final else float("nan")
        fs_s = f"{fs:10.4f}" if fs == fs else f"{'--':>10s}"
        print(f"{r.frame:14s} {r.protein:12s} {r.n_structures:3d} {r.n_pairs:5d} "
              f"{r.rmsd_mean:10.4f} {fs_s} {r.rmsd_p95:8.4f} {r.rmsd_max:8.4f} "
              f"{r.max_per_residue:8.4f}  {r.note}")
    if has_final and df.rmsd_final_step.notna().any():
        print("\nFor the recycle frames quote 'final step', not 'mean RMSD': the mean includes")
        print("the early, far-from-converged steps and so overstates the floor.")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
