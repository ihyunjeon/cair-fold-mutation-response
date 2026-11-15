"""ESMFold through the public ESM Atlas fold API -- DEV-SUBSET FALLBACK ONLY.

Used for subset_v0 because the laptop cannot hold ESMFold locally (fp32 = ~14 GB RAM;
fp16 on Apple MPS crashes on a mixed-dtype matmul). The full run uses run_esmfold.py on a
Hoffman2 GPU. Output format matches run_esmfold.py so summarize.py treats both alike;
`backend` in the JSON tells them apart. The API returns no pTM (recorded as NaN) and its
recycle count / precision are server-side and not reported.

Usage:  python run_esmfold_api.py --run-dir RUN_DIR [--wt-reps 3]
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"


def fold(seq: str, tries: int = 4) -> str:
    for i in range(tries):
        try:
            req = urllib.request.Request(URL, data=seq.encode(), method="POST",
                                         headers={"User-Agent": "cair-fold-mutation-response"})
            with urllib.request.urlopen(req, timeout=120) as r:
                txt = r.read().decode()
            if "ATOM" in txt:
                return txt
        except Exception as exc:  # transient server errors
            print(f"  retry {i + 1}: {exc}", flush=True)
        time.sleep(5 * (i + 1))
    raise RuntimeError("ESM Atlas API failed")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--wt-reps", type=int, default=3)
    args = ap.parse_args()
    run = Path(args.run_dir)
    outdir = run / "preds" / "esmfold"
    outdir.mkdir(parents=True, exist_ok=True)
    for _, s in pd.read_csv(run / "sequences.csv").iterrows():
        for k in range(1, (args.wt_reps if s.variant == "WT" else 1) + 1):
            stem = f"{s.seq_id}__rep{k}"
            if (outdir / f"{stem}.pdb").exists():
                continue
            t0 = time.time()
            pdb = fold(s.sequence)
            lines, ca_b = [], []
            for ln in pdb.splitlines():
                if ln.startswith("ATOM"):
                    b = float(ln[60:66])
                    b = b * 100.0 if b <= 1.0 else b    # API writes pLDDT as 0-1
                    ln = f"{ln[:60]}{b:6.2f}{ln[66:]}"
                    if ln[12:16].strip() == "CA":
                        ca_b.append(b)
                lines.append(ln)
            (outdir / f"{stem}.pdb").write_text("\n".join(lines) + "\n")
            meta = dict(seq_id=s.seq_id, rep=k, predictor="esmfold", model="esmfold_v1",
                        backend="esm_atlas_api", ptm=float("nan"),
                        mean_plddt=sum(ca_b) / len(ca_b), num_recycles=None, device="remote",
                        esm_dtype=None, runtime_s=round(time.time() - t0, 2))
            (outdir / f"{stem}.json").write_text(json.dumps(meta, indent=1))
            print(f"{stem}  pLDDT={meta['mean_plddt']:.1f}  {meta['runtime_s']}s", flush=True)
            time.sleep(1)  # be polite to a shared public server


if __name__ == "__main__":
    main()
