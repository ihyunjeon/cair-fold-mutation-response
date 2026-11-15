"""Fold every sequence in sequences.csv with ESMFold (HF transformers port of esmfold_v1).

WT sequences are folded --wt-reps times to measure the run-to-run noise floor. ESMFold
has no dropout or sampling at inference, so replicates may come out identical; that is
recorded, not assumed.

Writes  <out>/preds/esmfold/<seq_id>__rep<k>.pdb  (pLDDT 0-100 in B-factor)
        <out>/preds/esmfold/<seq_id>__rep<k>.json (ptm, mean pLDDT, runtime, settings)

Usage:  HF_HOME=... python run_esmfold.py --run-dir RUN_DIR [--wt-reps 3] [--device cpu]
Works on CPU (laptop) and CUDA (Hoffman2).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

NUM_RECYCLES = 4  # esmfold_v1 default


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--wt-reps", type=int, default=3)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--lm-dtype", choices=["fp16", "bf16", "fp32"], default=None,
                    help="ESM-2 language-model precision (default fp16 on cuda, fp32 elsewhere)")
    args = ap.parse_args()

    torch.set_num_threads(args.threads)
    run = Path(args.run_dir)
    outdir = run / "preds" / "esmfold"
    outdir.mkdir(parents=True, exist_ok=True)
    seqs = pd.read_csv(run / "sequences.csv")

    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    # ESM-2 3B language model in fp16 (as in the ESMFold reference code); folding trunk and
    # heads in fp32. Loading in fp16 first keeps peak RAM ~7 GB instead of ~14 GB.
    esm_dtype = args.lm_dtype or ("fp16" if args.device == "cuda" else "fp32")
    dt = dict(fp16=torch.float16, bf16=torch.bfloat16, fp32=torch.float32)[esm_dtype]
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", torch_dtype=dt)
    # promote everything outside the language model (incl. the bare esm_s_combine parameter)
    for name, t in list(model.named_parameters()) + list(model.named_buffers()):
        if not name.startswith("esm.") and t.is_floating_point():
            t.data = t.data.float()
    model = model.to(args.device).eval()

    for _, s in seqs.iterrows():
        reps = args.wt_reps if s.variant == "WT" else 1
        for k in range(1, reps + 1):
            stem = f"{s.seq_id}__rep{k}"
            if (outdir / f"{stem}.pdb").exists():
                continue
            t0 = time.time()
            ids = tok([s.sequence], return_tensors="pt", add_special_tokens=False)["input_ids"].to(args.device)
            with torch.no_grad():
                out = model(ids, num_recycles=NUM_RECYCLES)
            pdb = model.output_to_pdb(out)[0]
            # HF writes pLDDT in 0-1 into B-factor; rescale to 0-100 like ESM/AF2 do
            plddt_ca = out["plddt"][0, :, 1].float().cpu()  # atom37 index 1 = CA
            scale = 100.0 if float(plddt_ca.max()) <= 1.0 else 1.0
            lines = []
            for ln in pdb.splitlines():
                if ln.startswith("ATOM"):
                    b = float(ln[60:66]) * scale
                    ln = f"{ln[:60]}{b:6.2f}{ln[66:]}"
                lines.append(ln)
            (outdir / f"{stem}.pdb").write_text("\n".join(lines) + "\n")
            meta = dict(seq_id=s.seq_id, rep=k, predictor="esmfold", model="esmfold_v1",
                        ptm=float(out["ptm"]), mean_plddt=float(plddt_ca.mean()) * scale,
                        num_recycles=NUM_RECYCLES, device=args.device, esm_dtype=esm_dtype,
                        runtime_s=round(time.time() - t0, 2), torch=torch.__version__,
                        # full precision; the PDB B-factor column is rounded to 1 pLDDT unit
                        plddt_ca=[round(float(x) * scale, 3) for x in plddt_ca])
            (outdir / f"{stem}.json").write_text(json.dumps(meta, indent=1))
            print(f"{stem}  pLDDT={meta['mean_plddt']:.1f}  pTM={meta['ptm']:.3f}  {meta['runtime_s']}s", flush=True)


if __name__ == "__main__":
    main()
