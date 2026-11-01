"""Superposition and per-residue displacement between a wild-type and a mutant model.

Global RMSD is reported because prior work reports it, not because it answers our
question. The quantity that does is `per_residue`: how far each residue moved. That is
what distinguishes a real, local response to a substitution from numerical noise spread
over the whole chain, and it is what the PyMOL renders are colored by.

Both ESMFold and AlphaFold write per-residue pLDDT into the B-factor column, so
`load_ca` returns it alongside the coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from Bio.PDB import PDBParser


@dataclass
class CAModel:
    """C-alpha trace of one chain, keyed by residue number."""

    resids: np.ndarray          # (N,) int residue numbers
    coords: np.ndarray          # (N, 3) float
    plddt: np.ndarray           # (N,) float, from the B-factor column
    source: str = ""

    def __len__(self) -> int:
        return len(self.resids)


def load_ca(pdb_path: str, chain: str | None = None) -> CAModel:
    """Read C-alpha coordinates and pLDDT from a predicted structure.

    `chain=None` takes the first chain, which is what single-chain predictions give.
    """
    structure = PDBParser(QUIET=True).get_structure("s", pdb_path)
    model = next(structure.get_models())
    ch = model[chain] if chain is not None else next(model.get_chains())

    resids, coords, plddt = [], [], []
    for res in ch:
        if "CA" not in res:
            continue                      # skip waters, ligands, incomplete residues
        ca = res["CA"]
        resids.append(res.id[1])
        coords.append(ca.get_coord())
        plddt.append(ca.get_bfactor())

    return CAModel(
        resids=np.asarray(resids, dtype=int),
        coords=np.asarray(coords, dtype=float),
        plddt=np.asarray(plddt, dtype=float),
        source=pdb_path,
    )


@dataclass
class Displacement:
    """Result of superposing a mutant on its wild type."""

    rmsd: float                 # global C-alpha RMSD after superposition, angstrom
    resids: np.ndarray          # residue numbers common to both models
    per_residue: np.ndarray     # (N,) displacement per residue, angstrom
    delta_plddt: np.ndarray     # (N,) mutant pLDDT minus wild-type pLDDT
    n_matched: int = 0
    meta: dict = field(default_factory=dict)

    @property
    def max_displacement(self) -> float:
        return float(self.per_residue.max()) if len(self.per_residue) else float("nan")

    @property
    def mean_delta_plddt(self) -> float:
        return float(self.delta_plddt.mean()) if len(self.delta_plddt) else float("nan")


def superpose(ref: CAModel, mut: CAModel) -> Displacement:
    """Superpose `mut` onto `ref` on shared residues and measure what moved.

    Residues are matched by number, not by order, so a model missing a terminal
    residue does not silently shift the whole comparison by one.
    """
    common = np.intersect1d(ref.resids, mut.resids)
    if len(common) < 3:
        raise ValueError(
            f"only {len(common)} residues in common between {ref.source!r} and "
            f"{mut.source!r}; cannot superpose"
        )

    ri = np.searchsorted(ref.resids, common)
    mi = np.searchsorted(mut.resids, common)
    ref_xyz, mut_xyz = ref.coords[ri], mut.coords[mi]

    rot, tran = _kabsch(mut_xyz, ref_xyz)
    moved = mut_xyz @ rot + tran

    diff = np.linalg.norm(moved - ref_xyz, axis=1)
    rmsd = float(np.sqrt((diff ** 2).mean()))

    return Displacement(
        rmsd=rmsd,
        resids=common,
        per_residue=diff,
        delta_plddt=mut.plddt[mi] - ref.plddt[ri],
        n_matched=len(common),
        meta={"ref": ref.source, "mut": mut.source,
              "n_ref": len(ref), "n_mut": len(mut)},
    )


def _kabsch(mobile: np.ndarray, target: np.ndarray):
    """Optimal rotation+translation taking `mobile` onto `target`."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    m, t = mobile - mc, target - tc
    u, _, vt = np.linalg.svd(m.T @ t)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    corr = np.diag([1.0, 1.0, d])
    rot = (u @ corr @ vt)
    return rot, tc - mc @ rot


def locality_profile(disp: Displacement, ref: CAModel, site: int,
                     shell: float = 8.0) -> dict:
    """How concentrated is the motion around the mutated residue?

    Returns mean displacement inside and outside a distance shell around `site`,
    and their ratio. A ratio near 1 means the motion is spread over the chain --
    i.e. noise. A ratio well above 1 means the model responded locally, which is
    what a physically meaningful response looks like.

    `site` is a residue number, in the same numbering as the structure.
    """
    if site not in set(disp.resids.tolist()):
        raise ValueError(f"residue {site} not among the matched residues")

    idx = int(np.searchsorted(ref.resids, site))
    dist = np.linalg.norm(ref.coords - ref.coords[idx], axis=1)
    keep = np.searchsorted(ref.resids, disp.resids)
    dist = dist[keep]

    near, far = dist <= shell, dist > shell
    mean_near = float(disp.per_residue[near].mean()) if near.any() else float("nan")
    mean_far = float(disp.per_residue[far].mean()) if far.any() else float("nan")

    return {
        "site": site,
        "shell_angstrom": shell,
        "n_near": int(near.sum()),
        "n_far": int(far.sum()),
        "mean_displacement_near": mean_near,
        "mean_displacement_far": mean_far,
        "locality_ratio": mean_near / mean_far if mean_far else float("nan"),
    }


def noise_floor(models: list[CAModel]) -> dict:
    """Spread across repeated folds of the SAME sequence.

    Run this before interpreting any mutant. Every pairwise RMSD among replicate
    wild-type predictions; the resulting distribution is the threshold a mutant
    effect has to clear to be real.
    """
    if len(models) < 2:
        raise ValueError("need at least 2 replicate models to estimate a noise floor")

    rmsds, maxima = [], []
    for i in range(len(models)):
        for j in range(i + 1, len(models)):
            d = superpose(models[i], models[j])
            rmsds.append(d.rmsd)
            maxima.append(d.max_displacement)

    rmsds = np.asarray(rmsds)
    return {
        "n_replicates": len(models),
        "n_pairs": len(rmsds),
        "rmsd_mean": float(rmsds.mean()),
        "rmsd_max": float(rmsds.max()),
        "rmsd_p95": float(np.percentile(rmsds, 95)),
        "max_per_residue_displacement": float(np.max(maxima)),
    }
