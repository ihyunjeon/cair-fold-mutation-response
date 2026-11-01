"""Variant parsing, mutant sequence construction, and the substitution assertion.

The assertion is the load-bearing part. DMS papers, UniProt and PDB files routinely
disagree about residue numbering -- usually by a constant offset from a cleaved signal
peptide or an expression tag. That mismatch never raises an error on its own. It
silently produces a complete set of plausible mutant sequences that answer a question
nobody asked, and you find out only when the results make no sense in week 9.

So: build mutants through `apply_variant`, and run `assert_variant_set` over the whole
panel before anything is folded.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

_VARIANT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

AA20 = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class Variant:
    """A single substitution, in 1-based numbering as written in the source table."""

    wt: str
    pos: int
    mut: str

    def __str__(self) -> str:
        return f"{self.wt}{self.pos}{self.mut}"


def parse_variant(token: str) -> Variant:
    """Parse 'A42G' into Variant('A', 42, 'G').

    Raises ValueError on anything that is not a single canonical substitution --
    including synonymous tokens like 'A42A', which are usually a sign that an
    unfiltered DMS table leaked through.
    """
    token = token.strip().upper()
    m = _VARIANT_RE.match(token)
    if not m:
        raise ValueError(f"not a single substitution token: {token!r}")
    wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
    if wt not in AA20 or mut not in AA20:
        raise ValueError(f"non-canonical amino acid in {token!r}")
    if wt == mut:
        raise ValueError(f"synonymous token {token!r}; filter these out upstream")
    if pos < 1:
        raise ValueError(f"position must be 1-based and positive: {token!r}")
    return Variant(wt, pos, mut)


def apply_variant(seq: str, variant: str | Variant, offset: int = 0) -> str:
    """Return `seq` with the substitution applied, asserting the wild-type residue first.

    Parameters
    ----------
    seq
        The wild-type sequence this project folds, as a plain string.
    variant
        Either 'A42G' or a Variant.
    offset
        Added to the variant position before indexing. Use this when the DMS table
        numbers residues differently from `seq` -- for example `offset=-22` when the
        table numbers the precursor and `seq` is the mature chain. Determine it once
        per protein, write it down in data/README.md, and never guess it per variant.

    Raises
    ------
    ValueError
        If the position falls outside the sequence, or if the residue found there is
        not the wild-type residue the variant claims. Do not catch this and continue.
    """
    v = parse_variant(variant) if isinstance(variant, str) else variant
    idx = v.pos + offset - 1  # to 0-based

    if not (0 <= idx < len(seq)):
        raise ValueError(
            f"{v}: position {v.pos} with offset {offset} -> index {idx}, "
            f"outside a sequence of length {len(seq)}"
        )

    found = seq[idx]
    if found != v.wt:
        raise ValueError(
            f"{v}: expected wild-type {v.wt} at position {v.pos} "
            f"(offset {offset}, index {idx}) but found {found}. "
            f"This is a numbering mismatch, not a bad variant -- fix the offset."
        )

    return seq[:idx] + v.mut + seq[idx + 1 :]


def assert_variant_set(seq: str, variants, offset: int = 0) -> dict:
    """Check a whole panel before folding anything. Returns a summary dict.

    Reports every failure rather than stopping at the first, because a systematic
    offset error fails on nearly all variants at once and you want to see that
    pattern rather than one example of it.
    """
    failures = []
    for token in variants:
        try:
            apply_variant(seq, token, offset=offset)
        except ValueError as exc:
            failures.append((str(token), str(exc)))

    total = len(list(variants))
    summary = {
        "total": total,
        "passed": total - len(failures),
        "failed": len(failures),
        "failures": failures,
    }

    if failures and len(failures) == total:
        summary["diagnosis"] = (
            "every variant failed -- this is almost certainly a constant numbering "
            "offset, not bad data. Try inferring it before changing anything else."
        )
    return summary


def infer_offset(seq: str, variants, search=range(-50, 51)) -> list[int]:
    """Return offsets under which every variant's wild-type residue matches.

    A diagnostic, not a fix. If this returns exactly one offset, record it in
    data/README.md with the reason (signal peptide, tag, construct numbering) and pass
    it explicitly from then on. If it returns several, the panel is too small to
    disambiguate -- do not pick one arbitrarily.
    """
    ok = []
    for off in search:
        if all(_matches(seq, t, off) for t in variants):
            ok.append(off)
    return ok


def _matches(seq: str, token, offset: int) -> bool:
    try:
        v = parse_variant(token) if isinstance(token, str) else token
    except ValueError:
        return False
    idx = v.pos + offset - 1
    return 0 <= idx < len(seq) and seq[idx] == v.wt
