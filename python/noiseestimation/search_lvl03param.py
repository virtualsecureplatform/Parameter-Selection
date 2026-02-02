#!/usr/bin/env python3
"""
Parameter sweep helper for lvl03 (blind rotation to lvl3).

This script focuses on the constraints discussed in the repo:
- Keep `nbit`, `k`, and fresh noise `α` fixed (lvl3 target).
- Use DD/bivariate limb cover for TRGSW rows: `l̅ * B̅gbit = 128`.
- Enforce FFT/double mantissa safety (with margin bits):
    Bgbit + B̅gbit + nbit + margin < 53

The noise model is `python/noiseestimation/keyvariation.py` (variance), and we
rank candidates by the user cost proxy:
    (k*la + l) + (k*l̅a + l̅)
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

# Allow running as a script (without installing as a package).
_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT))

from python.noiseestimation.keyvariation import brnoisecalc  # noqa: E402
from python.noiseestimation.params.λ128bit import lvl0param  # noqa: E402


def _log2_std_over_q(var: float, q: int) -> float:
    return math.log2(math.sqrt(var) / q)


def _ok_fft(Bgbit: int, Bbarbit: int, nbit: int, margin: int) -> bool:
    return Bgbit + Bbarbit + nbit + margin < 53


def _eval_lvl03(
    *,
    nbit: int,
    k: int,
    alpha_pow2: float,
    l: int,
    la: int,
    Bgbit: int,
    Bga_bit: int,
    lbar: int,
    Bbarbit: int,
) -> tuple[float, int]:
    q = 2**128
    alpha = q * alpha_pow2
    sigma = alpha**2

    attrs = {
        "nbit": nbit,
        "k": k,
        "n": 1 << nbit,
        "q": q,
        "l": l,
        "la": la,
        "Bbit": Bgbit,
        "Babit": Bga_bit,
        "B": 2**Bgbit,
        "Ba": 2**Bga_bit,
        # DD/bivariate limb parameters (TRGSW rows)
        "l̅": lbar,
        "l̅a": lbar,
        "B̅gbit": Bbarbit,
        "B̅gabit": Bbarbit,
        # Fresh noise
        "α": alpha,
        "σ": sigma,
        # Secret distribution (ternary / {-1,0,+1})
        "variance_key_coefficient": 2.0 / 3.0,
        "expectation_key_coefficient": 0.0,
    }

    P = type("P", (), attrs)
    lvl03 = type("lvl03", (), {"domainP": lvl0param, "targetP": P})
    var = brnoisecalc(lvl03)
    log2std = _log2_std_over_q(var, q)
    cost = (k * la + l) + (k * lbar + lbar)
    return log2std, cost


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-log2-std-over-q", type=float, default=-64.0)
    ap.add_argument("--margin-bits", type=int, default=3)
    ap.add_argument("--lbar", type=int, default=8)
    ap.add_argument("--bbarbit", type=int, default=16)
    ap.add_argument("--l-max", type=int, default=6)
    ap.add_argument("--bgbit-min", type=int, default=8)
    ap.add_argument("--bgbit-max", type=int, default=21)
    args = ap.parse_args()

    nbit = 12
    k = 1
    alpha_pow2 = 2.0**-105  # TFHEpp lvl3param α (normalized)

    if args.lbar * args.bbarbit != 128:
        raise SystemExit("require full limb cover: l̅ * B̅gbit must equal 128")

    rows: list[tuple[int, float, int, int, int, int]] = []
    for l in range(1, args.l_max + 1):
        for la in range(1, args.l_max + 1):
            for Bgbit in range(args.bgbit_min, args.bgbit_max + 1):
                for Bga_bit in range(args.bgbit_min, args.bgbit_max + 1):
                    if not _ok_fft(Bgbit, args.bbarbit, nbit, args.margin_bits):
                        continue
                    if not _ok_fft(Bga_bit, args.bbarbit, nbit, args.margin_bits):
                        continue

                    log2std, cost = _eval_lvl03(
                        nbit=nbit,
                        k=k,
                        alpha_pow2=alpha_pow2,
                        l=l,
                        la=la,
                        Bgbit=Bgbit,
                        Bga_bit=Bga_bit,
                        lbar=args.lbar,
                        Bbarbit=args.bbarbit,
                    )
                    if log2std <= args.min_log2_std_over_q:
                        rows.append((cost, log2std, l, la, Bgbit, Bga_bit))

    rows.sort()
    print(f"candidates: {len(rows)}")
    for cost, log2std, l, la, Bgbit, Bga_bit in rows[:30]:
        print(
            f"cost={cost:2d}  log2(std/q)={log2std:8.3f}  l={l} la={la}  Bgbit={Bgbit} Bgₐbit={Bga_bit}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
