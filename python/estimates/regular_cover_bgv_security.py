#!/usr/bin/env sage -python
"""LWE-proxy estimate for the ordinary Binary-NTT RLWE source of regular-cover BGV."""

from __future__ import annotations

import argparse
import importlib
import math
from pathlib import Path
import sys


_SCRIPT = Path(__file__).resolve()
_PYTHON = _SCRIPT.parents[1]
sys.path.insert(0, str(_PYTHON))
estimator = importlib.import_module(".estimator", "lattice-estimator")


def minimum_rop_bits(costs) -> tuple[str, float]:
    best_name = ""
    best = math.inf
    for name, cost in costs.items():
        rop = cost.get("rop", math.inf)
        try:
            bits = float(rop.log(2))
        except AttributeError:
            bits = math.log2(float(rop))
        if bits < best:
            best_name, best = name, bits
    return best_name, best


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--qbits", type=int, required=True)
    parser.add_argument("--sigma", type=float, required=True)
    parser.add_argument("--samples", type=int, default=None)
    parser.add_argument("--target-bits", type=float, default=128.0)
    parser.add_argument("--mode", choices=("rough", "full"), default="rough")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sample_count = args.samples if args.samples is not None else math.inf
    parameters = estimator.LWE.Parameters(
        n=args.degree,
        q=2**args.qbits,
        Xs=estimator.ND.Binary,
        Xe=estimator.ND.DiscreteGaussian(stddev=args.sigma),
        m=sample_count,
        tag="regular-cover-binary-ntt-source-proxy",
    )
    if args.mode == "rough":
        costs = estimator.LWE.estimate.rough(parameters)
    else:
        costs = estimator.LWE.estimate(
            parameters, red_cost_model=estimator.RC.BDGL16, quiet=True
        )
    attack, bits = minimum_rop_bits(costs)
    status = "PASS_ESTIMATE" if bits >= args.target_bits else "FAIL_ESTIMATE"
    print("Regular-cover ordinary-source estimate")
    print(
        f"  n={args.degree} log2(q)={args.qbits} sigma={args.sigma} "
        f"m={sample_count}"
    )
    print(f"  best attack={attack} log2(rop)={bits:.2f} {status}")
    print("  caveat: this flattens Binary-NTT RLWE to an LWE attack-cost proxy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
