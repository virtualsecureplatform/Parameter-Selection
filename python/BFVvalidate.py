#!/usr/bin/env python3
"""Validation checks for the BFV average-case estimator."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent
sys.path[:0] = [str(_REPO_ROOT), str(_SCRIPT_DIR)]

try:
    from python.noiseestimation.bfv import (
        add_many,
        fresh,
        mul,
        mul_dependent,
        noise_budget_log2,
    )
    from python.noiseestimation.params.bfv import openfhe_paper
except ModuleNotFoundError as exc:
    if exc.name != "python":
        raise
    from noiseestimation.bfv import (
        add_many,
        fresh,
        mul,
        mul_dependent,
        noise_budget_log2,
    )
    from noiseestimation.params.bfv import openfhe_paper


TABLE7_EXPECTED = {
    "enc": {
        "label": "Encryption",
        "q_bits": 60,
        "expected": {12: 32.0, 13: 31.5, 14: 31.0, 15: 30.5},
    },
    "add": {
        "label": "Addition",
        "q_bits": 60,
        "expected": {12: 31.5, 13: 31.0, 14: 30.5, 15: 30.0},
    },
    "mul": {
        "label": "Multiplication",
        "q_bits": 120,
        "expected": {12: 65.0, 13: 63.6, 14: 62.1, 15: 60.5},
    },
}


TABLE10_EXPECTED = {
    "add-same": {
        "label": "Same Addition",
        "q_bits": 60,
        "expected": {12: 31.0, 13: 30.5, 14: 30.0, 15: 29.5},
    },
    "mul-same": {
        "label": "Same Multiplication",
        "q_bits": 120,
        "expected": {12: 64.5, 13: 63.1, 14: 61.6, 15: 60.0},
    },
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Compare estimator output with 600.pdf OpenFHE validation values."
    )
    ap.add_argument(
        "--tolerance",
        type=float,
        default=0.06,
        help="allowed absolute error in bits after reproducing rounded table values",
    )
    return ap.parse_args()


def estimate_budget(op: str, nbit: int, q_bits: int) -> float:
    params = openfhe_paper(nbit=nbit, q_bits=q_bits)
    clean = fresh(params)
    if op == "enc":
        state = clean
    elif op == "add":
        state = clean.add(clean)
    elif op == "mul":
        state = mul(clean, clean, params)
    elif op == "add-same":
        state = add_many([clean, clean], dependent=True, label="same-add")
    elif op == "mul-same":
        state = mul_dependent(clean, clean, params)
    else:
        raise ValueError(f"unsupported operation: {op}")
    return noise_budget_log2(state.log2_variance, d=6.0)


def main() -> int:
    args = parse_args()
    print("OpenFHE BFV validation against 600.pdf Table 7")
    print("  parameters: t=65537, sigma=3.19, chi_s=chi_u=U3, Hybrid/HPSPOVERQ")
    print("  target column: paper average-case 'our' maximum-value budget")
    header = "operation       nbit  q_bits  estimated  paper  delta   status"
    print(header)
    print("-" * len(header))

    ok = True
    for op, spec in TABLE7_EXPECTED.items():
        for nbit, expected in spec["expected"].items():
            estimated = estimate_budget(op, nbit, spec["q_bits"])
            delta = estimated - expected
            passed = abs(delta) <= args.tolerance
            ok = ok and passed
            print(
                f"{spec['label']:<14} {nbit:4d} {spec['q_bits']:7d} "
                f"{estimated:10.2f} {expected:6.1f} {delta:7.2f} "
                f"{'OK' if passed else 'FAIL'}"
            )

    print("\nDependent-input checks from 600.pdf Table 10")
    for op, spec in TABLE10_EXPECTED.items():
        for nbit, expected in spec["expected"].items():
            estimated = estimate_budget(op, nbit, spec["q_bits"])
            delta = estimated - expected
            passed = abs(delta) <= args.tolerance
            ok = ok and passed
            print(
                f"{spec['label']:<19} {nbit:4d} {spec['q_bits']:7d} "
                f"{estimated:10.2f} {expected:6.1f} {delta:7.2f} "
                f"{'OK' if passed else 'FAIL'}"
            )

    if ok:
        print(f"\nPASS: all Table 7/10 checks within {args.tolerance:.2f} bits")
        return 0
    print(f"\nFAIL: at least one Table 7/10 check exceeded {args.tolerance:.2f} bits")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
