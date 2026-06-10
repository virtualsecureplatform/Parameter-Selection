#!/bin/python3
"""Check TFHEpp lvl6 CKKS security for the 896-bit bootstrap modulus."""

import importlib
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
try:
    estimator = importlib.import_module(".estimator", "lattice-estimator")
except ModuleNotFoundError as exc:
    if exc.name == "sage":
        raise SystemExit(
            "Run this check with Sage's Python, for example: "
            "sage -python estimates/CKKS_lvl6.py"
        ) from None
    raise

TARGET_BITS = 128.0
N = 32768
LOG_Q = 896
ALPHA_LOG2 = -886
NOISE_STDDEV_LOG2 = LOG_Q + ALPHA_LOG2


def rop_log2(result):
    return float(math.log(result["rop"], 2))


def check_param(label, xs):
    param = estimator.lwe_parameters.LWEParameters(
        n=N,
        q=2**LOG_Q,
        Xs=xs,
        Xe=estimator.nd.DiscreteGaussian(stddev=2**NOISE_STDDEV_LOG2),
        tag=f"TFHEpp CKKS lvl6 {label}",
    )
    attacks = [
        ("primal_usvp", estimator.LWE.primal_usvp),
        ("primal_bdd", estimator.LWE.primal_bdd),
        ("dual", estimator.LWE.dual),
        ("dual_hybrid", estimator.LWE.dual_hybrid),
    ]

    print(
        f"{label}: n={N} q=2^{LOG_Q} alpha=2^{ALPHA_LOG2} "
        f"sigma=2^{NOISE_STDDEV_LOG2}",
        flush=True,
    )
    weakest = None
    for name, attack in attacks:
        bits = rop_log2(attack(param, red_cost_model=estimator.RC.BDGL16))
        weakest = bits if weakest is None else min(weakest, bits)
        print(f"  {name}: {bits:.3f} bits", flush=True)

    print(f"  weakest: {weakest:.3f} bits", flush=True)
    return weakest


def main():
    failures = []
    for label, xs in [
        ("dense-ternary", estimator.nd.UniformMod(3)),
        ("sparse-H16", estimator.nd.SparseTernary(8, 8, n=N)),
    ]:
        weakest = check_param(label, xs)
        if weakest < TARGET_BITS:
            failures.append((label, weakest))
    if failures:
        details = ", ".join(
            f"{label}: {weakest:.3f} < {TARGET_BITS:.1f}"
            for label, weakest in failures
        )
        raise SystemExit(f"Below target: {details}")
    print("Passed", flush=True)


if __name__ == "__main__":
    main()
