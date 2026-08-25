#!/usr/bin/env sage -python
"""Focused QH-SS-RLWE-as-LWE proxy without slow non-limiting attacks."""
import argparse
import importlib
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
estimator = importlib.import_module(".estimator", "lattice-estimator")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ring-n", type=int, required=True)
parser.add_argument("--qbits", type=int, required=True)
parser.add_argument("--sigma", type=float, required=True)
args = parser.parse_args()

params = estimator.lwe_parameters.LWEParameters(
    n=args.ring_n,
    q=2**args.qbits,
    Xs=estimator.nd.Ternary,
    Xe=estimator.nd.DiscreteGaussian(stddev=args.sigma),
    m=float("inf"),
    tag="algorithm3-focused-ring-proxy",
)
results = estimator.LWE.estimate(
    params,
    red_cost_model=estimator.RC.BDGL16,
    deny_list=("arora-gb", "bkw", "bdd_hybrid", "bdd_mitm_hybrid", "dual"),
    catch_exceptions=True,
)
security = min(math.log2(result["rop"]) for result in results.values())
print(f"focused ring proxy security={security:.1f} bits")
