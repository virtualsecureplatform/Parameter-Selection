#!/usr/bin/env sage -python
"""Security proxy for the no-boundary 128-bit DD/QH Algorithm-3 candidate."""

import importlib
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
estimator = importlib.import_module(".estimator", "lattice-estimator")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sigma", type=float, default=2**23)
parser.add_argument("--qbits", type=int, default=128)
parser.add_argument("--dimension", type=int, default=4096)
args = parser.parse_args()

params = estimator.lwe_parameters.LWEParameters(
    n=args.dimension,
    q=2**args.qbits,
    Xs=estimator.nd.Ternary,
    Xe=estimator.nd.DiscreteGaussian(stddev=args.sigma),
    m=float("inf"),
    tag="algorithm3-dd-qh-ss-rlwe-proxy",
)
results = estimator.LWE.estimate(
    params, red_cost_model=estimator.RC.BDGL16, catch_exceptions=True
)
security = min(math.log2(result["rop"]) for result in results.values())
print(f"DD/QH RLWE-as-LWE proxy security={security:.1f} bits")
if security < 128:
    raise SystemExit("DD/QH RLWE proxy does not meet 128 bits")
