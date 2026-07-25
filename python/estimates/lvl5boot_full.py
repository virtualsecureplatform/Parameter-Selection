#!/usr/bin/env python3
"""Full lattice-estimator run for the final BFV bootstrap parameters:
n = 2^15, Q = 2^576, sparse ternary secret h = 64, sigma = 2^23."""
import sys, importlib

sys.path.insert(0, "/work")
est = importlib.import_module(".estimator", "lattice-estimator")
ND = est.nd
LWE = est.LWE
LWEParameters = est.lwe_parameters.LWEParameters

param = LWEParameters(
    n=32768,
    q=2**576,
    Xs=ND.SparseTernary(32, 32, 32768),
    Xe=ND.DiscreteGaussian(stddev=2**23),
    m=32768,
    tag="lvl5bootparam-final",
)
print(param, flush=True)
res = LWE.estimate(param, red_cost_model=est.RC.BDGL16)
print("=" * 70, flush=True)
for k, v in res.items():
    print(k, "::", v, flush=True)
