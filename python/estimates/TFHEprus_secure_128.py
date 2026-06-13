#!/bin/python3
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib

estimator = importlib.import_module(".estimator", "lattice-estimator")

GOLDILOCKS_Q = 0xffff_ffff_0000_0001
TFHEPRUS_STDDEV = 2 ** 14

params = [
    estimator.lwe_parameters.LWEParameters(
        n=2048,
        q=GOLDILOCKS_Q,
        Xs=estimator.nd.UniformMod(2),
        Xe=estimator.nd.DiscreteGaussian(stddev=TFHEPRUS_STDDEV),
        tag="TFHEprus-secure-128-input-binary",
    ),
    estimator.lwe_parameters.LWEParameters(
        n=2048,
        q=GOLDILOCKS_Q,
        Xs=estimator.nd.UniformMod(3),
        Xe=estimator.nd.DiscreteGaussian(stddev=TFHEPRUS_STDDEV),
        tag="TFHEprus-secure-128-glwe-ternary-estimate",
    ),
]

for param in params:
    print("=" * 60)
    print(f"Estimating: {param.tag} (n={param.n})")
    print("=" * 60)
    estimator.LWE.estimate(param, red_cost_model=estimator.RC.BDGL16)
    print()
