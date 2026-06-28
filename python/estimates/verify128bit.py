#!/bin/python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

params = [
    estimator.lwe_parameters.LWEParameters(
        n=636,
        q=2 ** 16,
        Xs=estimator.nd.Binary,
        Xe=estimator.nd.DiscreteGaussian(stddev=0.000_092_511_997_467_675_6 * 2 ** 16),
        tag="TFHE636",
    ),
    estimator.lwe_parameters.LWEParameters(
        n=1024,
        q=2 ** 32,
        Xs=estimator.nd.UniformMod(3),
        Xe=estimator.nd.DiscreteGaussian(stddev=2 ** (32-25)),
        tag="TFHE1024",
    ),
    estimator.lwe_parameters.LWEParameters(
        n=2048,
        q=2 ** 64,
        Xs=estimator.nd.UniformMod(3),
        Xe=estimator.nd.DiscreteGaussian(stddev=2**(64-51)),
        tag="TFHE2048",
    ),
    estimator.lwe_parameters.LWEParameters(
        n=4096,
        q=2 ** 128,
        Xs=estimator.nd.UniformMod(3),
        Xe=estimator.nd.DiscreteGaussian(stddev=2**(128-105)),
        tag="TFHE4096",
    ),
]

for param in params:
    print("=" * 60)
    print(f"Estimating: {param.tag} (n={param.n})")
    print("=" * 60)
    r = estimator.LWE.estimate(param, red_cost_model=estimator.RC.BDGL16)
    print()
