#!/bin/python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

param = estimator.lwe_parameters.LWEParameters(
    n=760,
    q=2 ** 32,
    Xs=estimator.nd.Binary,
    # Xe=estimator.nd.DiscreteGaussian(stddev= 2 ** 7),
    Xe=estimator.nd.DiscreteGaussian(stddev= 2 ** (32-17)),
    tag="TFHE760",
)
print(param.n)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
