#!/bin/python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

param = estimator.lwe_parameters.LWEParameters(
   n=4096,
   q=2 ** 128,
   Xs=estimator.nd.UniformMod(3),
   Xe=estimator.nd.DiscreteGaussian(stddev=2**(128-105)),
   tag="TFHE4096",
)
print(param.n)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
