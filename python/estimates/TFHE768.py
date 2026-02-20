#!/bin/python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

param = estimator.lwe_parameters.LWEParameters(
    n=768,
    q= (2 ** 16 + 1)*(2 **8 +1),
    Xs=estimator.nd.UniformMod(3),
    # Xe=estimator.nd.DiscreteGaussian(stddev= 2 ** 7),
    Xe=estimator.nd.DiscreteGaussian(stddev= 8),
    tag="TFHE768",
)
print(param.n)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
