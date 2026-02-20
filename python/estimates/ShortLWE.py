#!/bin/python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

param = estimator.lwe_parameters.LWEParameters(
    n=128*3,
    q=2 ** 8,
    Xs=estimator.nd.UniformMod(2),
   #  Xe=estimator.nd.DiscreteGaussian(stddev=2 ** (16-10.5)),
    Xe=estimator.nd.CenteredBinomial(24),
    # m = 384,
    tag="ShortLWE",
)
print(param.n)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
