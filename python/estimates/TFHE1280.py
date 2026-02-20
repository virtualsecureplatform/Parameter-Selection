#!/bin/python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

param = estimator.lwe_parameters.LWEParameters(
   n=256*5,
   q=2 ** 64,
   Xs=estimator.nd.UniformMod(3),
#    Xe=estimator.nd.DiscreteGaussian(stddev=0.0000000000034525330484572114*2**64),
   Xe=estimator.nd.DiscreteGaussian(stddev=2**-30*2**(64)),
   tag="TFHE1280"
)
print(param.n)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
