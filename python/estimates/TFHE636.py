#!/bin/python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

param = estimator.lwe_parameters.LWEParameters(
   n=636,
   q=2 ** 16,
   Xs=estimator.nd.Binary,
   Xe=estimator.nd.DiscreteGaussian(stddev=0.000_092_511_997_467_675_6 * 2 ** 16),
   tag="TFHE636",
)
print(param.n)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
