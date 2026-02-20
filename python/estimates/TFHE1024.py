#!/bin/python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

param = estimator.lwe_parameters.LWEParameters(
    n=1024,
    # q= 5**4*2**16+1,
   #  q= 3**4*2**16+1,
    q=2 ** 32,
   #  q=2 ** 25,
    Xs=estimator.nd.UniformMod(3),
    # Xe=estimator.nd.DiscreteGaussian(stddev= 2 ** 7),
    # Xe=estimator.nd.DiscreteGaussian(stddev= 4.2),
   #  Xe=estimator.nd.CenteredBinomial(2),
    # Xe=estimator.nd.DiscreteGaussian(stddev= 2*4/4),
    # Xe=estimator.nd.DiscreteGaussian(stddev= 0.000_000_034_233_878_701_836_9 * 2 ** 32),
    Xe=estimator.nd.DiscreteGaussian(stddev= 2 ** (32-25)),
    tag="TFHE1024",
)
print(param.n)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
