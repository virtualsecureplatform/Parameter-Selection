#!/bin/python3
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")
# param = estimator.schemes.TFHE630
param = estimator.lwe_parameters.LWEParameters(
    n=586,
    q=2 ** 32,
    Xs=estimator.nd.NoiseDistribution.UniformMod(2),
    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=0.000_092_511_997_467_675_6 * 2 ** 32),
    tag="TFHE586",
)
r = estimator.LWE.estimate(param)