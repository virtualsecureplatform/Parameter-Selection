#!/bin/python3
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")
# param = estimator.schemes.TFHE630
# param = estimator.lwe_parameters.LWEParameters(
#     n=586,
#     q=2 ** 32,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=0.000_092_511_997_467_675_6 * 2 ** 32),
#     tag="TFHE586",
# )
param = estimator.lwe_parameters.LWEParameters(
   n=636,
   q=2 ** 32,
   Xs=estimator.nd.NoiseDistribution.UniformMod(2),
   Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=0.000_092_511_997_467_675_6 * 2 ** 32),
   tag="TFHE636",
)
# param = estimator.lwe_parameters.LWEParameters(
#     n=1024,
#     q=2 ** 32,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(3),
#     # Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 2 ** 7),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 0.000_000_034_233_878_701_836_9 * 2 ** 32),
#     tag="TFHE1024",
# )
r = estimator.LWE.estimate(param)
