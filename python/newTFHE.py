#!/bin/python3
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")
param = estimator.lwe_parameters.LWEParameters(
    n=511,
    q=2 ** 16,
    Xs=estimator.nd.NoiseDistribution.UniformMod(2),
    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=2 ** (12-10)),
    # Xe=estimator.nd.NoiseDistribution.UniformMod(8),
    tag="ShortLWE",
)
# param = estimator.schemes.LightSaber
# param = estimator.lwe_parameters.LWEParameters(
#     n=586,
#     q=2 ** 32,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=0.000_092_511_997_467_675_6 * 2 ** 32),
#     tag="TFHE586",
# )
# param = estimator.lwe_parameters.LWEParameters(
#    n=636,
#    q=2 ** 32,
#    Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=0.000_092_511_997_467_675_6 * 2 ** 32),
#    tag="TFHE636",
# )
# param = estimator.lwe_parameters.LWEParameters(
#     n=777,
#     q=2 ** 32,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=2 ** (32-17.6)),
#     tag="TFHE586",
# )
# param = estimator.lwe_parameters.LWEParameters(
#     n=1024,
#     q=2 ** 32,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(3),
#     # Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 2 ** 7),
#     # Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 0.000_000_034_233_878_701_836_9 * 2 ** 32),
#     tag="TFHE1024",
# )
# param = estimator.lwe_parameters.LWEParameters(
#    n=512*3,
#    q=2 ** 64,
#    Xs=estimator.nd.NoiseDistribution.UniformMod(3),
#    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=0.0000000000034525330484572114*2**64),
#    tag="TFHE1536"
# )
param = estimator.lwe_parameters.LWEParameters(
   n=2048,
   q=2 ** 64,
   Xs=estimator.nd.NoiseDistribution.UniformMod(2),
   Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=2**(64-47)),
   tag="TFHE2048",
)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
