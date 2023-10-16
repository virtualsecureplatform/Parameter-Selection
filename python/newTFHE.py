#!/bin/python3
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")
# param = estimator.lwe_parameters.LWEParameters(
#     n=511,
#     q=2 ** 16,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=2 ** (16-10.5)),
#     # Xe=estimator.nd.NoiseDistribution.UniformMod(8),
#     tag="ShortLWE",
# )
# param = estimator.lwe_parameters.LWEParameters(
#     n=480,
#     q=2 ** 8,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=3.19),
#     # Xe=estimator.nd.NoiseDistribution.UniformMod(8),
#     tag="ShortLWE",
# )
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
#    n=640,
#    q=2 ** 32,
#    Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=0.000_092_511_997_467_675_6 * 2 ** 32),
#    tag="TFHE640",
# )
# param = estimator.lwe_parameters.LWEParameters(
#     n=777,
#     q=2 ** 32,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=2 ** (32-17.6)),
#     tag="TFHE586",
# )
# param = estimator.lwe_parameters.LWEParameters(
#     n=760,
#     q=2 ** 32,
#     Xs=estimator.nd.NoiseDistribution.UniformMod(2),
#     # Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 2 ** 7),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 2 ** (32-17)),
#     tag="TFHE760",
# )
param = estimator.lwe_parameters.LWEParameters(
    n=1024,
    q= 5**4*2**16+1,
    # q=2 ** 32,
   #  q=2 ** 25,
    Xs=estimator.nd.NoiseDistribution.UniformMod(3),
    # Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 2 ** 7),
    # Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 4.2),
    # Xe=estimator.nd.NoiseDistribution.CenteredBinomial(4),
    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 2*4/4),
    # Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 0.000_000_034_233_878_701_836_9 * 2 ** 32),
    tag="TFHE1024",
)
# param = estimator.lwe_parameters.LWEParameters(
#     n=768,
#     q= (2 ** 16 + 1)*(2 **8 +1),
#     Xs=estimator.nd.NoiseDistribution.UniformMod(3),
#     # Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 2 ** 7),
#     Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev= 8),
#     tag="TFHE1024",
# )
# param = estimator.lwe_parameters.LWEParameters(
#    n=512*3,
#    q=2 ** 64,
#    Xs=estimator.nd.NoiseDistribution.UniformMod(3),
# #    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=0.0000000000034525330484572114*2**64),
#    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=2**(64-39)),
#    tag="TFHE1536"
# )
# param = estimator.lwe_parameters.LWEParameters(
#    n=2048,
# #    q=2 ** 64,
#    q=2 ** 48,
#    Xs=estimator.nd.NoiseDistribution.UniformMod(3),
# #    Xe=estimator.nd.NoiseDistribution.DiscreteGaussian(stddev=2**(64-47)),
#    Xe=estimator.nd.NoiseDistribution.CenteredBinomial(3),
#    tag="TFHE2048",
# )
print(param.n)
r = estimator.LWE.estimate(param,red_cost_model = estimator.RC.BDGL16)
