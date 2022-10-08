#!/bin/python3

from humanfriendly import round_number
import  numpy as np
from scipy.special import erfc
import gmpy2
from gmpy2 import mpfr

gmpy2.get_context().precision=200

# 128bit TFHE's parameter.
lvl0_μ = lvl1_μ = 2**29;
lvl2_μ = 2**61;
q = 2**32

class lvl0param:
    n = 586
    α = 0.000_092_511_997_467_675_6 * q
    # binary
    variance_key_coefficient = 1./4
    expectation_key_coefficient = 1./2

class lvl1param:
    nbit = 9
    n = 2**nbit
    k = 2
    l = 2
    ℬbit = 8
    ℬ = 2**ℬbit
    α = 0.0000000342338787018369 * q
    # ternary
    variance_key_coefficient = 2./3
    expectation_key_coefficient = 0

class lvl10param:
    t = 3
    basebit = 4

# class lvl0param:
#     n = 636
#     α = 0.000_092_511_997_467_675_6 * q
#     # binary
#     variance_key_coefficient = 1./4
#     expectation_key_coefficient = 1./2

# class lvl1param:
#     nbit = 10
#     n = 2**nbit
#     k = 1
#     l = 3
#     ℬbit = 6
#     ℬ = 2**ℬbit
#     α = 2 ** (32-25)
#     # ternary
#     variance_key_coefficient = 2./3
#     expectation_key_coefficient = 0

# class lvl10param:
#     t = 4
#     basebit = 3



def brnoisecalc(lowP,highP):
    res1 = highP.l * (highP.k + 1.) * highP.n * (highP.ℬ**2 + 2.) / 12. * highP.α**2
    res2 = (q**2-highP.ℬ**(2*highP.l)) / (24 * highP.ℬ**(2*highP.l)) * (1. + highP.k * highP.n * (lowP.variance_key_coefficient + lowP.expectation_key_coefficient**2)) + highP.k * highP.n/8 * lowP.variance_key_coefficient  + 1 / 16. * (1. - highP.k * highP.n * lowP.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return lowP.n * (res1+res2)

def iksnoisecalc(lowP,highP,funcP):
    # return 1/12*highP.k*highP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*highP.k*highP.n*(lowP.α**2)
    roundwidth = 2**(-funcP.basebit*funcP.t-1) * q
    round_variance = roundwidth**2/12 - 1/12
    round_expectation = -1./2
    return highP.k*highP.n*((round_variance*highP.variance_key_coefficient+round_variance*highP.expectation_key_coefficient**2+round_expectation**2 * highP.variance_key_coefficient)+funcP.t*(lowP.α**2))

brnoise = brnoisecalc(lvl0param,lvl1param)
print(brnoise)
print(np.sqrt(brnoise)/q)

iksnoise = iksnoisecalc(lvl0param,lvl1param,lvl10param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/q)
print(erfc((q/16)/np.sqrt(2*iksnoise)))
print(erfc((q/16)/np.sqrt(2*brnoise)))

print(erfc((q/16)/np.sqrt(2*(iksnoise+brnoise))))
