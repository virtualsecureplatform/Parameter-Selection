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
Q = 2**64

class lvl0param:
  n = 500
  k = 1
  α = 2.44e-5 * q
  # binary
  variance_key_coefficient = 1./4
  expectation_key_coefficient = 1./2

# class lvl0param:
#   n = 777
#   k = 1
#   α = 5.033523219195911e-06 * q
#   # binary
#   variance_key_coefficient = 1./4
#   expectation_key_coefficient = 1./2

# class lvl0param:
#     nbit = 9
#     n = 2**nbit
#     k = 2
#     l = 2
#     ℬbit = 8
#     ℬ = 2**ℬbit
#     α = 0.0000000342338787018369 * q
#     # ternary
#     variance_key_coefficient = 2./3
#     expectation_key_coefficient = 0

class lvl1param:
    nbit = 10
    n = 2**nbit
    k = 1
    l = 2
    ℬbit = 10
    ℬ = 2**ℬbit
    α = 3.73e-9 * q
    # ternary
    variance_key_coefficient = 1./4
    expectation_key_coefficient = 1./2

class lvl10param:
    t = 7
    basebit = 2

class lvl2param:
    nbit = 10
    k = 2
    n = 2**nbit
    l = 4
    ℬbit = 9
    ℬ = 2**ℬbit
    α = 2**-47
    ε = 1/(2*(ℬ**l))
    β = ℬ/2
    variance_key_coefficient = 1./4
    expectation_key_coefficient = 1./2

class lvl21param:
    t = 24
    basebit = 1
    domainP = lvl2param
    targetP = lvl1param

class lvl20param:
    t = 7
    basebit =  2
    domainP = lvl2param
    targetP = lvl0param

class lvl21mrlweparam:
    nbit = lvl1param.nbit
    n = 2**nbit
    k = lvl1param.k
    l = 4
    ℬbit = 6
    ℬ = 2**ℬbit
    α = lvl1param.α
    # ternary
    variance_key_coefficient = lvl1param.variance_key_coefficient
    expectation_key_coefficient = lvl1param.expectation_key_coefficient



def brnoisecalc(lowP,highP):
    res1 = highP.l * (highP.k + 1.) * highP.n * (highP.ℬ**2 + 2.) / 12. * highP.α**2
    res2 = (q**2-highP.ℬ**(2*highP.l)) / (24 * highP.ℬ**(2*highP.l)) * (1. + highP.k * highP.n * (lowP.variance_key_coefficient + lowP.expectation_key_coefficient**2)) + highP.k * highP.n/8 * lowP.variance_key_coefficient  + 1 / 16. * (1. - highP.k * highP.n * lowP.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return lowP.k * lowP.n * (res1+res2)

def unrollbrnoisecalc(lowP,highP,m):
    res1 = highP.l * (highP.k + 1.) * highP.n * (highP.ℬ**2 + 2.) / 12. * (2**m - 1)*(highP.α)**2
    res2 = (q**2-highP.ℬ**(2*highP.l)) / (24 * highP.ℬ**(2*highP.l)) * (1. + highP.k * highP.n * (lowP.variance_key_coefficient + lowP.expectation_key_coefficient**2)) + highP.k * highP.n/8 * lowP.variance_key_coefficient  + 1 / 16. * (1. - highP.k * highP.n * lowP.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return lowP.k * lowP.n/m * (res1+res2)

def iksnoisecalc(lowP,highP,funcP):
    # return 1/12*highP.k*highP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*highP.k*highP.n*(lowP.α**2)
    roundwidth = 2**(-funcP.basebit*funcP.t-1) * q
    round_variance = roundwidth**2/12 - 1/12
    round_expectation = -1./2
    return highP.k*highP.n*((round_variance*highP.variance_key_coefficient+round_variance*highP.expectation_key_coefficient**2+round_expectation**2 * highP.variance_key_coefficient)+funcP.t*(lowP.α**2))

def privksnoisecalc(lowP,highP,funcP):
    # return 1/12*highP.k*highP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*highP.k*highP.n*(lowP.α**2)
    roundwidth = 2**(-funcP.basebit*funcP.t-1) * q
    round_variance = roundwidth**2/12 - 1/12
    round_expectation = -1./2
    return (highP.k*highP.n+1)*((round_variance*highP.variance_key_coefficient+round_variance*highP.expectation_key_coefficient**2+round_expectation**2 * highP.variance_key_coefficient)+funcP.t*(lowP.α**2))

def mrlweikscalc(lowP,highP):
    res1 = lowP.l * (lowP.k + 1.) * lowP.n * (lowP.ℬ**2 + 2.) / 12. * lowP.α**2
    roundwidth = 2**(-lowP.ℬbit*lowP.l) * q
    round_variance = roundwidth**2/12 - 1/12
    round_expectation = -1./2
    res2 = (q**2-lowP.ℬ**(2*lowP.l)) / (24 * lowP.ℬ**(2*lowP.l)) * (1. + lowP.k * lowP.n * (highP.variance_key_coefficient + highP.expectation_key_coefficient**2)) + lowP.k * lowP.n/8 * highP.variance_key_coefficient  + 1 / 16. * (1. - lowP.k * lowP.n * highP.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return res1+lowP.n*res2

def cmuxnoisecalc(P,α):
    res1 = P.l * (P.k + 1.) * P.n * (P.ℬ**2 + 2.) / 12. * P.α**2
    res2 = (q**2-P.ℬ**(2*P.l)) / (24 * P.ℬ**(2*P.l)) * (1. + P.k * P.n) # Last Part seems to be integer representation specific.
    return res1 + res2

print("Gate")
print("BR noise")
brnoise = brnoisecalc(lvl0param,lvl1param)
print(brnoise)
print(np.sqrt(brnoise)/q)
print(erfc((q/16)/np.sqrt(2*brnoise)))

iksnoise = iksnoisecalc(lvl0param,lvl1param,lvl10param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/q)
print(erfc((q/16)/np.sqrt(2*iksnoise)))

print("Gate Error")
print(erfc((q/16)/np.sqrt(2*(iksnoise+brnoise))))

print("m = 4 BR noise")
brnoise = unrollbrnoisecalc(lvl0param,lvl1param,4)
print(brnoise)
print(np.sqrt(brnoise)/q)
print(erfc((q/16)/np.sqrt(2*brnoise)))

print("Gate Error")
print(erfc((q/16)/np.sqrt(2*(iksnoise+brnoise))))

print("CB")
brnoise = brnoisecalc(lvl0param,lvl2param)
print(brnoise)
print(np.sqrt(brnoise)/Q)

privksnoise = privksnoisecalc(lvl1param,lvl2param,lvl21param)
print("PrivIKS noise")
print(privksnoise)
print(np.sqrt(privksnoise)/Q)
print(erfc((Q/16)/np.sqrt(2*privksnoise)))
print(erfc((Q/16)/np.sqrt(2*brnoise)))

print(erfc((Q/16)/np.sqrt(2*(privksnoise+brnoise))))

print("MRLWE IKS noise")
mrlweiks = mrlweikscalc(lvl21mrlweparam,lvl2param)
print(mrlweiks)
print(privksnoise)
print(erfc((q/16)/np.sqrt(2*mrlweiks)))

def cbnoisecalc(domainP,middleP,targetP,privksP):
    return brnoisecalc(domainP,middleP)+privksnoisecalc(targetP,middleP,privksP)

def romnoisecalc(addressP,dataP,middleP,ikP,privksP,ROMaddress):
    # return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))+iksnoisecalc(addressP,dataP,ikP)
    return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))

ROMaddress = 256 # 4 word block
print("TFHE ROM CMUX noise")
romnoise = romnoisecalc(lvl0param,lvl1param,lvl2param,lvl10param,lvl21param,ROMaddress)
print(romnoise)
iksnoisecalc
print("TFHE ROM error prob")
print(erfc(q/(16*np.sqrt(2*romnoise))))

RAMaddress = 7
RAMwordbit = 8

def ramnoisecalc(addressP,dataP,middleP,ikP,privksP,RAMaddress):
    return brnoisecalc(addressP,dataP)+RAMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))+iksnoisecalc(addressP,dataP,ikP)

print("RAM Read Noise")
rnoise = ramnoisecalc(lvl0param,lvl1param,lvl2param,lvl10param,lvl21param,RAMaddress)
print(rnoise)
print("RAM Read error prob")
print(erfc(1/(16*np.sqrt(2*rnoise))))
print("RAM Read error prob in 3000 cycle")
print(1-(1-mpfr(erfc(1/(16*np.sqrt(2*rnoise)))))**(3000*RAMwordbit))

print("multi-bit LUT")
print("BR noise")
brnoise = brnoisecalc(lvl0param,lvl2param)
print(brnoise)
print(np.sqrt(brnoise)/Q)
print(erfc((Q/64)/np.sqrt(2*brnoise)))

iksnoise = iksnoisecalc(lvl0param,lvl2param,lvl20param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/Q)
print(erfc((q/64)/np.sqrt(2*iksnoise)))

print("Gate Error")
print(erfc((q/64)/np.sqrt(2*(iksnoise+brnoise))))
