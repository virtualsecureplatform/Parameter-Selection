#!/bin/python3

import  numpy as np
from scipy.special import erfc
import gmpy2
from gmpy2 import mpfr

gmpy2.get_context().precision=200

# 128bit TFHE's parameter.
lvl0_μ = lvl1_μ = 2**29;
lvl2_μ = 2**61;

# class lvl0param:
#   n = 512
#   k = 1
#   q = 2**16
#   α = 2**-10.5 * q
#   # binary
#   variance_key_coefficient = 1./4
#   expectation_key_coefficient = 1./2

# class lvl0param:
#   n = 636
#   k = 1
#   q = 2**32
#   α = 0.000_092_511_997_467_675_6 * q
#   # binary
#   variance_key_coefficient = 1./4
#   expectation_key_coefficient = 1./2

class lvl0param:
  n = 128
  k = 5
  q = 2**32
  α = 0.000_092_511_997_467_675_6 * q
  # binary
  variance_key_coefficient = 1./4
  expectation_key_coefficient = 1./2

# class lvl0param:
#   n = 777
#   k = 1
#   q = 2**32
#   α = 0.000003725679281679651 * q
#   # binary
#   variance_key_coefficient = 1./4
#   expectation_key_coefficient = 1./2

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
    nbit = 9
    # q = 2**32
    # q = 5**4*2**16+1
    q = 2**25
    n = 2**nbit
    k = 2
    l = 2
    ℬbit = 8
    ℬ = 2**ℬbit
    # α = 0.0000000342338787018369 * q
    α = 2*4/4
    # ternary
    variance_key_coefficient = 2./3
    expectation_key_coefficient = 0

# class lvl1param:
#     nbit = 9
#     n = 2**nbit
#     k = 2
#     l = 1
#     ℬbit = 16
#     ℬ = 2**ℬbit
#     α = 0.0000000000034525330484572114 * q
#     # ternary
#     variance_key_coefficient = 2./3
#     expectation_key_coefficient = 0

# class lvl1param:
#     nbit = 8
#     n = 2**nbit
#     k = 4
#     l = 2
#     ℬbit = 7
#     ℬ = 2**ℬbit
#     q = (2 ** 16 + 1)*(2 **8 +1)
#     α =  4.2
#     # ternary
#     variance_key_coefficient = 2./3
#     expectation_key_coefficient = 0

# class lvl1param:
#     nbit = 10
#     n = 2**nbit
#     q = 2**32
#     k = 1
#     l = 3
#     ℬbit = 6
#     ℬ = 2**ℬbit
#     α = 2**-25 * q
#     # ternary
#     variance_key_coefficient = 2./3
#     expectation_key_coefficient = 0

class lvl10param:
    t = 2
    basebit = 12

# class lvl10param:
#     t = 4
#     basebit = 3

class lvl2param:
    nbit = 11
    k = 1
    n = 2**nbit
    q = 2**64
    l = 4
    ℬbit = 9
    ℬ = 2**ℬbit
    α =  q * 2**-47
    ε = 1/(2*(ℬ**l))
    β = ℬ/2
    variance_key_coefficient = 1./4
    expectation_key_coefficient = 1./2

# class lvl2param:
#     nbit = 9
#     k = 3
#     n = 2**nbit
#     q = 2**64
#     l = 1
#     ℬbit = 18
#     ℬ = 2**ℬbit
#     α = 0.0000000000034525330484572114 * q
#     ε = 1/(2*(ℬ**l))
#     β = ℬ/2
#     variance_key_coefficient = 1./4
#     expectation_key_coefficient = 1./2

# class lvl2param:
#     nbit = 9
#     k = 3
#     n = 2**nbit
#     q = 2**64
#     l = 3
#     ℬbit = 9
#     ℬ = 2**ℬbit
#     α = 2**-39 * q
#     ε = 1/(2*(ℬ**l))
#     β = ℬ/2
#     variance_key_coefficient = 2./3
#     expectation_key_coefficient = 0.

class lvl11param:
    t = 5
    basebit = 5
    domainP = lvl2param
    targetP = lvl2param

class lvl21param:
    t = 8
    basebit = 3
    domainP = lvl2param
    targetP = lvl1param

class lvl22param:
    t = 38
    basebit = 1
    domainP = lvl2param
    targetP = lvl2param

class lvl20param:
    t = 7
    basebit =  2
    domainP = lvl2param
    targetP = lvl0param

class lvl21mrlweparam:
    nbit = lvl1param.nbit
    n = 2**nbit
    k = lvl1param.k
    q = 2**32
    l = 4
    ℬbit = 6
    ℬ = 2**ℬbit
    α = lvl1param.α
    # ternary
    variance_key_coefficient = lvl1param.variance_key_coefficient
    expectation_key_coefficient = lvl1param.expectation_key_coefficient

def extpnoisecalc(P,α):
    res1 = P.l * (P.k + 1.) * P.n * (P.ℬ**2 + 2.) / 12. * α**2
    res2 = P.ℬ**2/2
    res3 = (P.q**2-P.ℬ**(2*P.l)) / (24 * P.ℬ**(2*P.l)) * (1. + P.k * P.n * (P.variance_key_coefficient + P.expectation_key_coefficient**2))
    res4 = P.k * P.n/8 * P.variance_key_coefficient
    res5 = 1 / 16. * (1. - P.k * P.n * P.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return res1 + res2 + res3 + res4 + res5

def brnoisecalc(lowP,highP):
    return lowP.k * lowP.n * extpnoisecalc(highP,highP.α)

def unrollbrnoisecalc(lowP,highP,m):
    res1 = highP.l * (highP.k + 1.) * highP.n * (highP.ℬ**2 + 2.) / 12. * (2**m - 1)*(highP.α)**2
    res2 = (highP.q**2-highP.ℬ**(2*highP.l)) / (24 * highP.ℬ**(2*highP.l)) * (1. + highP.k * highP.n * (lowP.variance_key_coefficient + lowP.expectation_key_coefficient**2)) + highP.k * highP.n/8 * lowP.variance_key_coefficient  + 1 / 16. * (1. - highP.k * highP.n * lowP.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return lowP.k * lowP.n/m * (res1+res2)

def iksnoisecalc(lowP,highP,funcP):
    # return 1/12*highP.k*highP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*highP.k*highP.n*(lowP.α**2)
    roundwidth = 2**(-funcP.basebit*funcP.t-1) * lowP.q
    round_variance = (2*roundwidth)**2/12 - 1/12
    round_expectation = -1./2
    return highP.k*highP.n*((round_variance*highP.variance_key_coefficient+round_variance*highP.expectation_key_coefficient**2+round_expectation**2 * highP.variance_key_coefficient)+funcP.t*(lowP.α**2))

def privksnoisecalc(lowP,highP,funcP):
    # return 1/12*highP.k*highP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*highP.k*highP.n*(lowP.α**2)
    roundwidth = 2**(-funcP.basebit*funcP.t-1) * lowP.q
    round_variance = roundwidth**2/12 - 1/12
    round_expectation = -1./2
    return (highP.k*highP.n+1)*((round_variance*highP.variance_key_coefficient+round_variance*highP.expectation_key_coefficient**2+round_expectation**2 * highP.variance_key_coefficient)+funcP.t*(lowP.α**2))

def mrlweikscalc(lowP,highP):
    res1 = lowP.l * (lowP.k + 1.) * lowP.n * (lowP.ℬ**2 + 2.) / 12. * lowP.α**2
    roundwidth = 2**(-lowP.ℬbit*lowP.l) * lowP.q
    round_variance = roundwidth**2/12 - 1/12
    round_expectation = -1./2
    res2 = (lowP.q**2-lowP.ℬ**(2*lowP.l)) / (24 * lowP.ℬ**(2*lowP.l)) * (1. + lowP.k * lowP.n * (highP.variance_key_coefficient + highP.expectation_key_coefficient**2)) + lowP.k * lowP.n/8 * highP.variance_key_coefficient  + 1 / 16. * (1. - lowP.k * lowP.n * highP.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return res1+lowP.n*res2

def cmuxnoisecalc(P,α):
    res1 = P.l * (P.k + 1.) * P.n * (P.ℬ**2 + 2.) / 12. * α**2
    res2 = P.ℬ**2/2
    res3 = (P.q**2-P.ℬ**(2*P.l)) / (24 * P.ℬ**(2*P.l)) * (1. + P.k * P.n * (P.variance_key_coefficient + P.expectation_key_coefficient**2))
    res4 = P.k * P.n/8 * P.variance_key_coefficient
    res5 = 1 / 16. * (1. - P.k * P.n * P.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return res1 + res2 + res3 + res4 + res5

def brroundnoise(domainP,targetP):
    roundwidth = targetP.q/(4*targetP.n)
    round_variance = (2*roundwidth)**2/12 - 1/12
    round_expectation = -1./2
    return domainP.k*domainP.n*(round_variance*domainP.variance_key_coefficient+round_variance*domainP.expectation_key_coefficient**2+round_expectation**2 * domainP.variance_key_coefficient)
    # return domainP.k*domainP.n*1/4*targetP.n*domainP.q

def annihilaterecursive(P,nbit,α):
    if(nbit==0):
        return α
    else:
        return 2*annihilaterecursive(P,nbit-1,α)+extpnoisecalc(P,P.α)

def annihilatecalc(P,α):
    return annihilaterecursive(P,P.nbit,α)

print("Gate")
print("BR noise")
brnoise = brnoisecalc(lvl0param,lvl1param)
print(brnoise)
print(np.sqrt(brnoise)/lvl1param.q)
print(erfc((lvl1param.q/16)/np.sqrt(2*brnoise)))

iksnoise = iksnoisecalc(lvl0param,lvl1param,lvl10param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/lvl0param.q)
print(erfc((lvl0param.q/8)/np.sqrt(2*iksnoise)))

print("Gate Error")
print(erfc((lvl0param.q/8)/np.sqrt(2*(iksnoise+2*brnoise*((lvl0param.q/lvl1param.q)**2)))))

print("BR Round Noise")
roundnoise = brroundnoise(lvl0param,lvl1param)
print(roundnoise)
print(erfc((lvl0param.q/8)/np.sqrt(2*(roundnoise+iksnoise+2*brnoise*(lvl0param.q/lvl1param.q)**2))))

print("m = 2 BR noise")
brnoise = unrollbrnoisecalc(lvl0param,lvl1param,2)
print(brnoise)
print(np.sqrt(brnoise)/lvl0param.q)
print(erfc((lvl2param.q/16)/np.sqrt(2*brnoise)))

print("m=2 Gate Error")
print(erfc((lvl0param.q/16)/np.sqrt(2*(iksnoise+brnoise))))

print("lvl02 Gate")
print("BR noise")
brnoise = brnoisecalc(lvl0param,lvl2param)
print(brnoise)
print(np.sqrt(brnoise)/lvl2param.q)
print(erfc((lvl2param.q/16)/np.sqrt(2*brnoise)))

iksnoise = iksnoisecalc(lvl0param,lvl2param,lvl20param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/lvl0param.q)
print(erfc((lvl0param.q/16)/np.sqrt(2*iksnoise)))

print("Gate Error")
print(brnoise*(lvl0param.q/lvl2param.q)**2)
print(erfc((lvl0param.q/16)/np.sqrt(2*(iksnoise+brnoise*(lvl0param.q/lvl2param.q)**2))))

print("BRlvl02")
brnoise = brnoisecalc(lvl0param,lvl2param)
print(brnoise)
print(np.sqrt(brnoise)/lvl2param.q)
brnoise *= ((lvl1param.q/lvl2param.q)**2)
print(erfc((lvl1param.q/16)/np.sqrt(2*brnoise)))

privksnoise = privksnoisecalc(lvl1param,lvl2param,lvl21param)
print("PrivIKS noise")
print(privksnoise)
print(np.sqrt(privksnoise)/lvl1param.q)
print(erfc((lvl1param.q/16)/np.sqrt(2*privksnoise)))

print("CB")
print(erfc((lvl1param.q/16)/np.sqrt(2*(privksnoise+brnoise))))

print("MRLWE IKS noise")
mrlweiks = mrlweikscalc(lvl21mrlweparam,lvl2param)
print(mrlweiks)
print(privksnoise)
print(erfc((lvl0param.q/16)/np.sqrt(2*mrlweiks)))

def cbnoisecalc(domainP,middleP,targetP,privksP):
    return brnoisecalc(domainP,middleP)*((targetP.q/middleP.q)**2)+privksnoisecalc(targetP,middleP,privksP)
print(brnoise+privksnoise-cbnoisecalc(lvl0param,lvl2param,lvl1param,lvl21param))

def romnoisecalc(addressP,dataP,middleP,ikP,privksP,ROMaddress):
    # return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))+iksnoisecalc(addressP,dataP,ikP)
    return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))

ROMaddress = 40 # 4 word block
print("TFHE ROM CMUX noise")
romnoise = romnoisecalc(lvl0param,lvl1param,lvl2param,lvl10param,lvl21param,ROMaddress)
print(romnoise)
print(np.sqrt(romnoise)/lvl1param.q)
print("TFHE ROM error prob")
print(erfc(lvl1param.q/(4*np.sqrt(2*romnoise))))

brnoise = brnoisecalc(lvl0param,lvl1param)

privksnoise = privksnoisecalc(lvl1param,lvl1param,lvl11param)
print("PrivIKS lvl11 noise")
print(privksnoise)
print(np.sqrt(privksnoise)/lvl1param.q)
print(erfc((lvl1param.q/16)/np.sqrt(2*privksnoise)))

print("CBlvl01")
print(erfc((lvl1param.q/16)/np.sqrt(2*(privksnoise+brnoise))))
print("CBlvl01 extp")
extpnoise = extpnoisecalc(lvl1param,np.sqrt(privksnoise+brnoise))#+brnoisecalc(lvl0param,lvl1param)
print(erfc((lvl1param.q/16)/np.sqrt(2*(extpnoise))))

print("Annihilate CB lvl22")
brnoise = brnoisecalc(lvl0param,lvl2param)
annihilatenoise = annihilatecalc(lvl2param,brnoise)
print(annihilatenoise)
annihilatecbnoise = annihilatenoise + extpnoisecalc(lvl2param,lvl2param.α)
print(erfc((lvl2param.q/16)/np.sqrt(2*(annihilatecbnoise))))

print("TFHE ROM CMUX noise lvl22")
romnoise = romnoisecalc(lvl0param,lvl2param,lvl2param,lvl20param,lvl22param,ROMaddress)
print(romnoise)
print(np.sqrt(romnoise)/lvl2param.q)
print("TFHE ROM error prob lvl22")
print(erfc(lvl2param.q/(4*np.sqrt(2*romnoise))))

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
print(np.sqrt(brnoise)/lvl2param.q)
print(erfc((lvl2param.q/64)/np.sqrt(2*brnoise)))

iksnoise = iksnoisecalc(lvl0param,lvl2param,lvl20param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/lvl2param.q)
print(erfc((lvl0param.q/64)/np.sqrt(2*iksnoise)))

print("Gate Error")
print(erfc((lvl0param.q/64)/np.sqrt(2*(iksnoise+brnoise))))
