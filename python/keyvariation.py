#!/bin/python3

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
    n = 635
    α = 2**-15

class lvl1param:
    nbit = 10
    n = 2**nbit
    k = 1
    l = 3
    Bgbit = 6
    Bg = 2**Bgbit
    α = 2**-25
    ε = 1/(2*(Bg**l))
    β = Bg/2

class Annihilatelvl1param:
    nbit = 10
    n = 2**nbit
    l = 12
    Bgbit = 2
    Bg = 2**Bgbit
    α = 2**-25
    ε = 1/(2*(Bg**l))
    β = Bg/2

class lvl2param:
    nbit = 11
    n = 2**nbit
    k = 1
    l = 4
    Bgbit = 9
    Bg = 2**Bgbit
    α = 2**-44
    ε = 1/(2*(Bg**l))
    β = Bg/2

class lvl10param:
    t = 7
    basebit = 2

class lvl21param:
    t = 8
    basebit = 3

class lvl22param:
    t = 10
    basebit = 3

class lvl20param:
    t = 7
    basebit =  2

variance_key_coefficient = 1./4
expectation_key_coefficient = 1./2

def brnoisecalc(lowP,highP):
    res1 = highP.l * (highP.k + 1.) * highP.n * (highP.Bg**2 + 2.) / 12. * highP.α**2
    res2 = (highP.ϵ**2) / 3 * (1. + highP.k * highP.n * (variance_key_coefficient + expectation_key_coefficient**2)) + highP.k * highP.n * (highP.ϵ**2) * variance_key_coefficient # + (highP.ϵ**2) / 4. * (1. - highP.k * highP.n * expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return lowP.n * (res1+res2)

print(brnoisecalc(lvl0param,lvl1param))