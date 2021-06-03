#!/bin/python3

import  numpy as np
from scipy.special import erfc
import gmpy2
from gmpy2 import mpfr

gmpy2.get_context().precision=200

# 128bit TFHE's parameter.
lvl0_μ = lvl1_μ = 2**29;
lvl2_μ = 2**61;

class lvl0param:
    n = 635
    α = 2**-15

class lvl1param:
    nbit = 10
    n = 2**nbit
    l = 3
    Bgbit = 6
    Bg = 2**Bgbit
    α = 2**-25
    ε = 1/(2*(Bg**l))
    β = Bg/2

class lvl2param:
    nbit = 11
    n = 2**nbit
    l = 4
    Bgbit = 9
    Bg = 2**Bgbit
    α = 2**-44
    ε = 1/(2*(Bg**l))
    β = Bg/2

class lvl10param:
    t = 8
    basebit = 2

class lvl21param:
    t = 7
    basebit = 4

class lvl22param:
    t = 8
    basebit = 4

class lvl20param:
    t = 7
    basebit =  2

# Rounding error from Decomposition of External Product should be treated as Irwin-Hall like Keyswitching.
def cmuxnoisecalc(P,α):
    # return 2*P.l*P.n*(P.β**2)*(α**2)+1/12*(1+P.n)*((2*P.ε)**2)
    return 2*P.l*P.n*(P.β**2)*(α**2)+(1+P.n)*(P.ε**2)

def brnoisecalc(lowP,highP):
    return lowP.n*cmuxnoisecalc(highP,highP.α)

# https://tches.iacr.org/index.php/TCHES/article/view/8793
def iknoisecalc(lowP,highP,funcP):
    return 1/12*highP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*highP.n*(lowP.α**2)

def gbnoisecalc(lowP,highP,funcP):
    return brnoisecalc(lowP,highP)+iknoisecalc(lowP,highP,funcP)

# Processor parameter

ROMaddress = 7 # 4 word block
RAMaddress = 9
RAMwordbit = 8

print("TFHE IK noise")
print(iknoisecalc(lvl0param,lvl1param,lvl10param))

print("TFHE Gate Bootstrapping Noise")

gbnoise  = gbnoisecalc(lvl0param,lvl1param,lvl10param)

print(gbnoise)

print("TFHE GB error prob")
print(erfc(1/(16*np.sqrt(2*gbnoise))))

# https://tches.iacr.org/index.php/TCHES/article/view/8793
def privksnoisecalc(domainP,targetP,privksP):
    return 1/12*(domainP.n+1)*(2**(-2*(privksP.basebit*privksP.t)))+privksP.t*(domainP.n+1)*(targetP.α**2)

def cbnoisecalc(domainP,middleP,targetP,privksP):
    return brnoisecalc(domainP,middleP)+privksnoisecalc(middleP,targetP,privksP)

manyLUTrounded = 2*gbnoise + lvl0param.n*(1/(2*lvl2param.n)**2)
print("manyLUTrounded")
print(manyLUTrounded)
print(erfc(1/(4*np.sqrt(2*manyLUTrounded))))