#!/bin/python3

import  numpy as np
from scipy.special import erf,erfc
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
    l = 3
    Bgbit = 11
    Bg = 2**Bgbit
    α = 2**-44
    ε = 1/(2*(Bg**l))
    β = Bg/2

class lvl10param:
    t = 8
    basebit = 2

class lvl21param:
    t = 10
    basebit = 3

class lvl22param:
    t = 10
    basebit = 3

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
    return highP.n*(2**(-2*(funcP.basebit*funcP.t)+1))+funcP.t*highP.n*(lowP.α**2)

def gbnoisecalc(lowP,highP,funcP):
    return brnoisecalc(lowP,highP)+iknoisecalc(lowP,highP,funcP)

# FA parameter

d = 1600
m = 9

fanoise = d*cmuxnoisecalc(lvl1param,lvl1param.α)
print(fanoise)
print(erfc(1/(2*2*2*np.sqrt(2*fanoise))))
print(np.prod([(1-mpfr(erfc(1/(2*2*2*np.sqrt(2* (i+1) * cmuxnoisecalc(lvl1param,lvl1param.α))))))**(m*lvl1param.n) for i in range(d)]))