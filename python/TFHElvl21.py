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
    t = 24
    basebit = 1

class lvl22param:
    t = 10
    basebit = 3

class lvl20param:
    t = 7
    basebit =  2

# Rounding error from Decomposition of External Product should be treated as Irwin-Hall(though not upper bound)
# https://tches.iacr.org/index.php/TCHES/article/view/8793
def cmuxnoisecalc(P,α):
    # return 2*P.l*P.n*(P.β**2)*(α**2)+1/12*(1+P.n)*((2*P.ε)**2)
    return 2*P.l*P.n*(P.β**2)*(α**2)+(1+P.n)*(P.ε**2)

def brnoisecalc(lowP,highP):
    return lowP.n*cmuxnoisecalc(highP,highP.α)

def iknoisecalc(lowP,highP,funcP):
    # return 1/12*highP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*highP.n*(lowP.α**2)
    return highP.n*(2**(-2*(funcP.basebit*funcP.t+1)))+funcP.t*highP.n*(lowP.α**2)

def gbnoisecalc(lowP,highP,funcP):
    return brnoisecalc(lowP,highP)+iknoisecalc(lowP,highP,funcP)

# Processor parameter

ROMaddress = 7 # 4 word block
RAMaddress = 7
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
    return (domainP.n+1)*(2**(-2*(privksP.basebit*privksP.t+1)))+privksP.t*(domainP.n+1)*(targetP.α**2)

def cbnoisecalc(domainP,middleP,targetP,privksP):
    return brnoisecalc(domainP,middleP)+privksnoisecalc(middleP,targetP,privksP)

cbnoise = cbnoisecalc(lvl0param,lvl2param,lvl1param,lvl21param)

print("TFHE lvl20 iks noise")
print(iknoisecalc(lvl0param,lvl2param,lvl20param))

print("BR lvl 02")
print(brnoisecalc(lvl0param,lvl2param))
print("privks lvl21")
print(privksnoisecalc(lvl2param,lvl1param,lvl21param))
print("TFHE Circuit Bootstrapping lvl21 Noise")
print(cbnoise)

def romnoisecalc(addressP,dataP,middleP,ikP,privksP,ROMaddress):
    return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))+iknoisecalc(addressP,dataP,ikP)

print("TFHE ROM CMUX noise")
romnoise = romnoisecalc(lvl0param,lvl1param,lvl2param,lvl10param,lvl21param,ROMaddress)
print(romnoise)

print("TFHE ROM error prob")
print(erfc(1/(16*np.sqrt(2*romnoise))))

def ramnoisecalc(addressP,dataP,middleP,ikP,privksP,RAMaddress):
    return brnoisecalc(addressP,dataP)+RAMaddress*np.sqrt(cmuxnoisecalc(dataP,cbnoisecalc(addressP,middleP,dataP,privksP)))+iknoisecalc(addressP,dataP,ikP)

print("RAM Read Noise")
rnoise = ramnoisecalc(lvl0param,lvl1param,lvl2param,lvl10param,lvl21param,RAMaddress)
print(rnoise)
print("RAM Read error prob")
print(erfc(1/(16*np.sqrt(2*rnoise))))
print("RAM Read error prob in 3000 cycle")
print(1-(1-mpfr(erfc(1/(16*np.sqrt(2*rnoise)))))**(3000*RAMwordbit))

# print("Packing Switch noise")


# psnoise = DEF_n*2*DEF_l*DEF_N*(DEF_β**2)*(DEF_αbk**2)+DEF_n*(1+DEF_N)*(DEF_ε**2) + DEF_N*(2**(-2*(DEF_basebit*DEF_t+1)))+8*DEF_t*DEF_N*(DEF_αbk**2)
# print(psnoise)

# print("significant term")
# print(DEF_n*2*DEF_lbar*DEF_nbar*(DEF_βbar**2)*(DEF_αbklvl02**2))
# print(DEF_n*(1+DEF_nbar)*(DEF_εbar**2))
# print(DEF_nbar*(2**(-2*(DEF_basebitlvl21*DEF_tbar+1))))
# print(DEF_tbar*DEF_nbar*(DEF_αprivks**2))
