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

def externalproduct(P):
    return 2*P.l*P.n*(P.β**2)*(P.α**2)+(1+P.n)*P.n*(P.ε**2)

def exprivkscalc(domainP,targetP,privksP):
    return 1/12*(2**(-2*(privksP.basebit*privksP.t))) + externalproduct(targetP)

def brnoisecalc(lowP,highP):
    return lowP.n*cmuxnoisecalc(highP,highP.α)

# https://tches.iacr.org/index.php/TCHES/article/view/8793
def iknoisecalc(domainP,targetP,funcP):
    return 1/12*domainP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*domainP.n*(targetP.α**2)

def gbnoisecalc(lowP,highP,funcP):
    return brnoisecalc(lowP,highP)+iknoisecalc(highP,lowP,funcP)

# Processor parameter

ROMaddress = 7 # 4 word block
RAMaddress = 9
RAMwordbit = 8

print("TFHE IK noise")
print(iknoisecalc(lvl1param,lvl0param,lvl10param))

print("BR noise")
print(brnoisecalc(lvl0param,lvl1param))

print("TFHE Gate Bootstrapping Noise")

gbnoise  = gbnoisecalc(lvl0param,lvl1param,lvl10param)

print(gbnoise)

print("TFHE GB error prob")
print(erfc(1/(16*np.sqrt(2*gbnoise))))

def annihilatecalc(P):
    return annihilaterecursive(P,P.nbit)

def annihilaterecursive(P,nbit):
    if(nbit==0):
        return 0
    else:
        return 2*annihilaterecursive(P,nbit-1)+P.l*P.n*(P.β**2)*(P.α**2)+P.n*(P.ε**2)

print("Annihilate")
print(annihilatecalc(Annihilatelvl1param))
# print(annihilatecalc(lvl2param))
print(erfc(1/(4*np.sqrt(2*100*(annihilatecalc(Annihilatelvl1param)+brnoisecalc(lvl0param,lvl1param))))))

# https://tches.iacr.org/index.php/TCHES/article/view/8793
def privksnoisecalc(domainP,targetP,privksP):
    return 1/12*(domainP.n+1)*(2**(-2*(privksP.basebit*privksP.t)))+privksP.t*(domainP.n+1)*(targetP.α**2)

def cbnoisecalc(domainP,middleP,targetP,privksP):
    return brnoisecalc(domainP,middleP)+privksnoisecalc(middleP,targetP,privksP)

cbnoise = cbnoisecalc(lvl0param,lvl2param,lvl2param,lvl22param)

print("TFHE lvl20 iks noise")
print(iknoisecalc(lvl2param,lvl0param,lvl20param))

print("TFHE lvl22 privks noise")
print(privksnoisecalc(lvl2param,lvl2param,lvl22param))
print(externalproduct(lvl2param))
print(exprivkscalc(lvl0param,lvl2param,lvl22param))

print("TFHE Circuit Bootstrapping lvl22 Noise")
print(cbnoise)

def romnoisecalc(addressP,dataP,middleP,ikP,privksP,ROMaddress):
    return dataP.α+ROMaddress*cmuxnoisecalc(dataP,cbnoisecalc(addressP,middleP,dataP,privksP))+iknoisecalc(dataP,addressP,ikP)

print("TFHE ROM CMUX noise")
romnoise = romnoisecalc(lvl0param,lvl2param,lvl2param,lvl20param,lvl22param,ROMaddress)
print(romnoise)

print("TFHE ROM error prob")
print(erfc(1/(16*np.sqrt(2*romnoise))))

def ramnoisecalc(addressP,dataP,middleP,ikP,privksP,RAMaddress):
    return brnoisecalc(addressP,dataP)+RAMaddress*cmuxnoisecalc(dataP,cbnoisecalc(addressP,middleP,dataP,privksP))+iknoisecalc(dataP,addressP,ikP)

print("RAM Read Noise")
rnoise = ramnoisecalc(lvl0param,lvl2param,lvl2param,lvl20param,lvl22param,RAMaddress)
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
