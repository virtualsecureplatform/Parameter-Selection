#!/bin/python3
from sympy import integrate,Symbol,pi,oo,log,diff
from sympy.functions import exp, sqrt
from sympy.printing.julia import JuliaCodePrinter
from sympy.printing.octave import print_octave_code

import numpy as np

# 128bit TFHE's parameter.

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
    domainP = lvl1param
    targetP = lvl0param

class lvl11param:
    t = 6
    basebit = 4
    domainP = lvl1param
    targetP = lvl1param

class lvl21param:
    t = 10
    basebit = 3
    domainP = lvl2param
    targetP = lvl1param

class lvl22param:
    t = 8
    basebit = 4
    domainP = lvl2param
    targetP = lvl2param


class lvl20param:
    t = 7
    basebit =  2
    domainP = lvl2param
    targetP = lvl0param

class lvl01param:
    domainP = lvl0param
    targetP = lvl1param

class lvl02param:
    domainP = lvl0param
    targetP = lvl2param


ROMaddress = 7 # 4 word block
RAMaddress = 9
RAMwordbit = 8

def ccfunc(μ,dists):
    x = Symbol('x')
    t = Symbol('t')
    func = -t*μ
    # func = exp(-t*μ)
    if "uniform" in dists:
        for interval, num in dists["uniform"].items():
            func += num*log(integrate(exp(t*x)/(2*interval),(x,-interval,interval)))
            # func *= integrate(exp(t*x)/(2*interval),(x,-interval,interval))**num
    if "normal" in dists:
        # func += log(integrate(exp(t*x)/(sqrt(2*pi*dists["normal"]))*exp(-(x**2)/(2*dists["normal"])),(x,-oo,oo)))
        func += (t**2)*dists["normal"]/2
        # func *= exp((t**2)*dists["normal"]/2)
    # print(func)
    def numccfunc(numtarr):
        return np.array([func.subs([(t,numt)]).evalf() for numt in numtarr],dtype=np.float64)
        #  return func.subs([(t,numtarr)]).evalf()
    def diffccfunc(numtarr):
        return np.array([diff(func,t).subs([(t,numt)]).evalf() for numt in numtarr],dtype=np.float64)
        # return diff(func,t).subs([(t,numtarr)]).evalf()
    return numccfunc, diffccfunc

def cmux(P,cmuxdists,dists):
    dists['normal'] += 2*P.l*P.n*(P.β**2)*cmuxdists["normal"]
    if P.ε in dists['uniform']:
        dists['uniform'][P.ε] += 1+P.n
    else:
        dists['uniform'][P.ε] = 1+P.n
    for key,value in cmuxdists["uniform"].items():
        if P.β*key in dists['uniform']:
            dists['uniform'][P.β*key] += 2*P.l*P.n*value
        else:
            dists['uniform'][P.β*key] = 2*P.l*P.n*value

def blindrotate(P,dists):
    cmuxdists = {'normal': P.targetP.α**2,'uniform':{}}
    for i in range(P.domainP.n):
        cmux(P.targetP,cmuxdists,dists)

def identitiykeyswithing(P,dists):
    dists['normal'] += P.t*P.domainP.n*(P.targetP.α**2)
    a = 2**(-P.basebit*P.t-1)
    if a in dists['uniform']:
        dists["uniform"][a] += P.domainP.n
    else:
        dists["uniform"][a] = P.domainP.n

def privatekeyswitching(P,dists):
    dists['normal'] += P.t*(P.domainP.n+1)*(P.targetP.α**2)
    a = 2**(-P.basebit*P.t-1)
    if a in dists['uniform']:
        dists["uniform"][a] += P.domainP.n + 1
    else:
        dists["uniform"][a] = P.domainP.n + 1

def annihilatekeyswitching(P,dists):
    annihilaterecursive(P,P.nbit,dists)

def annihilaterecursive(P,nbit,dists):
    if(nbit != 0):
        annihilaterecursive(P,nbit-1,dists)
        dists["normal"] *= 2
        for key in dists["uniform"].keys():
            dists["uniform"][key] *= 2
        cmuxdists = {'normal': P.α**2,'uniform':{}}
        cmux(P,cmuxdists,dists)


def GateBootstrapping(brP,ikP,dists):
    blindrotate(brP,dists)
    identitiykeyswithing(ikP,dists)

def CircuitBootstrapping(brP,privksP,dists):
    blindrotate(brP,dists)
    privatekeyswitching(privksP,dists)

def romnoisecalc(brP,ikP,privksP,ROMaddress):
    dists = {'normal': 0,'uniform':{}}
    dists['normal'] += privksP.targetP.α**2

    cbdists = {'normal': 0,'uniform':{}}
    CircuitBootstrapping(brP,privksP,cbdists)
    for i in range(ROMaddress):
        cmux(privksP.targetP,cbdists,dists)
    identitiykeyswithing(ikP,dists)
    return dists

dists = {'normal': 0,'uniform':{}}

# blindrotate(lvl01param,dists)
# blindrotate(lvl02param,dists)
# GateBootstrapping(lvl01param,lvl10param,dists)
# privatekeyswitching(lvl11param,dists)
# privatekeyswitching(lvl22param,dists)
# CircuitBootstrapping(lvl02param,lvl21param,dists)
# CircuitBootstrapping(lvl02param,lvl22param,dists)
# dists = romnoisecalc(lvl01param,lvl10param,lvl11param,ROMaddress)
# dists = romnoisecalc(lvl02param,lvl10param,lvl21param,ROMaddress)
dists = romnoisecalc(lvl02param,lvl20param,lvl22param,ROMaddress)

print(dists)
print([((2*key)**2)*value/12 for key,value in dists["uniform"].items()])
print(dists["normal"]+sum([((2*key)**2)*value/12 for key,value in dists["uniform"].items()]))

import math
numccfunc, diffccfunc = ccfunc(1/16,dists)

from scipy.optimize import minimize

result = minimize(fun = numccfunc,x0 = np.array([7]), jac = diffccfunc, method = 'Newton-CG')
print(result['x'])
print(result['fun'])
print(2*math.exp(result['fun']))