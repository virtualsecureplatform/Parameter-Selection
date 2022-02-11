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
    nbit = 9
    k = 3
    n = 2**nbit
    l = 2
    Bgbit = 8
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
    t = 3
    basebit = 4
    domainP = lvl1param
    targetP = lvl0param

class lvl11param:
    t = 6
    basebit = 4
    domainP = lvl1param
    targetP = lvl1param

class lvl21param:
    t = 8
    basebit = 3
    domainP = lvl2param
    targetP = lvl1param

class lvl22param:
    t = 42
    basebit = 1
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

class lvl12param:
    domainP = lvl1param
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
    def funwithjac(numtarr):
        return numccfunc(numtarr),diffccfunc(numtarr)
    return numccfunc, diffccfunc, funwithjac

def externalproduct(P,squaresum,linf,trgswdists,dists):
    dists['normal'] = squaresum*dists['normal']+P.k*P.l*P.n*(P.β**2)*trgswdists["normal"]
    if P.ε in dists['uniform']:
        if linf*P.ε in dists['uniform']:
            dists['uniform'][linf*P.ε] += (1+P.n)*dists['uniform'][P.ε]
        else:
            dists['uniform'][linf*P.ε] = (1+P.n)*dists['uniform'][P.ε]
        del dists['uniform'][P.ε]
    else:
        if linf*P.ε in dists['uniform']:
            dists['uniform'][linf*P.ε] += 1+P.n
        else:
            dists['uniform'][linf*P.ε] = 1+P.n
    for key,value in trgswdists["uniform"].items():
        if P.β*key in dists['uniform']:
            dists['uniform'][P.β*key] += P.k*P.l*P.n*value
        else:
            dists['uniform'][P.β*key] = P.k*P.l*P.n*value

def cmux(P,cmuxdists,dists):
    dists['normal'] += P.k*P.l*P.n*(P.β**2)*cmuxdists["normal"]
    if P.ε in dists['uniform']:
        dists['uniform'][P.ε] += 1+P.n
    else:
        dists['uniform'][P.ε] = 1+P.n
    for key,value in cmuxdists["uniform"].items():
        if P.β*key in dists['uniform']:
            dists['uniform'][P.β*key] += P.k*P.l*P.n*value
        else:
            dists['uniform'][P.β*key] = P.k*P.l*P.n*value

def blindrotate(P,dists):
    cmuxdists = {'normal': P.targetP.α**2,'uniform':{}}
    for i in range(P.domainP.n):
        cmux(P.targetP,cmuxdists,dists)

def identitiykeyswithing(P,dists):
    dists['normal'] += P.t*(P.domainP.k-1)*P.domainP.n*(P.targetP.α**2)
    a = 2**(-P.basebit*P.t-1)
    if a in dists['uniform']:
        dists["uniform"][a] += (P.domainP.k-1)*P.domainP.n
    else:
        dists["uniform"][a] = (P.domainP.k-1)*P.domainP.n

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
        currentdists = dists
        trgswdists = {'normal': P.α**2,'uniform':{}}
        cmux(P,trgswdists,dists)
        dists["normal"] += currentdists["normal"]
        for key,value in currentdists["uniform"].items():
            if key in dists["uniform"]:
                dists["uniform"][key] += value
            else:
                dists["uniform"][key] = value


def GateBootstrapping(brP,ikP,dists):
    blindrotate(brP,dists)
    identitiykeyswithing(ikP,dists)

def CircuitBootstrapping(brP,privksP,dists):
    blindrotate(brP,dists)
    privatekeyswitching(privksP,dists)

def ChensPackingCircuitBootstrapping(brP,dists):
    blindrotate(brP,dists)
    annihilatekeyswitching(brP.targetP,dists)

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
# identitiykeyswithing(lvl10param,dists)
GateBootstrapping(lvl01param,lvl10param,dists)
# GateBootstrapping(lvl12param,lvl21param,dists)
# privatekeyswitching(lvl11param,dists)
# privatekeyswitching(lvl21param,dists)
# privatekeyswitching(lvl22param,dists)
# annihilatekeyswitching(lvl2param,dists)
# CircuitBootstrapping(lvl02param,lvl21param,dists)
# CircuitBootstrapping(lvl12param,lvl21param,dists)
# CircuitBootstrapping(lvl02param,lvl22param,dists)
# ChensPackingCircuitBootstrapping(lvl02param,dists)
# dists = romnoisecalc(lvl01param,lvl10param,lvl11param,ROMaddress)
# dists = romnoisecalc(lvl02param,lvl10param,lvl21param,ROMaddress)
# dists = romnoisecalc(lvl02param,lvl20param,lvl22param,ROMaddress)

# dists = {'normal': 0.0006479548101205879,'uniform': {9.5367431640625e-07: 1000000000}} # bounds=[(1e-3,1e5)],initial_temp=1e4
# dists = {'normal': 0.0006479548101205879,'uniform': {2.3283064365386963e-10: 559583539}} # bounds=[(1e-4,1e5)],initial_temp=1e4
print(dists)
print([((2*key)**2)*value/12 for key,value in dists["uniform"].items()])
conventionalvariance = dists["normal"]+sum([((2*key)**2)*value/12 for key,value in dists["uniform"].items()])
print(conventionalvariance)
from scipy.special import erfc
print(erfc(1/(16*np.sqrt(2*dists["normal"]))))
print(erfc(1/(16*np.sqrt(2*conventionalvariance))))

import math
numccfunc, diffccfunc, funwithjac = ccfunc(1/16,dists)

from scipy.optimize import minimize,shgo,dual_annealing

# result = minimize(fun = numccfunc,x0 = np.array([1e-6]), jac = diffccfunc, method = 'Newton-CG')
# print(result['x'])
# print(result['fun'])
# print(2*math.exp(result['fun']))

result = shgo(numccfunc,bounds=[(1e-6,None)],minimizer_kwargs={'method': "SLSQP", 'jac':diffccfunc})
# result = dual_annealing(numccfunc,bounds=[(1e-3,1e5)],initial_temp=1e4)

print(result.x)
print(result.fun)
print(2*math.exp(result.fun))