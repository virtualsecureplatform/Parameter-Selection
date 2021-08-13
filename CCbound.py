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
    l = 3
    Bgbit = 11
    Bg = 2**Bgbit
    α = 2**-44
    ε = 1/(2*(Bg**l))
    β = Bg/2

class lvl10param:
    t = 8
    basebit = 2
    domainP = lvl1param
    targetP = lvl0param

class lvl21param:
    t = 10
    basebit = 3
    domainP = lvl2param
    targetP = lvl1param

class lvl22param:
    t = 10
    basebit = 3
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
        return np.array([func.subs([(t,numt)]).evalf() for numt in numtarr])
        #  return func.subs([(t,numtarr)]).evalf()
    def diffccfunc(numtarr):
        return np.array([diff(func,t).subs([(t,numt)]).evalf() for numt in numtarr])
        # return diff(func,t).subs([(t,numtarr)]).evalf()
    return numccfunc, diffccfunc

dists = {'normal': 0,'uniform':{}}

def cmux(P,α,dists):
    dists['normal'] += 2*P.l*P.n*(P.β**2)*(α**2)
    if P.ε in dists['uniform']:
        dists['uniform'][P.ε] += 1+P.n
    else:
        dists['uniform'][P.ε] = 1+P.n

def blindrotate(domainP,targetP,dists):
    for i in range(domainP.n):
        cmux(targetP,targetP.α,dists)

def identitiykeyswithing(P,dists):
    dists['normal'] += P.t*P.domainP.n*(P.targetP.α**2)
    a = 2**(-P.basebit*P.t-1)
    if a in dists['uniform']:
        dists["uniform"][a] += P.domainP.n
    else:
        dists["uniform"][a] = P.domainP.n

def GateBootstrapping(brP,ikP,dists):
    blindrotate(brP.domainP,brP.targetP,dists)
    identitiykeyswithing(ikP,dists)

GateBootstrapping(lvl01param,lvl10param,dists)

print(dists)
print(dists["normal"]+((2*lvl1param.ε)**2)*dists["uniform"][lvl1param.ε]/12+((2*2**(-lvl10param.basebit*lvl10param.t-1))**2)*dists["uniform"][2**(-lvl10param.basebit*lvl10param.t-1)]/12)

import math
numccfunc, diffccfunc = ccfunc(1/16,dists)

from scipy.optimize import minimize

result = minimize(fun = numccfunc,x0 = [1e3], jac = diffccfunc, method = 'Nelder-Mead')
print(result['x'])
print(result['fun'])
print(2*math.exp(result['fun']))