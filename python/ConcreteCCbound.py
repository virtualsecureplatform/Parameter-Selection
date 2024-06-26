#!/bin/python3
from sympy import integrate,Symbol,pi,oo,log,diff
from sympy.functions import exp, sqrt
from sympy.printing.julia import JuliaCodePrinter
from sympy.printing.octave import print_octave_code

import numpy as np

# 128bit TFHE's parameter.

class lvl0param:
    n = 636
    α = 0.0000925119974676756
# class lvl0param:
#     n = 512
#     α = 2**-12

# class lvl1param:
#     nbit = 8
#     k = 4
#     n = 2**nbit
#     l = 2
#     Bgbit = 8
#     Bg = 2**Bgbit
#     α = 0.0000000342338787018369
#     ε = 1/(2*(Bg**l))
#     β = Bg/2

# class lvl1param:
#     nbit = 9
#     k = 2
#     n = 2**nbit
#     l = 2
#     Bgbit = 8
#     Bg = 2**Bgbit
#     α = 0.0000000342338787018369
#     ε = 1/(2*(Bg**l))
#     β = Bg/2

class lvl1param:
    nbit = 10
    k = 1
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
    k = 1
    l = 4
    Bgbit = 9
    Bg = 2**Bgbit
    α = 2**-44
    ε = 1/(2*(Bg**l))
    β = Bg/2

# class lvl2param:
#    nbit = 9
#    n = 2**nbit
#    k = 3
#    l = 3
#    Bgbit = 9
#    Bg = 2**Bgbit
#    α = 0.0000000000034525330484572114
#    ε = 1/(2*(Bg**l))
#    β = Bg/2

# class lvl10param:
#     t = 3
#     basebit = 4
#     domainP = lvl1param
#     targetP = lvl0param

class lvl10param:
    t = 7
    basebit = 2
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
    t = 32
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


ROMaddress = 9 # 4 word block
RAMaddress = 9
RAMwordbit = 8

def ccfunc(μ,dists):
    x = Symbol('x')
    t = Symbol('t')
    func = -t*μ
    # func = exp(-t*μ)
    if "uniform" in dists:
        for interval, num in dists["uniform"].items():
            # func += num*log(integrate(exp(t*x)/(2*interval),(x,-interval,interval)))
            func += num * log((exp(t*interval)-exp(t*-interval))/(2*t*interval))
            # func *= integrate(exp(t*x)/(2*interval),(x,-interval,interval))**num
    if "normal" in dists:
        # https://people.math.wisc.edu/~roch/grad-prob/gradprob-notes7.pdf
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

from scipy.special import logsumexp
def precccfunc(ta,μ,dists):
    t= ta[0]
    funclist = [-t*μ]
    # func = math.exp(-t*μ)
    if "uniform" in dists:
        for interval, num in dists["uniform"].items():
            # funclist.append(num * (logsumexp([t*interval,t*-interval],b=[1,-1]) - math.log(2*t*interval)))
            funclist.append(logsumexp([num*np.log(np.expm1(2*t*interval)/(2*interval*t)),math.log(num*interval*t)],b=[1,-1]))
            # funclist.append(math.log((np.expm1(2*t*interval)/(2*interval*t))**num-(num*interval*t)))
            # print([math.log(math.exp(t*interval)-math.exp(t*-interval)),logsumexp([t*interval,t*-interval],b=[1,-1]),math.log(2*t*interval),logsumexp([t*interval,t*-interval],b=[1,-1])-math.log(2*t*interval),logsumexp([t*interval- math.log(2*t*interval),t*-interval- math.log(2*t*interval)],b=[1,-1]),num * (logsumexp([t*interval - math.log(2*t*interval),t*-interval - math.log(2*t*interval)],b=[1,-1]))])
    if "normal" in dists:
        # https://people.math.wisc.edu/~roch/grad-prob/gradprob-notes7.pdf
        # func += log(integrate(exp(t*x)/(sqrt(2*pi*dists["normal"]))*exp(-(x**2)/(2*dists["normal"])),(x,-oo,oo)))
        # func *= math.exp((t**2)*dists["normal"]/2)
        funclist.append((t**2)*dists["normal"]/2)
    # return func
    return math.fsum(funclist)
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
    dists['normal'] += P.t*P.domainP.k*P.domainP.n*(P.targetP.α**2)
    a = 2**(-P.basebit*P.t-1)
    if a in dists['uniform']:
        dists["uniform"][a] += P.domainP.k*P.domainP.n
    else:
        dists["uniform"][a] = P.domainP.k*P.domainP.n

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

def BRround(brP,manybit,dists):
    roundwidth = 2**manybit/(4*brP.targetP.n)
    if roundwidth in dists["uniform"]:
        dists["uniform"][roundwidth] += brP.domainP.n
    else:
        dists["uniform"][roundwidth] = brP.domainP.n

def ChensPackingCircuitBootstrapping(brP,dists):
    blindrotate(brP,dists)
    annihilatekeyswitching(brP.targetP,dists)

def romnoisecalc(brP,ikP,privksP,ROMaddress):
    dists = {'normal': privksP.targetP.α**2,'uniform':{}}

    cbdists = {'normal': 0,'uniform':{}}
    CircuitBootstrapping(brP,privksP,cbdists)
    for i in range(ROMaddress):
        cmux(privksP.targetP,cbdists,dists)
    identitiykeyswithing(ikP,dists)
    return dists

def ChensPackingromnoisecalc(brP,ikP,privksP,ROMaddress):
    dists = {'normal': privksP.targetP.α**2,'uniform':{}}
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
# GateBootstrapping(lvl01param,lvl10param,dists)
# GateBootstrapping(lvl12param,lvl21param,dists)
# BRround(lvl01param,0,dists)
# privatekeyswitching(lvl11param,dists)
privatekeyswitching(lvl21param,dists)
# privatekeyswitching(lvl22param,dists)
# annihilatekeyswitching(lvl2param,dists)
# CircuitBootstrapping(lvl01param,lvl11param,dists)
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
print("conventional std")
print(np.sqrt(conventionalvariance))
print("conventional variance")
print(conventionalvariance)

messagedistance = 1/8
from scipy.special import erfc
print("conventional error rate")
print(2*erfc(messagedistance/(2*np.sqrt(2*conventionalvariance))))

import math
numccfunc, diffccfunc, funwithjac = ccfunc(1/4,dists)



from scipy.optimize import minimize,shgo,dual_annealing,direct

# result = minimize(fun = numccfunc,x0 = np.array([1e-6]), jac = diffccfunc, method = 'Newton-CG')
# print(result['x'])
# print(result['fun'])
# print(2*math.exp(result['fun']))

# result = shgo(numccfunc,bounds=[(1e-6,None)],minimizer_kwargs={'method': "SLSQP", 'jac':diffccfunc})
# result = dual_annealing(numccfunc,bounds=[(1e-3,1e5)],initial_temp=1e4)
# result = direct(numccfunc,bounds=[(1e-3,1e5)])

# result = direct(precccfunc,bounds=[(1e-8,1e4)],args=[messagedistance/2,dists])
result = dual_annealing(precccfunc,bounds=[(1e-6,1e4)],args=[messagedistance/2,dists])

print("Result")
print(precccfunc(result.x,messagedistance/2,dists))
print(result.x)
# print(1/(4*dists["normal"]))
print(result.fun)
# print(precccfunc([1/(4*dists["normal"])],messagedistance/2,dists))
# print(2*math.exp(-1/(2*4**2*dists["normal"])))
print(2*math.exp(result.fun))
