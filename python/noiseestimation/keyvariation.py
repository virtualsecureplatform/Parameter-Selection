import  numpy as np

# https://eprint.iacr.org/2021/729
def extpnoisecalc(P,α,β,exp,var):
    # Step 1
    res1 = (P.lₐ * P.k * P.n * (P.ℬₐ**2 + 2.) / 12. + P.l * P.n * (P.ℬ**2 + 2.) / 12.) * α
    # Step 2

    noncevar = (P.q**2-P.ℬₐ**(2*P.lₐ)) / (12 * P.ℬₐ**(2*P.lₐ)) * (P.k * P.n * (P.variance_key_coefficient + P.expectation_key_coefficient**2)) + P.k * P.n/4 * P.variance_key_coefficient
    nonnoncevar = (P.q**2-P.ℬ**(2*P.l)) / (12 * P.ℬ**(2*P.l))
    trgswvar = β+noncevar+nonnoncevar
    squareexp = 1 / 4. * (1. - P.k * P.n * P.expectation_key_coefficient)**2
    # Last Part seems to be integer representation specific.
    return res1 + (var+exp**2)*trgswvar + var*squareexp

def brnoisecalc(lowP,highP = None, σ = 0):
    if highP is None:
        highP = lowP.targetP
        lowP = lowP.domainP
    return lowP.k * lowP.n * extpnoisecalc(highP,highP.σ,σ,lowP.expectation_key_coefficient,lowP.variance_key_coefficient)

def unrollbrnoisecalc(lowP,highP,m):
    res1 = highP.l * (highP.k + 1.) * highP.n * (highP.ℬ**2 + 2.) / 12. * (2**m - 1)*(highP.α)**2
    res2 = (highP.q**2-highP.ℬ**(2*highP.l)) / (24 * highP.ℬ**(2*highP.l)) * (1. + highP.k * highP.n * (lowP.variance_key_coefficient + lowP.expectation_key_coefficient**2)) + highP.k * highP.n/8 * lowP.variance_key_coefficient  + 1 / 16. * (1. - highP.k * highP.n * lowP.expectation_key_coefficient)**2; # Last Part seems to be integer representation specific.
    return lowP.k * lowP.n/m * (res1+res2)

def iksnoisecalc(lowP,highP=None,funcP=None):
    if highP is None:
        funcP = lowP
        highP = lowP.domainP
        lowP = lowP.targetP
    # return 1/12*highP.k*highP.n*(2**(-2*(funcP.basebit*funcP.t)))+funcP.t*highP.k*highP.n*(lowP.α**2)
    roundwidth = 2**(-funcP.basebit*funcP.t-1) * lowP.q
    round_variance = (lowP.q**2-2**(funcP.basebit*funcP.t+1))/(12*2**(2*funcP.basebit*funcP.t+1))
    round_expectation = -1./2
    return highP.k*highP.n*((round_variance*highP.variance_key_coefficient+round_variance*highP.expectation_key_coefficient**2+round_expectation**2 * highP.variance_key_coefficient)+funcP.t*(lowP.α**2))

def privksnoisecalc(lowP,highP=None,funcP=None):
    if highP is None:
        funcP = lowP
        highP = lowP.domainP
        lowP = lowP.targetP
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

def cmuxnoisecalc(P,α,β,γ,exp,var):
    return extpnoisecalc(P,α,max(β,γ),exp,var)

def brroundnoise(domainP,targetP,ϑ=0):
    if targetP == None:
        targetP = domainP.targetP
        domainP = domainP.domainP
    roundwidth = domainP.q/(4*(targetP.n/2**ϑ))
    round_variance = (2*roundwidth)**2/12 - 1/12
    round_expectation = -1./2
    return domainP.k*domainP.n*(round_variance*domainP.variance_key_coefficient+round_variance*domainP.expectation_key_coefficient**2+round_expectation**2 * domainP.variance_key_coefficient)
    # return domainP.k*domainP.n*1/4*targetP.n*domainP.q

def cbnoisecalc(brP,privksP):
    assert(brP.targetP == privksP.domainP)
    return brnoisecalc(brP)*((privksP.targetP.q/brP.targetP.q)**2)+privksnoisecalc(privksP)

# https://eprint.iacr.org/2021/729
def autoextpnoiseccalc(P):
    # Step 1
    res1 = (P.lₐ * P.k * P.n * (P.ℬₐ**2 + 2.) / 12) * P.σ
    # res2 = P.ℬ**2/2
    # Step 2

    noncevar = (P.q**2-P.ℬₐ**(2*P.lₐ)) / (12 * P.ℬₐ**(2*P.lₐ)) * (P.k * P.n * (P.variance_key_coefficient + P.expectation_key_coefficient**2)) + P.k * P.n/4 * P.variance_key_coefficient
    trgswvar = noncevar
    squareexp = 1 / 4. * (- P.k * P.n * P.expectation_key_coefficient)**2
    # Last Part seems to be integer representation specific.
    return res1 + (P.variance_key_coefficient+P.expectation_key_coefficient**2)*trgswvar + P.variance_key_coefficient*squareexp

def annihilaterecursive(P,nbit,α):
    if(nbit==0):
        return α
    else:
        # return extpnoisecalc(P,P.σ,annihilaterecursive(P,nbit-1,α),P.expectation_key_coefficient,P.variance_key_coefficient) # Our algorithm divides the ciphertext in each loop
        return annihilaterecursive(P,nbit-1,α) + autoextpnoiseccalc(P)# Our algorithm divides the ciphertext in each loop

def annihilatecalc(P,α):
    return annihilaterecursive(P,P.nbit,α)

def annihilatecbnoisecalc(domainP,targetP):
    return annihilatecalc(targetP,brnoisecalc(domainP,targetP))

def romnoisecalc(brP,privksP,ROMaddress):
    # return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))+iksnoisecalc(addressP,dataP,ikP)
    σ = privksP.targetP.α**2
    for i in range(ROMaddress):
        σ = cmuxnoisecalc(privksP.targetP,cbnoisecalc(brP,privksP),σ,σ,1,0)
    return σ

def annihilateromnoisecalc(brP,ROMaddress):
    # return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))+iksnoisecalc(addressP,dataP,ikP)
    σ = brP.targetP.α**2
    for i in range(ROMaddress):
        σ = cmuxnoisecalc(brP.targetP,annihilatecbnoisecalc(brP.domainP,brP.targetP),σ,σ,1,0)
    return σ

def ramnoisecalc(brP,iksP,cbbrP,privksP,RAMaddress):
    assert(cbbrP.targetP == privksP.domainP)
    assert(brP.targetP == privksP.targetP)
    assert(iksP.targetP == brP.domainP)
    σ  = brnoisecalc(brP)
    for i in range(RAMaddress):
        σ = cmuxnoisecalc(privksP.targetP,cbnoisecalc(cbbrP,privksP),σ,σ,1,0)
    return σ+iksnoisecalc(iksP)