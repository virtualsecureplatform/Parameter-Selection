import  numpy as np

# Helpers for Double Decomposition (TFHEpp "DD"/bivariate) parameters.
# If a parameter set does not define DD fields, we treat it as "no DD".
def _dd_levels(P, nonce: bool) -> int:
    return int(getattr(P, "l̅ₐ" if nonce else "l̅", 1))

def _dd_basebit(P, nonce: bool) -> int:
    return int(getattr(P, "B̅gₐbit" if nonce else "B̅gbit", 0))

def _qbit_from_q(q: int) -> int:
    # q is expected to be a power of two in all parameter sets used here.
    return int(q).bit_length() - 1

def _decomp_round_variance_pow2(q: int, basebit: int, levels: int, lbar: int, bbarbit: int) -> float:
    """
    Rounding variance for approximate gadget decomposition of a torus integer of width log2(q),
    when keeping `levels*basebit + (lbar-1)*bbarbit` most significant bits.

    This matches TFHEpp's DD constraint comment (e.g. 128-bit lvl3param uses:
      levels*basebit + (lbar-1)*bbarbit = 128).
    For lbar==1, the DD term disappears.
    """
    qbit = _qbit_from_q(q)
    covered_bits = int(levels) * int(basebit) + max(0, int(lbar) - 1) * int(bbarbit)
    remaining_bits = qbit - covered_bits
    if remaining_bits <= 0:
        return 0.0
    roundwidth = float(2 ** remaining_bits)
    return roundwidth * roundwidth / 12.0 - 1.0 / 12.0

# https://eprint.iacr.org/2021/729
def extpnoisecalc(P,α,β,exp,var):
    lbar = _dd_levels(P, nonce=False)
    lbara = _dd_levels(P, nonce=True)
    # Step 1
    res1 = ((P.lₐ * lbara) * P.k * P.n * (P.ℬₐ**2 + 2.) / 12. +
            (P.l * lbar) * P.n * (P.ℬ**2 + 2.) / 12.) * α
    # Step 2

    nonce_roundvar = _decomp_round_variance_pow2(
        P.q, int(P.ℬₐbit), int(P.lₐ), lbara, _dd_basebit(P, nonce=True)
    )
    nonnonce_roundvar = _decomp_round_variance_pow2(
        P.q, int(P.ℬbit), int(P.l), lbar, _dd_basebit(P, nonce=False)
    )
    noncevar = nonce_roundvar * (
        P.k * P.n * (P.variance_key_coefficient + P.expectation_key_coefficient**2)
    ) + P.k * P.n/4 * P.variance_key_coefficient
    nonnoncevar = nonnonce_roundvar
    trgswvar = β+noncevar+nonnoncevar
    squareexp = 1 / 4. * (1. - P.k * P.n * P.expectation_key_coefficient)**2
    # Last Part seems to be integer representation specific.
    return res1 + (var+exp**2)*trgswvar + var*squareexp

# https://eprint.iacr.org/2025/809
def brnoisecalc(lowP,highP = None, σ = 0):
    if highP is None:
        highP = lowP.targetP
        lowP = lowP.domainP
    acc = σ
    for i in range(lowP.k*lowP.n):
        acc = cmuxnoisecalc(highP,highP.σ,acc,acc,lowP.expectation_key_coefficient,lowP.variance_key_coefficient)
    return np.float64(acc)

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

# In CMUX, max(β,γ) must be retain since one of the input is selected.
def cmuxnoisecalc(P,α,β,γ,exp,var):
    # must be binominal distribution
    assert(abs(exp)<=1)
    assert(abs(var)<=1/2.)
    lbar = _dd_levels(P, nonce=False)
    lbara = _dd_levels(P, nonce=True)
    # Step 1
    res1 = ((P.lₐ * lbara) * P.k * P.n * (P.ℬₐ**2 + 2.) / 12. +
            (P.l * lbar) * P.n * (P.ℬ**2 + 2.) / 12.) * α
    # Step 2

    nonce_roundvar = _decomp_round_variance_pow2(
        P.q, int(P.ℬₐbit), int(P.lₐ), lbara, _dd_basebit(P, nonce=True)
    )
    nonnonce_roundvar = _decomp_round_variance_pow2(
        P.q, int(P.ℬbit), int(P.l), lbar, _dd_basebit(P, nonce=False)
    )
    noncevar = nonce_roundvar * (P.k * P.n * (P.variance_key_coefficient + P.expectation_key_coefficient**2)) + P.k * P.n/4 * P.variance_key_coefficient
    nonnoncevar = nonnonce_roundvar
    trgswvar = noncevar+nonnoncevar
    squareexp = 1 / 4. * (1. - P.k * P.n * P.expectation_key_coefficient)**2
    # Last Part seems to be integer representation specific.
    return res1 + (var+exp**2)*trgswvar + var*squareexp + max(β,γ)

    

def brroundnoise(domainP,targetP=None,ϑ=0):
    if targetP == None:
        targetP = domainP.targetP
        domainP = domainP.domainP
    roundwidth = domainP.q/(4*(targetP.n/2**ϑ))
    round_variance = (2*roundwidth)**2/12 - 1/12
    round_expectation = -1./2
    # return domainP.k*domainP.n*(round_variance*domainP.variance_key_coefficient+round_variance*domainP.expectation_key_coefficient**2+round_expectation**2 * domainP.variance_key_coefficient)
    return round_variance*(domainP.k*domainP.n*domainP.variance_key_coefficient+1)

def cbnoisecalc(brP,privksP):
    assert(brP.targetP == privksP.domainP)
    return brnoisecalc(brP)*((privksP.targetP.q/brP.targetP.q)**2)+privksnoisecalc(privksP)

# https://eprint.iacr.org/2021/729
def autoextpnoiseccalc(P):
    # Step 1
    res1 = (P.lₐ * P.k * P.n * (P.ℬₐ**2 + 2.) / 12) * P.σ
    # Step 2

    noncevar = (P.q**2-P.ℬₐ**(2*P.lₐ)) / (12 * P.ℬₐ**(2*P.lₐ)) * (P.k * P.n * (P.variance_key_coefficient + P.expectation_key_coefficient**2)) + P.k * P.n/4 * P.variance_key_coefficient
    trgswvar = noncevar
    squareexp = 1 / 4. * (- P.k * P.n * P.expectation_key_coefficient)**2
    # Last Part seems to be integer representation specific.
    return res1 + P.n*(P.variance_key_coefficient+P.expectation_key_coefficient**2)*trgswvar + P.n*P.variance_key_coefficient*squareexp

def annihilaterecursive(P,nbit,α):
    if(nbit==0):
        return α
    else:
        # return extpnoisecalc(P,P.σ,annihilaterecursive(P,nbit-1,α),P.expectation_key_coefficient,P.variance_key_coefficient) # Our algorithm divides the ciphertext in each loop
        return annihilaterecursive(P,nbit-1,α) + autoextpnoiseccalc(P)# Our algorithm divides the ciphertext in each loop

def annihilatecalc(P,α):
    return annihilaterecursive(P,P.nbit,α)

# https://eprint.iacr.org/2024/1318
def annihilatecbnoisecalc(domainP,targetP,ahP=None):
    if ahP == None:
        ahP = targetP
        targetP = domainP.targetP
        domainP = domainP.domainP
    # return extpnoisecalc(ahP,targetP.σ,annihilatecalc(ahP,brnoisecalc(domainP,targetP)),targetP.n*targetP.expectation_key_coefficient,targetP.n*targetP.variance_key_coefficient)
    # the brnoise will be annihilated except for the constant term
    return (targetP.variance_key_coefficient+targetP.expectation_key_coefficient**2)*brnoisecalc(domainP,targetP)+extpnoisecalc(ahP,targetP.σ,annihilatecalc(ahP,0),np.sqrt(targetP.n)*targetP.expectation_key_coefficient,targetP.n*targetP.variance_key_coefficient)

def romnoisecalc(brP,privksP,ROMaddress):
    # return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))+iksnoisecalc(addressP,dataP,ikP)
    σ = privksP.targetP.α**2
    for i in range(ROMaddress):
        σ = cmuxnoisecalc(privksP.targetP,cbnoisecalc(brP,privksP),σ,σ,1,0)
    return σ

def annihilateromnoisecalc(brP,ahP,ROMaddress):
    # return dataP.α**2+ROMaddress*cmuxnoisecalc(dataP,np.sqrt(cbnoisecalc(addressP,middleP,dataP,privksP)))+iksnoisecalc(addressP,dataP,ikP)
    σ = brP.targetP.α**2
    for i in range(ROMaddress):
        σ = cmuxnoisecalc(brP.targetP,annihilatecbnoisecalc(brP.domainP,brP.targetP,ahP),σ,σ,1,0)
    return σ

def ramnoisecalc(brP,iksP,cbbrP,privksP,RAMaddress,σ=None):
    assert(cbbrP.targetP == privksP.domainP)
    assert(brP.targetP == privksP.targetP)
    assert(iksP.targetP == brP.domainP)
    if σ == None:
        σ  = brnoisecalc(brP)
    for i in range(RAMaddress):
        σ = cmuxnoisecalc(privksP.targetP,cbnoisecalc(cbbrP,privksP),σ,σ,1,0)
    return σ+iksnoisecalc(iksP)

def annihilateramnoisecalc(brP,iksP,ahP,RAMaddress,σ=None):
    assert(iksP.targetP == brP.domainP)
    if σ == None:
        σ  = brnoisecalc(brP)
    for i in range(RAMaddress):
        σ = cmuxnoisecalc(brP.targetP,annihilatecbnoisecalc(brP,ahP),σ,σ,1,0)
    return σ+iksnoisecalc(iksP)

from scipy.special import erfc
def erfccalc(P,intereval,σ):
    return erfc((P.q/intereval)/np.sqrt(2*σ))
