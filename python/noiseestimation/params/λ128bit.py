class lvl0param:
  n = 636
  k = 1
  q = 2**32
  α = 0.000_092_511_997_467_675_6 * q
  σ = α**2
  # binary
  variance_key_coefficient = 1./4
  expectation_key_coefficient = 1./2

class lvl1param:
    nbit = 10
    q = 2**32
    n = 2**nbit
    k = 1
    l = 2
    lₐ = l
    ℬbit = 8
    ℬₐbit = ℬbit
    ℬ = 2**ℬbit
    ℬₐ = 2**ℬₐbit
    α = 0.0000000342338787018369 * q
    σ = α**2
    # ternary
    variance_key_coefficient = 2./3
    expectation_key_coefficient = 0

class lvl2param:
    nbit = 11
    k = 1   
    n = 2**nbit
    q = 2**64
    l = 4
    lₐ = 4
    ℬbit = 9
    ℬₐbit = 9
    ℬ = 2**ℬbit
    ℬₐ = 2**ℬₐbit
    # Double Decomposition (TFHEpp DD / bivariate) parameters (trivial here)
    l̅ = 1
    l̅ₐ = 1
    B̅gbit = 0
    B̅gₐbit = 0
    α =  q * 2**-51
    σ = α**2
    variance_key_coefficient = 2./3
    expectation_key_coefficient = 0.

class lvl3param:
    # Proposed lvl3 parameter set (128-bit torus + DD/bivariate representation).
    nbit = 12
    k = 1
    n = 2**nbit
    q = 2**128
    # Primary (TRGSW) decomposition
    # Chosen so that approximation noise is meaningfully below lvl02:
    # keep_bits = l*ℬbit = 4*19 = 76 (std/q ≈ 2^-68 in this estimator).
    l = 4
    lₐ = 4
    ℬbit = 19
    ℬₐbit = 19
    ℬ = 2**ℬbit
    ℬₐ = 2**ℬₐbit
    # Double Decomposition (auxiliary) parameters
    # Conservative full-limb cover: l̅ * B̅gbit = 128 (and same for nonce part).
    # Chosen to satisfy FFT safety constraint: ℬbit + B̅gbit + nbit + 3 < 53.
    l̅ = 8
    l̅ₐ = 8
    B̅gbit = 16
    B̅gₐbit = 16
    # Fresh noise (TFHEpp uses α normalized; python stores α*q and σ=α^2)
    α = q * 2**-105
    σ = α**2
    # Secret distribution in TFHEpp is uniform over {-1, 0, +1} (key_value_min=-1, max=+1).
    variance_key_coefficient = 2.0/3.0
    expectation_key_coefficient = 0.0

class annihilatelvl2param:
    nbit = lvl2param.nbit
    k = lvl2param.k
    n = lvl2param.n
    q = lvl2param.q
    l = 5
    lₐ = 5
    ℬbit = 8
    ℬₐbit = 8
    ℬ = 2**ℬbit
    ℬₐ = 2**ℬₐbit
    α = lvl2param.α
    σ = lvl2param.σ
    variance_key_coefficient = lvl2param.variance_key_coefficient
    expectation_key_coefficient = lvl2param.expectation_key_coefficient

class lvl10param:
    t = 4
    basebit = 3
    domainP = lvl1param
    targetP = lvl0param

class lvl11param:
    t = 5
    basebit = 5
    domainP = lvl2param
    targetP = lvl2param

class lvl21param:
    t = 9
    basebit = 3
    domainP = lvl2param
    targetP = lvl1param

class lvl22param:
    t = 38
    basebit = 1
    domainP = lvl2param
    targetP = lvl2param

class lvl20param:
    t = 7
    basebit =  2
    domainP = lvl2param
    targetP = lvl0param

class lvl21mrlweparam:
    nbit = lvl1param.nbit
    n = 2**nbit
    k = lvl1param.k
    q = 2**32
    l = 4
    ℬbit = 6
    ℬ = 2**ℬbit
    α = lvl1param.α
    # ternary
    variance_key_coefficient = lvl1param.variance_key_coefficient
    expectation_key_coefficient = lvl1param.expectation_key_coefficient

class lvl01param:
    domainP = lvl0param
    targetP = lvl1param
class lvl02param:
    domainP = lvl0param
    targetP = lvl2param

class lvl03param:
    domainP = lvl0param
    targetP = lvl3param
