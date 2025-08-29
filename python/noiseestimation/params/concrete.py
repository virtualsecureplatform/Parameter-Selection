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
    nbit = 9
    q = 2**32
    # q = 5**4*2**16+1
    # q = 3**4*2**16+1
    # q = 2**25
    n = 2**nbit
    k = 2
    l = 1
    lₐ = 2
    ℬbit = 10
    ℬₐbit = 8
    ℬ = 2**ℬbit
    ℬₐ = 2**ℬₐbit
    α = 0.0000000342338787018369 * q
    σ = α**2
    # α = 2*2/4
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
    α =  q * 2**-51
    σ = α**2
    ε = 1/(2*(ℬ**l))
    β = ℬ/2
    variance_key_coefficient = 2./3
    expectation_key_coefficient = 0.

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
