# To reproduce the estimate run this snippet on http://aleph.sagemath.org/
from sage.all import load, sqrt, RR, ZZ, pi, oo
load('https://bitbucket.org/malb/lwe-estimator/raw/HEAD/estimator.py')

n = 1024                # ciphertext dimension (also, key entropy)
sd = 2**(-29)            # noise standard deviation
alpha = sqrt(2*pi)*sd    # estimator defines noise rate = sqrt(2pi).stdev
q = 2**32                # for compatibility only
m = oo                   # the attacker can use as many samples he wishes 
secret_distribution = (0,1)
secret_distribution = uniform
success_probability = 0.99


# Chosen cost model 
# BKZ cost models: CLASSICAL - 0.292*beta + 16.4 + log(8*d,2) - primal
# i.e. BKZ.sieve =  lambda beta, d, B: ZZ(2)**RR(0.292*beta + 16.4 + log(8*d,2))
print("CLASSICAL PRIMAL")
print(primal_usvp(n, alpha, q, secret_distribution=secret_distribution, m=m, success_probability=success_probability, reduction_cost_model=BKZ.sieve))
# BKZ cost models: CLASSICAL - 0.292*beta + 16.4 + log(8*d,2) - dual
# i.e. BKZ.sieve =  lambda beta, d, B: ZZ(2)**RR(0.292*beta + 16.4 + log(8*d,2))
print("CLASSICAL DUAL")
print(dual_scale(n, alpha, q, secret_distribution=secret_distribution, m=m, success_probability=success_probability, reduction_cost_model=BKZ.sieve))


# For more conservative parameters, both classical and quantum  
# BKZ cost models: CLASSICAL - 0.292 beta - primal
reduction_cost_model =  lambda beta, d, B: ZZ(2)**RR(0.292*beta)
print("CLASSICAL PRIMAL (conservative)")
print(primal_usvp(n, alpha, q, secret_distribution=secret_distribution, m=m, success_probability=success_probability, reduction_cost_model=reduction_cost_model))
# BKZ cost models: CLASSICAL - 0.292 beta - dual
print("CLASSICAL DUAL (conservative)")
print(dual_scale(n, alpha, q, secret_distribution=secret_distribution, m=m, success_probability=success_probability, reduction_cost_model=reduction_cost_model))
# BKZ cost models: QUANTUM - 0.265 beta - primal
reduction_cost_model =  lambda beta, d, B: ZZ(2)**RR(0.265*beta)
print("QUANTUM PRIMAL (conservative)")
print(primal_usvp(n, alpha, q, secret_distribution=secret_distribution, m=m, success_probability=success_probability, reduction_cost_model=reduction_cost_model))
# BKZ cost models: QUANTUM - 0.265 beta - dual
print("QUANTUM DUAL (conservative)")
print(dual_scale(n, alpha, q, secret_distribution=secret_distribution, m=m, success_probability=success_probability, reduction_cost_model=reduction_cost_model))
