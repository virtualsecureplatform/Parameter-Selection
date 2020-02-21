load("lwe-estimator/estimator.py")

n=500; p = 2^13; q = 2^32;
alpha=sqrt(((q/p)^2-1)/12)/q;
print(alpha)
print("lvl0")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1))
print("lvl0,quantum")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1), reduction_cost_model=BKZ.qsieve)

n=1024; p = 2^29; q = 2^32;
alpha=sqrt(((q/p)^2-1)/12)/q;
print("lvl1")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1))
print("lvl1,quantum")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1), reduction_cost_model=BKZ.qsieve)

n=2048; p = 2^61; q = 2^64;
alpha=sqrt(((q/p)^2-1)/12)/q;
print("lvl2")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1))
print("lvl2,quantum")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1), reduction_cost_model=BKZ.qsieve)