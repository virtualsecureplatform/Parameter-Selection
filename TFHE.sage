load("lwe-estimator/estimator.py")

n=500; alpha=2.44e-5; q = 2^32;
print("lvl0")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1))
print("lvl0,quantum")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1), reduction_cost_model=BKZ.qsieve)

n=1024; alpha = 3.73e-9; q = 2^32;
print("lvl1")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1))
print("lvl1,quantum")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1), reduction_cost_model=BKZ.qsieve)

n=1024; alpha = 2^(-31); q = 2^32;
print("lvl21KSK")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1))
print("lvl21KSK,quantum")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1), reduction_cost_model=BKZ.qsieve)

n=2048; alpha = 2^(-44); q = 2^64;
print("lvl2")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1))
print("lvl2,quantum")
costs = estimate_lwe(n, sqrt(2*pi)*alpha, q, secret_distribution=(0,1), reduction_cost_model=BKZ.qsieve)