#!/usr/bin/env python3
"""Security estimate for TFHEpp lvl5bootparam (and lvl5param baseline)."""
import sys, os, importlib

sys.path.insert(0, "/work")
est = importlib.import_module(".estimator", "lattice-estimator")
ND = est.nd
LWE = est.LWE
LWEParameters = est.lwe_parameters.LWEParameters

cases = []

# lvl5bootparam: n=2^14, Q=2^576, sparse ternary h=64, sigma = 2^(576-553) = 2^23
cases.append(LWEParameters(
    n=16384,
    q=2**576,
    Xs=ND.SparseTernary(32, 32, 16384),
    Xe=ND.DiscreteGaussian(stddev=2**23),
    m=16384,
    tag="lvl5boot(h64,q576,s23)",
))

# lvl5param baseline: dense ternary, Q=2^448, sigma=2^23
cases.append(LWEParameters(
    n=16384,
    q=2**448,
    Xs=ND.UniformMod(3),
    Xe=ND.DiscreteGaussian(stddev=2**23),
    m=16384,
    tag="lvl5(dense,q448,s23)",
))

for p in cases:
    print("=" * 70)
    print(p.tag, flush=True)
    try:
        r = LWE.estimate.rough(p)
        for k, v in r.items():
            print(f"  {k}: rop=2^{float(v['rop']).__format__('.1f')}" if hasattr(v, '__getitem__') else f"  {k}: {v}")
    except Exception as e:
        print("  rough estimate failed:", repr(e))
    print(flush=True)
