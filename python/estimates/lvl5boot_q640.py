#!/usr/bin/env python3
"""Hybrid estimate at q=2^640 for the 10-limb candidate (h=96, sigma=2^33)."""
import sys, importlib
sys.path.insert(0, "/work")
est = importlib.import_module(".estimator", "lattice-estimator")
ND = est.nd
LWEParameters = est.lwe_parameters.LWEParameters
p = LWEParameters(n=32768, q=2**640,
                  Xs=ND.SparseTernary(48, 48, 32768),
                  Xe=ND.DiscreteGaussian(stddev=2**33),
                  m=32768, tag="q640 s33 h96")
print(p, flush=True)
r = est.LWE.primal_hybrid(p, mitm=True, babai=True, red_cost_model=est.RC.BDGL16)
print("bdd_mitm_hybrid:", repr(r), flush=True)
r2 = est.LWE.primal_hybrid(p, mitm=False, babai=False, red_cost_model=est.RC.BDGL16)
print("bdd_hybrid:", repr(r2), flush=True)
