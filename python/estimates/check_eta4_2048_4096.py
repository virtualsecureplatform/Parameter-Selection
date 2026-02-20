#!/bin/python3
"""
Check CBD(4) (eta=4) security for n=2048 and n=4096 at various q.
"""
import sys, os
from math import log2
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

def min_rop(result):
    best = None
    for attack, res in result.items():
        rop_val = res.get("rop", 2**1000)
        try:
            rop = float(rop_val.log(2))
        except AttributeError:
            rop = log2(float(rop_val))
        if best is None or rop < best:
            best = rop
    return best

def estimate_security(n, Xs, log2q, eta):
    param = estimator.lwe_parameters.LWEParameters(
        n=n,
        q=2 ** log2q,
        Xs=Xs,
        Xe=estimator.nd.CenteredBinomial(eta),
        tag=f"n{n}_q2^{log2q}_cbd{eta}",
    )
    r = estimator.LWE.estimate(param, red_cost_model=estimator.RC.BDGL16)
    sec = min_rop(r)
    return sec

TARGET = 128
Xs = estimator.nd.UniformMod(3)

# n=2048, eta=4: original q=2^64, test q=2^32..2^56
print("=" * 70, flush=True)
print("n=2048, eta=4, Xs=UniformMod(3)", flush=True)
print("=" * 70, flush=True)
for log2q in [32, 36, 40, 44, 48, 52, 56]:
    sec = estimate_security(2048, Xs, log2q, 4)
    status = "PASS" if sec >= TARGET else "FAIL"
    print(f"  q=2^{log2q}, eta=4: {sec:.1f} bit  [{status}]", flush=True)

# n=4096, eta=4: original q=2^128, test q=2^64..2^112
print(flush=True)
print("=" * 70, flush=True)
print("n=4096, eta=4, Xs=UniformMod(3)", flush=True)
print("=" * 70, flush=True)
for log2q in [64, 72, 80, 88, 96, 104, 112]:
    sec = estimate_security(4096, Xs, log2q, 4)
    status = "PASS" if sec >= TARGET else "FAIL"
    print(f"  q=2^{log2q}, eta=4: {sec:.1f} bit  [{status}]", flush=True)

print("\nDone.", flush=True)
