#!/bin/python3
"""
Extra checks:
  n=636:  q=2^15, binary-search for min eta
  n=1024: q=2^25,2^26,2^27 with eta=4 only
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

# --- n=636, q=2^15: binary search for min eta ---
print("=" * 70, flush=True)
print("n=636, q=2^15: binary search for min eta >= 128-bit", flush=True)
print("=" * 70, flush=True)
n, Xs, log2q = 636, estimator.nd.UniformMod(2), 15
lo, hi = 1, 64

sec_hi = estimate_security(n, Xs, log2q, hi)
print(f"  eta={hi}: {sec_hi:.1f} bit", flush=True)
if sec_hi < TARGET:
    print(f"  NOT secure even at eta={hi}", flush=True)
else:
    sec_lo = estimate_security(n, Xs, log2q, lo)
    print(f"  eta={lo}: {sec_lo:.1f} bit", flush=True)
    if sec_lo >= TARGET:
        print(f"  => min eta={lo} ({sec_lo:.1f} bit)", flush=True)
    else:
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            sec_mid = estimate_security(n, Xs, log2q, mid)
            print(f"  eta={mid}: {sec_mid:.1f} bit", flush=True)
            if sec_mid >= TARGET:
                hi = mid
            else:
                lo = mid
        sec_final = estimate_security(n, Xs, log2q, hi)
        print(f"  => min eta={hi} ({sec_final:.1f} bit)", flush=True)

# --- n=1024, q=2^25..2^27 with eta=4 ---
print(flush=True)
print("=" * 70, flush=True)
print("n=1024, q=2^25..2^27, eta=4", flush=True)
print("=" * 70, flush=True)
Xs1024 = estimator.nd.UniformMod(3)
for log2q in [25, 26, 27]:
    sec = estimate_security(1024, Xs1024, log2q, 4)
    status = "PASS" if sec >= TARGET else "FAIL"
    print(f"  q=2^{log2q}, eta=4: {sec:.1f} bit  [{status}]", flush=True)

print("\nDone.", flush=True)
