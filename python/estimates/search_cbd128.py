#!/bin/python3
"""
Search for 128-bit secure CBD parameters for n=636 and n=1024.
Strategy: for each q (power of 2), binary-search for minimum eta that gives >=128-bit security.
"""
import sys, os
from math import log2
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import importlib
estimator = importlib.import_module(".estimator","lattice-estimator")

def min_rop(result):
    """Extract minimum rop (log2) from estimation result."""
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

# n=636: q around 2^16 (original TFHE636-16bit used q=2^16, CBD(40))
# n=1024: q around 2^32 (original TFHE1024 used q=2^32, Gaussian stddev=2^7)
configs = [
    (636,  estimator.nd.Binary, [12, 14, 16, 18, 20]),
    (1024, estimator.nd.UniformMod(3), [16, 20, 24, 28, 32]),
]

for n, Xs, q_bits_list in configs:
    print("=" * 70, flush=True)
    print(f"Searching: n={n}, target >= {TARGET}-bit security", flush=True)
    print("=" * 70, flush=True)

    results = []
    for log2q in q_bits_list:
        # Binary search for minimum eta that gives >= 128-bit security at this q.
        # Larger eta => more noise => more security.
        lo, hi = 1, 64

        # First check: does eta=hi give 128-bit security?
        sec_hi = estimate_security(n, Xs, log2q, hi)
        print(f"  q=2^{log2q}, eta={hi}: {sec_hi:.1f} bit", flush=True)
        if sec_hi < TARGET:
            print(f"  q=2^{log2q}: NOT secure even at eta={hi}, skipping", flush=True)
            results.append((log2q, None, sec_hi))
            continue

        # Check: does eta=lo already pass?
        sec_lo = estimate_security(n, Xs, log2q, lo)
        print(f"  q=2^{log2q}, eta={lo}: {sec_lo:.1f} bit", flush=True)
        if sec_lo >= TARGET:
            results.append((log2q, lo, sec_lo))
            print(f"  => q=2^{log2q}: min eta={lo} ({sec_lo:.1f} bit)", flush=True)
            print(flush=True)
            continue

        # Binary search: find minimum eta where security >= TARGET
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            sec_mid = estimate_security(n, Xs, log2q, mid)
            print(f"  q=2^{log2q}, eta={mid}: {sec_mid:.1f} bit", flush=True)
            if sec_mid >= TARGET:
                hi = mid
            else:
                lo = mid

        sec_hi = estimate_security(n, Xs, log2q, hi)
        results.append((log2q, hi, sec_hi))
        print(f"  => q=2^{log2q}: min eta={hi} ({sec_hi:.1f} bit)", flush=True)
        print(flush=True)

    print(f"\n--- Summary for n={n} ---", flush=True)
    print(f"{'q':>8s}  {'min eta':>8s}  {'security':>10s}", flush=True)
    for log2q, eta, sec in results:
        eta_str = str(eta) if eta else "N/A"
        print(f"{'2^'+str(log2q):>8s}  {eta_str:>8s}  {sec:10.1f}", flush=True)
    print(flush=True)
