#!/usr/bin/env python3
"""Security sweep for the BFV bootstrap parameter redesign (n=2^15 ring)."""
import sys, importlib

sys.path.insert(0, "/work")
est = importlib.import_module(".estimator", "lattice-estimator")
ND = est.nd
LWE = est.LWE
LWEParameters = est.lwe_parameters.LWEParameters


def case(n, qbits, sbits, secret, tag):
    if secret == "dense":
        Xs = ND.UniformMod(3)
    else:
        h = int(secret)
        Xs = ND.SparseTernary(h // 2, h - h // 2, n)
    return LWEParameters(
        n=n, q=2**qbits, Xs=Xs,
        Xe=ND.DiscreteGaussian(stddev=2**sbits),
        m=n, tag=tag,
    )


grid = [
    case(32768, 576, 23, "64",    "n32768 q576 s23 h64"),
    case(32768, 576, 23, "128",   "n32768 q576 s23 h128"),
    case(32768, 576, 23, "192",   "n32768 q576 s23 h192"),
    case(32768, 576, 23, "dense", "n32768 q576 s23 dense"),
    case(32768, 576, 33, "64",    "n32768 q576 s33 h64"),
    case(32768, 576, 33, "128",   "n32768 q576 s33 h128"),
    case(32768, 640, 33, "128",   "n32768 q640 s33 h128"),
    case(32768, 512, 23, "64",    "n32768 q512 s23 h64"),
]

for p in grid:
    print("=" * 70)
    print(p.tag, flush=True)
    try:
        LWE.estimate.rough(p)
    except Exception as e:
        print("  rough estimate failed:", repr(e))
    print(flush=True)
