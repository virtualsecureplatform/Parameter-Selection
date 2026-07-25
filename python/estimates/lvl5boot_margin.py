#!/usr/bin/env python3
"""Targeted hybrid-attack estimates to pick the security margin for
lvl5bootparam: vary sigma and h at n=2^15, q=2^576."""
import sys, importlib

sys.path.insert(0, "/work")
est = importlib.import_module(".estimator", "lattice-estimator")
ND = est.nd
LWE = est.LWE
LWEParameters = est.lwe_parameters.LWEParameters

variants = [
    (33, 64),
    (33, 96),
    (23, 96),
]

for sbits, h in variants:
    p = LWEParameters(
        n=32768, q=2**576,
        Xs=ND.SparseTernary(h // 2, h - h // 2, 32768),
        Xe=ND.DiscreteGaussian(stddev=2**sbits),
        m=32768, tag=f"s{sbits} h{h}",
    )
    print("=" * 60)
    print(p.tag, flush=True)
    for attack in ("bdd_mitm_hybrid", "bdd_hybrid"):
        try:
            f = getattr(est.lwe_primal.primal_hybrid, "__call__")
            if attack == "bdd_mitm_hybrid":
                r = est.LWE.primal_hybrid(p, mitm=True, babai=True,
                                          red_cost_model=est.RC.BDGL16)
            else:
                r = est.LWE.primal_hybrid(p, mitm=False, babai=False,
                                          red_cost_model=est.RC.BDGL16)
            print(f"  {attack}: {r!r}", flush=True)
        except Exception as e:
            print(f"  {attack} failed: {e!r}", flush=True)
