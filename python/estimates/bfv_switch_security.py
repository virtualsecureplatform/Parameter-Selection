#!/usr/bin/env sage -python
"""Reproduce the scalar TFHE/BFV chain's four-attack classical LWE screen.

This is an RLWE-as-LWE heuristic screen, not a circular/KDM security proof.
Run with SageMath and the repository's lattice-estimator checkout.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys

from sage.all import oo

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))
estimator = importlib.import_module(".estimator", "lattice-estimator")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--estimator-revision", help="revision of an exported checkout when git is unavailable")
    args = parser.parse_args()
    report = {
        "utc": datetime.now(timezone.utc).isoformat(),
        "cost_model": "RC.BDGL16 (classical)", "samples": "unbounded",
        "target_bits": 128,
        "scope": "Four-attack heuristic LWE screen; excludes ring structure and circular/KDM assumptions",
        "estimator_revision": args.estimator_revision or subprocess.check_output(
            ["git", "-C", str(PYTHON_ROOT / "lattice-estimator"), "rev-parse", "HEAD"], text=True).strip(),
        "cases": {},
    }
    digest = hashlib.sha256()
    for source in sorted((PYTHON_ROOT / "lattice-estimator/estimator").glob("*.py")):
        digest.update(source.name.encode())
        digest.update(source.read_bytes())
    report["estimator_sources_sha256"] = digest.hexdigest()
    header = PYTHON_ROOT.parents[1] / "TFHEpp/include/params/128bit.hpp"
    report["parameter_header_sha256"] = hashlib.sha256(header.read_bytes()).hexdigest()
    attacks = {"primal_usvp": estimator.LWE.primal_usvp,
               "primal_bdd": estimator.LWE.primal_bdd,
               "dual": estimator.LWE.dual,
               "dual_hybrid": estimator.LWE.dual_hybrid}
    for name, n, qbits, sigma, alphabet in (
        ("tfhe_io", 1024, 32, 2**7, 3),
        ("tfhe_half", 760, 32, 2**15, 2),
        ("bfv_ring", 4096, 128, 2**23, 3),
    ):
        params = estimator.lwe_parameters.LWEParameters(
            n=n, q=2**qbits, Xs=estimator.nd.Binary if alphabet == 2 else estimator.nd.Ternary,
            Xe=estimator.nd.DiscreteGaussian(stddev=sigma), m=oo, tag=name)
        case = dict(n=n, qbits=qbits, sigma=sigma, secret_alphabet=alphabet, attacks={})
        case["secret_bounds"] = [int(bound) for bound in params.Xs.bounds]
        report["cases"][name] = case
        for attack_name, attack in attacks.items():
            print(f"{name}: {attack_name}", flush=True)
            try:
                cost = attack(params, red_cost_model=estimator.RC.BDGL16)
                bits = float(cost["rop"].log2()) if hasattr(cost["rop"], "log2") else math.log2(cost["rop"])
                case["attacks"][attack_name] = dict(bits=bits, cost=str(cost))
                print(f"  {bits:.4f} bits", flush=True)
            except Exception as error:
                case["attacks"][attack_name] = {"error": repr(error)}
            case["complete"] = len(case["attacks"]) == len(attacks) and all(
                "bits" in item for item in case["attacks"].values())
            case["weakest_bits"] = min((item["bits"] for item in case["attacks"].values()
                                        if "bits" in item), default=None)
            case["passes_lwe_screen"] = case["complete"] and case["weakest_bits"] >= 128
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    return int(not all(case["passes_lwe_screen"] for case in report["cases"].values()))


if __name__ == "__main__":
    raise SystemExit(main())
