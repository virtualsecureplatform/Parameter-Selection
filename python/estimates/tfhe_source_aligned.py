#!/usr/bin/env python3
"""lattice-estimator proxies for the correlated source-aligned TFHE screen.

Run inside the repository's SageMath container, for example:

    singularity exec --no-home --pwd /workspace/Parameter-Selection \
      --env DOT_SAGE=/tmp/.sage \
      --bind /path/to/workspace:/workspace python/sagemath.sif \
      sage -python python/estimates/tfhe_source_aligned.py --mode full

The suffix case deliberately treats the known-prefix structured suffix-RLWE
source as ordinary LWE on the unknown coefficients.  It is a conservative
heuristic proxy, not an equivalence theorem and not a replacement for the
formal RLWE assumption.
"""

from __future__ import annotations

import argparse
import importlib
import json
import math
import sys
from pathlib import Path
from typing import Any

from sage.all import RR, log


REPOSITORY = Path(__file__).resolve().parents[2]
WORKSPACE = REPOSITORY.parent
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "python"))
sys.path.insert(0, str(REPOSITORY / "python/proof"))

estimator = importlib.import_module(".estimator", "lattice-estimator")
from tfhe_source_aligned_screen import analyse


def _log2_cost(value: Any) -> float:
    return float(log(RR(value), 2))


def _summarize(costs: dict[str, Any]) -> dict[str, Any]:
    algorithms = {}
    for name, cost in costs.items():
        algorithms[name] = {
            "rop_log2": _log2_cost(cost["rop"]),
            "summary": str(cost),
        }
    finite = [item["rop_log2"] for item in algorithms.values()]
    return {
        "minimum_rop_log2": min(finite) if finite else math.inf,
        "algorithms": algorithms,
    }


def _parameters(tag: str, n: int, q: int, m: int, sigma: int, secret: str):
    if secret == "binary":
        secret_distribution = estimator.nd.Binary
    elif secret == "ternary":
        secret_distribution = estimator.nd.UniformMod(3)
    else:
        raise ValueError(f"unsupported secret distribution: {secret}")
    return estimator.lwe_parameters.LWEParameters(
        n=n,
        q=q,
        Xs=secret_distribution,
        Xe=estimator.nd.DiscreteGaussian(stddev=sigma),
        m=m,
        tag=tag,
    )


def run(mode: str) -> dict[str, Any]:
    source = analyse(WORKSPACE / "TFHEpp")
    prefix = source["hardness_instances"]["prefix_first"]
    suffix = source["hardness_instances"]["suffix"]
    cases = [
        _parameters(
            "source-aligned-prefix-baseline",
            prefix["n"],
            prefix["q"],
            prefix["m"],
            prefix["error_sigma_baseline"],
            "binary",
        ),
        _parameters(
            "source-aligned-prefix-largest-screened-sigma",
            prefix["n"],
            prefix["q"],
            prefix["m"],
            prefix["error_sigma_largest_integer_passing_fresh_term_screen"],
            "binary",
        ),
        _parameters(
            "source-aligned-known-prefix-suffix-lwe-proxy",
            suffix["unknown_iid_ternary_suffix_coefficients"],
            suffix["q"],
            suffix["scalar_equations_if_all_negacyclic_rotations_are_exposed"],
            suffix["error_sigma"],
            "ternary",
        ),
    ]

    output = {
        "scope": {
            "mode": mode,
            "cost_model": "BDGL16",
            "warning": (
                "Heuristic lattice-estimator costs.  The suffix calculation "
                "forgets negacyclic structure and is not a formal reduction."
            ),
            "individual_128_bit_target": 128,
            "equal_allocation_target_for_three_factor_two_terms": 131,
        },
        "cases": {},
    }
    for parameters in cases:
        print(f"estimating {parameters.tag}: {parameters}", flush=True)
        if mode == "rough":
            costs = estimator.LWE.estimate.rough(parameters)
        else:
            costs = estimator.LWE.estimate(
                parameters, red_cost_model=estimator.RC.BDGL16
            )
        summary = _summarize(costs)
        summary["parameters"] = str(parameters)
        summary["meets_128_bit_cost"] = summary["minimum_rop_log2"] >= 128
        summary["meets_equal_reduction_allocation"] = (
            summary["minimum_rop_log2"] >= 131
        )
        output["cases"][parameters.tag] = summary
        print(
            f"  minimum log2(rop): {summary['minimum_rop_log2']:.3f}",
            flush=True,
        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("rough", "full"), default="rough")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = run(args.mode)
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print("=" * 72)
        for name, item in result["cases"].items():
            print(
                f"{name}: min log2(rop)={item['minimum_rop_log2']:.3f}, "
                f"128-bit={item['meets_128_bit_cost']}, "
                f"equal-allocation={item['meets_equal_reduction_allocation']}"
            )


if __name__ == "__main__":
    main()
