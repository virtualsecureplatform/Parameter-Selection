#!/usr/bin/env python3
"""Full lattice-estimator screen for the correlated lvl02 candidate."""

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
from tfhe_lvl02_correlated_candidate import analyse


def _log2_cost(value: Any) -> float:
    return float(log(RR(value), 2))


def _summarize(costs: dict[str, Any]) -> dict[str, Any]:
    algorithms = {
        name: {"rop_log2": _log2_cost(cost["rop"]), "summary": str(cost)}
        for name, cost in costs.items()
    }
    finite = [value["rop_log2"] for value in algorithms.values()]
    return {
        "minimum_rop_log2": min(finite) if finite else math.inf,
        "algorithms": algorithms,
    }


def _parameters(tag: str, instance: dict[str, Any], secret: str):
    if secret == "binary":
        secret_distribution = estimator.nd.Binary
    elif secret == "ternary":
        secret_distribution = estimator.nd.UniformMod(3)
    else:
        raise ValueError(f"unsupported secret distribution: {secret}")
    return estimator.lwe_parameters.LWEParameters(
        n=instance["n"],
        q=instance["q"],
        Xs=secret_distribution,
        Xe=estimator.nd.DiscreteGaussian(stddev=RR(instance["sigma"])),
        m=instance["m"],
        tag=tag,
    )


def run(mode: str) -> dict[str, Any]:
    source = analyse(WORKSPACE / "TFHEpp")
    instances = source["hardness_instances"]
    cases = [
        _parameters("candidate-prefix", instances["prefix_first"], "binary"),
        _parameters(
            "candidate-known-prefix-suffix-lwe-proxy",
            instances["suffix_proxy"],
            "ternary",
        ),
        _parameters("candidate-input-tlwe", instances["input_tlwe"], "binary"),
    ]

    output = {
        "scope": {
            "mode": mode,
            "cost_model": "BDGL16",
            "individual_target_bits": 128,
            "reduction_adjusted_target_bits": 131,
            "warning": (
                "Heuristic lattice-estimator costs.  The suffix calculation "
                "forgets negacyclic structure and is not a formal reduction."
            ),
        },
        "source_candidate": source["candidate"],
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
        summary["meets_reduction_adjusted_cost"] = (
            summary["minimum_rop_log2"] >= 131
        )
        output["cases"][parameters.tag] = summary
        print(
            f"  minimum log2(rop): {summary['minimum_rop_log2']:.3f}",
            flush=True,
        )
    output["all_cases_meet_reduction_adjusted_cost"] = all(
        case["meets_reduction_adjusted_cost"]
        for case in output["cases"].values()
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
                f"adjusted-target={item['meets_reduction_adjusted_cost']}"
            )
        print(
            "all adjusted targets: "
            f"{result['all_cases_meet_reduction_adjusted_cost']}"
        )


if __name__ == "__main__":
    main()
