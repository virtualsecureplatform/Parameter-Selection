#!/usr/bin/env python3
"""Sage/Lattice Estimator screens for the actual narrow DD port key moduli.

Uniform sparse proxies do not model the public one-nonzero-per-block constraint,
ring-specific attacks, or circular/KDM assumptions. A passing screen is not a
complete security claim. Run with `sage -python` and an estimator checkout.
"""

import argparse
import concurrent.futures
import datetime
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys


CASES = {
    "dd_dense_selectors": (8192, 184, None),
    "dd_dense_packed_setup": (8192, 217, None),
    "dd_structured_source_uniform_sparse_proxy": (8192, 25, 31),
    "dd_ll13_auxiliary_dense_selectors": (8192, 187, None),
    "dd_ll13_auxiliary_dense_factor_keys": (8192, 198, None),
    "dd_ll13_auxiliary_sparse_switch_uniform_proxy": (8192, 53, 31),
    "dd_ll13_fbs_dense_setup": (8192, 217, None),
    "dd_ll13_fbs_sparse_switch_uniform_proxy": (8192, 49, 31),
    "paper_ll13_dense_reference": (8192, 218, None),
    "paper_ll13_sparse_reference_proxy": (8192, 55, 31),
    "dd_ll14_fbs_dense_selectors": (16384, 387, None),
    "dd_ll14_fbs_dense_packed_setup": (16384, 435, None),
    "dd_ll14_dense_selectors": (16384, 379, None),
    # Direct440 setup: use --noise-stddev 102.00190833992261 to match
    # the integer RMS of continuous sigma=102.4; the default is not this setup.
    "dd_ll14_basis_direct_selectors": (16384, 440, None),
    # Current CPU replay: configured continuous sigma3.2 truncates to integer
    # RMS2.828301687834186; supply it with --noise-stddev.
    "dd_ll14_native_automorphism_keys": (16384, 435, None),
    "dd_ll14_direct_dd_master": (16384, 433, None),
    "dd_ll14_dense_packed_setup": (16384, 427, None),
    "dd_ll14_structured_source_uniform_sparse_proxy": (16384, 48, 31),
    "dd_ll14_auxiliary_sparse_switch_uniform_proxy": (16384, 76, 31),
    "dd_ll14_plain37_aux16_sparse_switch_uniform_proxy": (16384, 64, 31),
    "dd_ll14_plain37_aux8_sparse_switch_uniform_proxy": (16384, 56, 31),
    "dd_ht14_auxiliary_dense_packed_setup": (16384, 420, None),
    "dd_ht14_auxiliary_sparse_switch_uniform_proxy": (16384, 86, 31),
    "dd_ht15_auxiliary_dense_factor_keys": (32768, 842, None),
    "dd_ht15_auxiliary_sparse_switch_uniform_proxy": (32768, 262, 31),
    "dd_basis_repeatable_dense_packed_setup": (32768, 590, None),
    "dd_basis_repeatable_dense_factor_keys": (32768, 541, None),
    "dd_basis_repeatable_sparse_switch_uniform_proxy": (32768, 92, 31),
    "dd_lean_basis_dense_packed_setup": (32768, 572, None),
    "dd_lean_basis_dense_factor_keys": (32768, 523, None),
    "dd_ll14_plain37_dense_packed_setup": (16384, 436, None),
}
ATTACKS = ("primal_usvp", "primal_bdd", "dual", "dual_hybrid")
OPTIONAL_ATTACKS = ("mitm", "bdd_hybrid", "bdd_mitm_hybrid")
COST_MODELS = ("BDGL16", "MATZOV")


def configure_max_beta(value):
    # The pinned estimator binds conf.max_beta when these modules are imported.
    # Updating conf alone after importing estimator would leave stale limits.
    from estimator import conf

    if value is not None:
        conf.max_beta = value
        for name in ("lwe_primal", "lwe_dual"):
            importlib.import_module("estimator." + name).max_beta_global = value
    limits = {name: int(importlib.import_module("estimator." + name).max_beta_global)
              for name in ("lwe_primal", "lwe_dual")}
    if any(limit != int(conf.max_beta) for limit in limits.values()):
        raise RuntimeError("estimator module block-size limits differ from configuration")
    return int(conf.max_beta), limits


def estimate(task):
    root, name, attack, max_beta, noise_stddev, cost_model = task
    sys.path.insert(0, root)
    from estimator import LWE, ND, RC
    from sage.all import log

    _, limits = configure_max_beta(max_beta)
    reduction_cost = getattr(RC, cost_model)
    n, bits, h = CASES[name]
    secret = ND.UniformMod(3) if h is None else ND.SparseTernary((h + 1) // 2, h // 2, n=n)
    params = LWE.Parameters(n=n, q=2**bits, Xs=secret,
                            Xe=ND.DiscreteGaussian(noise_stddev), tag=name)
    try:
        if attack == "mitm":
            result = LWE.mitm(params)
        elif attack in ("bdd_hybrid", "bdd_mitm_hybrid"):
            use_mitm = attack == "bdd_mitm_hybrid"
            result = LWE.primal_hybrid(params, mitm=use_mitm, babai=use_mitm,
                                      red_cost_model=reduction_cost)
        else:
            result = getattr(LWE, attack)(params, red_cost_model=reduction_cost)
        cost = float(log(result["rop"], 2))
        if not math.isfinite(cost):
            # Preserve the optimizer output. An infinite search result is
            # distinct from an estimator exception and is not a security bound.
            kind = "nan" if math.isnan(cost) else "positive_infinity" if cost > 0 else "negative_infinity"
            return name, attack, {
                "error": "ValueError: nonfinite attack cost",
                "status": "nonfinite_optimizer_result",
                "nonfinite_kind": kind,
                "rop": str(result["rop"]),
                "cost": str(result),
                "estimator_module_max_beta": limits,
            }
        return name, attack, {"bits": cost, "cost": str(result),
                             "estimator_module_max_beta": limits}
    except Exception as exc:
        return name, attack, {"error": f"{type(exc).__name__}: {exc}"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--estimator", type=Path,
        default=Path(__file__).resolve().parents[1] / "lattice-estimator")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=CASES, action="append")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--noise-stddev", type=float, default=3.2,
        help="integer-noise RMS for the discrete-Gaussian proxy; does not model the exact sampler distribution")
    parser.add_argument("--cost-model", choices=COST_MODELS, default="BDGL16",
        help="classical lattice reduction cost model; historical default is preserved")
    parser.add_argument("--attack", choices=ATTACKS + OPTIONAL_ATTACKS, action="append",
        help="run selected attacks; default is the four original lattice screens")
    parser.add_argument("--max-beta", type=int,
        help="override the pinned estimator block-size search limit in both primal and dual modules")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("jobs must be positive")
    if not math.isfinite(args.noise_stddev) or args.noise_stddev <= 0:
        parser.error("noise standard deviation must be finite and positive")
    if args.max_beta is not None and args.max_beta < 41:
        parser.error("max-beta must be at least 41")
    root = args.estimator.resolve()
    try:
        revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except FileNotFoundError:
        # Minimal Sage images omit git; retain provenance from checkout metadata.
        git_dir = root / ".git"
        if git_dir.is_file():
            git_dir = (root / git_dir.read_text().strip().removeprefix("gitdir: ")).resolve()
        revision = (git_dir / "HEAD").read_text().strip()
        if revision.startswith("ref: "):
            ref = revision[5:]
            loose = git_dir / ref
            if loose.exists():
                revision = loose.read_text().strip()
            else:
                revision = next(line.split()[0] for line in
                    (git_dir / "packed-refs").read_text().splitlines()
                    if not line.startswith(("#", "^")) and line.split()[1] == ref)
    names = args.case or list(CASES)
    attacks = list(dict.fromkeys(args.attack or ATTACKS))
    sys.path.insert(0, str(root))
    configured_max_beta, module_limits = configure_max_beta(args.max_beta)
    report = {
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "estimator_source": "https://github.com/malb/lattice-estimator",
        "estimator_revision": revision,
        "cost_model": "RC." + args.cost_model + " classical",
        "cost_model_name": args.cost_model,
        "samples": "unbounded",
        "estimator_global_max_beta": configured_max_beta,
        "estimator_module_max_beta": module_limits,
        "estimator_max_beta_override": args.max_beta,
        "requested_attacks": attacks,
        "original_four_attacks": list(ATTACKS),
        "mitm_cost_model": "estimator analytical word-operation model" if "mitm" in attacks else None,
        "nonfinite_costs_establish_security": False,
        "scope": "Generic LWE and uniform-sparse proxies; excludes public block constraints, ring structure, and circular/KDM security",
        "security_validated": False,
        "noise_proxy": "DiscreteGaussian with supplied integer RMS; exact sampler distribution is not modeled",
        "noise_stddev": args.noise_stddev,
        "cases": {},
    }
    for name in names:
        dimension, bits, weight = CASES[name]
        report["cases"][name] = {"n": dimension, "log_q": bits, "sigma": args.noise_stddev,
            "secret": "dense ternary" if weight is None else "uniform sparse ternary proxy (16 positive, 15 negative)",
            "attacks": {}}
        if weight is not None:
            width = (dimension + weight - 1) // weight
            sizes = [min(width, dimension - x * width) for x in range(weight)]
            report["cases"][name]["public_block_secret_entropy_bits"] = sum(
                math.log2(2 * size) for size in sizes)
    tasks = [(str(root), name, attack, args.max_beta, args.noise_stddev, args.cost_model)
             for name in names for attack in attacks]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
        for future in concurrent.futures.as_completed([executor.submit(estimate, task) for task in tasks]):
            name, attack, result = future.result()
            report["cases"][name]["attacks"][attack] = result
            print(name, attack, result.get("bits", result.get("error")), flush=True)
            # Preserve partial authoritative results during long sparse screens.
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    for case in report["cases"].values():
        case["requested_attacks_complete"] = all(
            attack in case["attacks"] and "bits" in case["attacks"][attack] for attack in attacks)
        case["complete"] = all(attack in case["attacks"] for attack in ATTACKS) and case["requested_attacks_complete"]
        if case["requested_attacks_complete"]:
            case["weakest_requested_attack_bits"] = min(result["bits"] for result in case["attacks"].values())
        if case["complete"]:
            case["weakest_bits"] = case["weakest_requested_attack_bits"]
            case["passes_generic_screen"] = case["weakest_bits"] >= 128
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0 if all(case["requested_attacks_complete"] for case in report["cases"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
