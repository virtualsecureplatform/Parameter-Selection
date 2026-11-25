#!/usr/bin/env python3
"""Lattice Estimator screens for the keys a TFHEpp TUTTI schedule publishes.

Two secrets are exposed (ePrint 2026/2848, Section 3):

* the dense ternary key, at the largest evaluation-key modulus and degree N;
* the block-binary key, only at the base modulus. A subring key s'(X^R) is an
  instance of degree N/R; fixing s_0 = 1 removes one block, which leaves
  dimension N/R - B and weight h - 1.

The sparse instance is screened as a uniform sparse binary secret: the
estimator models neither the one-per-block constraint nor ring rotations,
and the paper's rotated primal hybrid caps are reported next to it as an
external input. Combinatorial attacks are costed separately by
`noiseestimation.block_binary_meet_lwe`. A passing screen is not a security
proof; circular security is assumed. Run with `sage -python`.
"""

import argparse
import concurrent.futures
import datetime
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from noiseestimation.tutti import PROFILES, budget, published_caps  # noqa: E402

DENSE_ATTACKS = ("primal_usvp", "primal_bdd", "dual", "dual_hybrid")
HYBRID_ATTACKS = ("bdd_hybrid", "bdd_mitm_hybrid", "dual_mitm_hybrid")


def secret(kind, dimension, weight):
    from estimator import ND

    if kind == "dense":
        if weight == 0:
            return ND.UniformMod(3)
        return ND.SparseTernary(weight // 2, weight - weight // 2, n=dimension)
    return ND.SparseBinary(weight, n=dimension)


def run(task):
    root, kind, dimension, weight, bits, sigma, attack, samples, max_beta = task
    sys.path.insert(0, root)
    from narrow_ckks_security import configure_max_beta
    from estimator import LWE, ND
    from estimator.lwe_dual import dual_hybrid
    from sage.all import log, oo

    configure_max_beta(max_beta)
    params = LWE.Parameters(n=dimension, q=2**bits, Xs=secret(kind, dimension, weight),
                            Xe=ND.DiscreteGaussian(sigma),
                            m=oo if samples is None else samples)
    try:
        if attack == "bdd_hybrid":
            result = LWE.primal_hybrid(params, mitm=False, babai=False)
        elif attack == "bdd_mitm_hybrid":
            result = LWE.primal_hybrid(params, mitm=True, babai=True)
        elif attack == "dual_mitm_hybrid":
            result = dual_hybrid(params, mitm_optimization=True)
        else:
            result = getattr(LWE, attack)(params)
        cost = float(log(result["rop"], 2))
        if not math.isfinite(cost):
            return task, {"status": "nonfinite_optimizer_result", "cost": str(result)}
        return task, {"bits": cost, "cost": str(result)}
    except Exception as exc:  # the estimator raises on some degenerate cases
        return task, {"error": f"{type(exc).__name__}: {exc}"}


def cases_for(args):
    cases = {}
    if args.profile:
        p = PROFILES[args.profile]
        exposure = budget(p)
        anchored = p.sparse_degree - p.sparse_degree // p.blocks
        cases["dense"] = ("dense", p.n, p.dense_weight,
                          exposure["maximum_dense_key_bits"], p.sigma)
        cases["sparse"] = ("sparse", anchored, p.blocks - 1, p.base_bits, p.sigma)
    for spec in args.dense or ():
        n, weight, bits = spec
        cases[f"dense_n{n}_w{weight}_q{bits}"] = ("dense", n, weight, bits, args.noise_stddev)
    for spec in args.sparse or ():
        degree, blocks, bits = spec
        cases[f"sparse_n{degree}_h{blocks}_q{bits}"] = (
            "sparse", degree - degree // blocks, blocks - 1, bits, args.noise_stddev)
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--estimator", type=Path,
                        default=Path(__file__).resolve().parents[1] / "lattice-estimator")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", choices=sorted(PROFILES))
    parser.add_argument("--dense", type=int, nargs=3, action="append",
                        metavar=("N", "WEIGHT", "LOGQ"),
                        help="extra dense case; weight 0 is uniform ternary")
    parser.add_argument("--sparse", type=int, nargs=3, action="append",
                        metavar=("DEGREE", "BLOCKS", "LOGQ"),
                        help="extra anchored block-binary case at degree N/R")
    parser.add_argument("--noise-stddev", type=float, default=3.2,
                        help="integer noise RMS of the extra cases")
    parser.add_argument("--one-sample", action="store_true",
                        help="give the attack one RLWE sample (m = dimension) instead of unlimited")
    parser.add_argument("--skip-hybrids", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--max-beta", type=int)
    args = parser.parse_args()
    cases = cases_for(args)
    if not cases:
        parser.error("give --profile, --dense or --sparse")
    root = str(args.estimator.resolve())
    git_dir = args.estimator / ".git"
    revision = None
    if git_dir.is_file():
        git_dir = (args.estimator / git_dir.read_text().strip().removeprefix("gitdir: ")).resolve()
    head = git_dir / "HEAD"
    if head.exists():
        revision = head.read_text().strip()
    report = {
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "scope": __doc__.splitlines()[0],
        "estimator_revision": revision,
        "cost_model": "estimator default (RC.MATZOV classical), rop",
        "samples": "one RLWE sample" if args.one_sample else "unbounded",
        "security_validated": False,
        "excluded": ["block constraint in lattice attacks", "ring rotations (rotated primal hybrid)",
                     "circular/KDM security", "exact error distribution"],
        "cases": {},
    }
    if args.profile:
        report["profile"] = args.profile
        report["published_caps"] = published_caps(PROFILES[args.profile])
    tasks = []
    for name, (kind, dimension, weight, bits, sigma) in cases.items():
        hybrids = not args.skip_hybrids and (kind == "sparse" or weight != 0)
        attacks = DENSE_ATTACKS + (HYBRID_ATTACKS if hybrids else ())
        report["cases"][name] = {
            "kind": kind, "dimension": dimension, "log_q": bits, "sigma": sigma,
            "secret": ("uniform ternary" if kind == "dense" and weight == 0 else
                       f"sparse ternary, weight {weight}" if kind == "dense" else
                       f"uniform sparse binary proxy, weight {weight}"),
            "requested_attacks": list(attacks), "attacks": {}}
        samples = dimension if args.one_sample else None
        tasks += [(root, kind, dimension, weight, bits, sigma, attack, samples, args.max_beta)
                  for attack in attacks]
    names = {(case[0], case[1], case[2], case[3]): name for name, case in cases.items()}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
        for future in concurrent.futures.as_completed([executor.submit(run, t) for t in tasks]):
            task, result = future.result()
            name = names[task[1:5]]
            report["cases"][name]["attacks"][task[6]] = result
            print(name, task[6], result.get("bits", result.get("error", result.get("status"))),
                  flush=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    complete = True
    for case in report["cases"].values():
        finished = [a for a in case["requested_attacks"] if "bits" in case["attacks"].get(a, {})]
        case["complete"] = len(finished) == len(case["requested_attacks"])
        complete = complete and case["complete"]
        if finished:
            weakest = min(finished, key=lambda a: case["attacks"][a]["bits"])
            case["weakest_attack"] = weakest
            case["weakest_bits"] = case["attacks"][weakest]["bits"]
            case["passes_128_screen"] = case["complete"] and case["weakest_bits"] >= 128
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
