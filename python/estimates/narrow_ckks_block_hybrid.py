#!/usr/bin/env python3
"""Exact public-block guessing counts with a bounded hybrid-cost sensitivity probe.

The count model matches independent uniform positions and signs in each public
block. Lattice geometry uses average marginal variance and the pinned estimator's
heuristics. A fixed-parameter probe is not a complete attack optimization or a
security certificate. Run --self-test with ordinary Python; run the probe in Sage.
"""

import argparse
import datetime
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys


def multiply_linear(coefficients, constant, linear):
    result = [0] * (len(coefficients) + 1)
    for index, value in enumerate(coefficients):
        result[index] += value * constant
        result[index + 1] += value * linear
    return result


def guessing_counts(widths, guessed):
    """Return cumulative support size and exact Pr[guessed weight <= w].

    Candidate polynomial: product (1 + 2*g_j*x), except a fully guessed
    block forces its one signed nonzero, giving 2*g_j*x. Probability polynomial:
    product ((width_j-g_j) + g_j*x), divided by product width_j. Signs cancel
    from that probability because each full position vector has 2**h signs.
    """
    if not widths or len(widths) != len(guessed):
        raise ValueError("nonempty matching width and guessed vectors required")
    candidates, masses, denominator = [1], [1], 1
    for width, count in zip(widths, guessed):
        if width < 1 or count < 0 or count > width:
            raise ValueError("invalid public block or guessed coordinate count")
        candidates = multiply_linear(candidates, int(count < width), 2 * count)
        masses = multiply_linear(masses, width - count, count)
        denominator *= width
    assert sum(masses) == denominator
    result, support, mass = [], 0, 0
    for weight, (count, numerator) in enumerate(zip(candidates, masses)):
        support += count
        mass += numerator
        result.append((weight, support, Fraction(mass, denominator)))
    return result


def conditional_moments(widths, guessed, threshold):
    """Exact residual second moments conditional on guessed weight <= threshold.

    For a residual coordinate in block j, its conditional squared coefficient
    expectation is B_j/A: A counts all accepted position vectors, while B_j
    counts accepted position vectors of the other blocks after fixing that
    coordinate. Independent signs retain zero means under this conditioning.
    """
    guessing_counts(widths, guessed)  # validate the geometry
    if not 0 <= threshold <= len(widths):
        raise ValueError("weight threshold outside blocks")
    masses = [1]
    for width, count in zip(widths, guessed):
        masses = multiply_linear(masses, width - count, count)
    accepted = sum(masses[:threshold + 1])
    if not accepted:
        return None
    variances = []
    for index, (width, count) in enumerate(zip(widths, guessed)):
        if width == count:
            variances.append(None)
            continue
        others = [1]
        for j, (other_width, other_count) in enumerate(zip(widths, guessed)):
            if j != index:
                others = multiply_linear(others, other_width - other_count, other_count)
        variances.append(Fraction(sum(others[:threshold + 1]), accepted))
    norm = sum((width - count) * variance for width, count, variance in
               zip(widths, guessed, variances) if variance is not None)
    independently = Fraction(sum((len(widths) - k) * mass for k, mass in
                                 enumerate(masses[:threshold + 1])), accepted)
    assert norm == independently
    remaining = sum(widths) - sum(guessed)
    return {"remaining_norm_squared": norm,
            "average_remaining_coordinate_variance": norm / remaining if remaining else None,
            "remaining_coordinate_variance_per_block": variances}


def layout(widths, zeta, strategy):
    n = sum(widths)
    if not 0 <= zeta <= n:
        raise ValueError("guessing dimension outside ring")
    if strategy == "prefix":
        counts = []
        remaining = zeta
        for width in widths:
            count = min(remaining, width)
            counts.append(count)
            remaining -= count
    elif strategy == "balanced":
        counts = [zeta * width // n for width in widths]
        remainder = zeta - sum(counts)
        order = sorted(range(len(widths)), key=lambda j: (zeta * widths[j] % n, -j), reverse=True)
        for index in order[:remainder]:
            counts[index] += 1
    else:
        raise ValueError("unknown layout")
    assert sum(counts) == zeta and all(0 <= c <= w for c, w in zip(counts, widths))
    return counts


def self_test():
    cases = 0
    # Independently enumerate secrets and distinct guessed subvectors, including
    # forced nonzeros in fully guessed blocks and unequal final-block widths.
    for widths in ([2, 3, 1], [1, 1], [3]):
        choices = [[(position, sign) for position in range(width) for sign in (-1, 1)]
                   for width in widths]
        secrets = list(itertools.product(*choices))
        for guessed in itertools.product(*(range(width + 1) for width in widths)):
            counts = guessing_counts(widths, guessed)
            projections = []
            for secret in secrets:
                projected = tuple((position, sign) if position < count else None
                                  for (position, sign), count in zip(secret, guessed))
                projections.append(projected)
            for threshold, support, probability in counts:
                accepted = [p for p in projections if sum(v is not None for v in p) <= threshold]
                assert len(set(accepted)) == support
                assert Fraction(len(accepted), len(secrets)) == probability
                moments = conditional_moments(widths, guessed, threshold)
                if not accepted:
                    assert moments is None
                    continue
                accepted_secrets = [s for s, p in zip(secrets, projections)
                                    if sum(v is not None for v in p) <= threshold]
                norm = sum(sum(position >= count for (position, _), count in zip(s, guessed))
                           for s in accepted_secrets)
                assert moments["remaining_norm_squared"] == Fraction(norm, len(accepted_secrets))
                for j, (width, count) in enumerate(zip(widths, guessed)):
                    for position in range(count, width):
                        actual = Fraction(sum(s[j][0] == position for s in accepted_secrets), len(accepted_secrets))
                        assert moments["remaining_coordinate_variance_per_block"][j] == actual
            cases += 1
        for zeta in range(sum(widths) + 1):
            for strategy in ("balanced", "prefix"):
                layout(widths, zeta, strategy)
    assert guessing_counts([2], [0]) == [(0, 1, Fraction(1)), (1, 1, Fraction(1))]
    print("exact public-block counts and conditional moments passed exhaustive checks:", cases, "partitions", flush=True)


def hash_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def estimator_revision(root):
    try:
        return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    except FileNotFoundError:
        metadata = root / ".git"
        if metadata.is_file():
            metadata = (root / metadata.read_text().strip().removeprefix("gitdir: ")).resolve()
        head = (metadata / "HEAD").read_text().strip()
        if not head.startswith("ref: "):
            return head
        ref = head[5:]
        if (metadata / ref).exists():
            return (metadata / ref).read_text().strip()
        return next(line.split()[0] for line in (metadata / "packed-refs").read_text().splitlines()
                    if not line.startswith(("#", "^")) and line.split()[1] == ref)


def cost_record(result, log):
    bits = float(log(result["rop"], 2))
    return {"bits": bits if math.isfinite(bits) else None, "cost": str(result),
            "parameters": {key: int(result[key]) for key in ("beta", "zeta", "d", "eta") if key in result}}


def probe(args):
    sys.path.insert(0, str(args.estimator.resolve()))
    from estimator import LWE, ND, RC
    from estimator.lwe_primal import PrimalHybrid
    from sage.all import QQ, log, sqrt
    from narrow_ckks_security import configure_max_beta

    configure_max_beta(args.max_beta)
    n, h = args.n, args.h
    width = (n + h - 1) // h
    widths = [min(width, n - j * width) for j in range(h)]
    if min(widths) < 1:
        raise ValueError("public ceil-width geometry has empty blocks")
    baseline_secret = ND.SparseTernary((h + 1) // 2, h // 2, n=n)
    baseline_params = LWE.Parameters(n=n, q=2**args.log_q, Xs=baseline_secret,
                                     Xe=ND.DiscreteGaussian(args.sigma))
    # The actual sign distribution has zero mean. Average marginal variance is
    # h/n; per-block and conditional variances are not represented in this GSA.
    marginal_secret = ND.NoiseDistribution(n=n, mean=0, stddev=sqrt(QQ(h) / n),
                                          bounds=(-1, 1), _density=QQ(h) / n)
    structured_params = LWE.Parameters(n=n, q=2**args.log_q, Xs=marginal_secret,
                                       Xe=ND.DiscreteGaussian(args.sigma))
    revision = estimator_revision(args.estimator.resolve())
    saved_baselines = None
    if args.baseline_report is not None:
        saved_baselines = json.loads(args.baseline_report.read_text())
        for key, value in {"n": n, "h": h, "log_q": args.log_q, "sigma": args.sigma,
                           "estimator_revision": revision}.items():
            if saved_baselines.get(key) != value:
                raise ValueError("saved baseline mismatch: " + key)
    report = {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "scope": __doc__, "source_sha256": hash_file(__file__),
              "estimator_revision": revision, "n": n, "h": h, "log_q": args.log_q,
              "sigma": args.sigma, "block_widths": widths, "max_beta": args.max_beta,
              "geometry_model": args.geometry,
              "baseline_report": str(args.baseline_report) if args.baseline_report else None,
              "baseline_report_sha256": hash_file(args.baseline_report) if args.baseline_report else None,
              "exact_count_self_test_passed": True, "attacks": {},
              "security_validated": False, "complete_attack_optimization": False,
              "limitations": ["Fixed uniform-proxy beta/zeta/d parameters; no optimization over block layouts or lattice parameters",
                              "Scalar average variance used in GSA; heterogeneous variances, conditional tails and correlations are not modeled",
                              "Estimator Babai/MITM probabilities and square-root speedup remain heuristic",
                              "Ring-specific and circular/correlated-key attacks are not modeled"]}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    for attack, mitm in (("bdd_hybrid", False), ("bdd_mitm_hybrid", True)):
        if saved_baselines is None:
            baseline = LWE.primal_hybrid(baseline_params, mitm=mitm, babai=mitm, red_cost_model=RC.BDGL16)
            baseline_record = cost_record(baseline, log)
        else:
            baseline_record = saved_baselines["attacks"][attack]["uniform_baseline"]
        if baseline_record["bits"] is None:
            raise RuntimeError("baseline optimization returned no finite probe point")
        point = baseline_record["parameters"]
        entry = {"uniform_baseline": baseline_record, "layouts": {}}
        report["attacks"][attack] = entry
        print(attack, "baseline", baseline_record, flush=True)
        for strategy in ("balanced", "prefix"):
            guessed = layout(widths, point["zeta"], strategy)
            best = None
            for weight, support, probability in guessing_counts(widths, guessed):
                if not support or not probability:
                    continue
                moments = conditional_moments(widths, guessed, weight)
                variance = moments["average_remaining_coordinate_variance"]
                parameters = structured_params
                if args.geometry == "conditional-mean":
                    if not variance:
                        raise ValueError("zero residual variance requires a separate exhaustive-search model")
                    scalar_variance = QQ(variance.numerator) / variance.denominator
                    # Only the projected secret's scalar variance is changed.
                    # params.n retains the original n so beta_params projects
                    # exactly n-zeta secret columns. This is a GSA sensitivity
                    # input, not a full distribution for the unprojected vector.
                    residual_geometry = ND.NoiseDistribution(n=n, mean=0,
                        stddev=sqrt(scalar_variance), bounds=(-1, 1), _density=scalar_variance)
                    parameters = LWE.Parameters(n=n, q=2**args.log_q, Xs=residual_geometry,
                                                Xe=ND.DiscreteGaussian(args.sigma))
                result = PrimalHybrid.cost(beta=point["beta"], params=parameters,
                    zeta=point["zeta"], d=point["d"], babai=mitm, mitm=mitm,
                    red_cost_model=RC.BDGL16, search_space=QQ(support),
                    hit_probability=QQ(probability.numerator) / probability.denominator)
                record = cost_record(result, log)
                if record["bits"] is not None and (best is None or record["bits"] < best["bits"]):
                    best = {**record, "maximum_guessed_weight": weight,
                            "conditional_remaining_norm_squared": str(moments["remaining_norm_squared"]),
                            "conditional_average_remaining_variance": str(variance),
                            "conditional_remaining_variance_per_block": [str(v) if v is not None else None for
                                                                          v in moments["remaining_coordinate_variance_per_block"]],
                            "candidate_count": str(support),
                            "hit_probability_numerator": str(probability.numerator),
                            "hit_probability_denominator": str(probability.denominator)}
            entry["layouts"][strategy] = {"guessed_per_block": guessed, "best_fixed_parameter_cost": best}
            print(attack, strategy, None if best is None else best["bits"], flush=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    report["fixed_parameter_probe_complete"] = True
    args.output.write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--estimator", type=Path, default=Path(__file__).resolve().parents[1] / "lattice-estimator")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--n", type=int, default=16384)
    parser.add_argument("--h", type=int, default=31)
    parser.add_argument("--log-q", type=int, default=76)
    parser.add_argument("--sigma", type=float, default=3.2)
    parser.add_argument("--max-beta", type=int, default=4096)
    parser.add_argument("--baseline-report", type=Path, help="reuse verified saved uniform baseline points")
    parser.add_argument("--geometry", choices=("unconditional-mean", "conditional-mean"), default="unconditional-mean",
                        help="scalar residual variance used in the fixed-parameter GSA sensitivity probe")
    args = parser.parse_args()
    if (args.n < 1 or not 1 <= args.h <= args.n or args.log_q < 1 or
            not math.isfinite(args.sigma) or args.sigma <= 0 or args.max_beta < 41):
        parser.error("positive valid dimension, weight, modulus and noise required")
    self_test()
    if args.self_test:
        return 0
    if args.output is None:
        parser.error("--output is required for the Sage probe")
    probe(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
