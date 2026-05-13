#!/usr/bin/env python3
"""Command-line BFV average-case noise estimator."""

from __future__ import annotations

import argparse
from dataclasses import replace
import math
import sys
from pathlib import Path

# Allow running as a script from either the repository root or python/.
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent
sys.path[:0] = [str(_REPO_ROOT), str(_SCRIPT_DIR)]

try:
    from python.noiseestimation.bfv import (  # noqa: E402
        BFVParams,
        NoiseState,
        correctness_d_for_failure,
        estimate_polyeval_bsgs,
        failure_log2_from_variance,
        find_bsgs_params,
        format_log2,
        fresh,
        log2_correctness_threshold,
    )
    from python.noiseestimation.bfv_polynomial import (  # noqa: E402
        lowest_digit_removal_polynomial_over_range,
        prime_from_prime_square,
    )
    from python.noiseestimation.params.bfv import PRESETS  # noqa: E402
except ModuleNotFoundError as exc:
    if exc.name != "python":
        raise
    from noiseestimation.bfv import (  # type: ignore[no-redef] # noqa: E402
        BFVParams,
        NoiseState,
        correctness_d_for_failure,
        estimate_polyeval_bsgs,
        failure_log2_from_variance,
        find_bsgs_params,
        format_log2,
        fresh,
        log2_correctness_threshold,
    )
    from noiseestimation.bfv_polynomial import (  # type: ignore[no-redef] # noqa: E402
        lowest_digit_removal_polynomial_over_range,
        prime_from_prime_square,
    )
    from noiseestimation.params.bfv import PRESETS  # type: ignore[no-redef] # noqa: E402


def parse_int_range(text: str) -> list[int]:
    values: list[int] = []
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = [int(p) for p in chunk.split(":")]
        if len(parts) == 1:
            values.append(parts[0])
        elif len(parts) in {2, 3}:
            start, stop = parts[0], parts[1]
            step = parts[2] if len(parts) == 3 else 1
            if step == 0:
                raise argparse.ArgumentTypeError("range step must be nonzero")
            if (stop - start) * step < 0:
                raise argparse.ArgumentTypeError(f"empty range: {chunk}")
            values.extend(range(start, stop + (1 if step > 0 else -1), step))
        else:
            raise argparse.ArgumentTypeError(f"bad range: {chunk}")
    if not values:
        raise argparse.ArgumentTypeError("empty range")
    return values


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Estimate average-case BFV invariant noise using 600.pdf formulas."
    )
    ap.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        default="tfhepp-lvl3simd-boot",
        help="parameter preset; boot uses TFHEpp PrimePower2Param with plaintext modulus p^2",
    )
    ap.add_argument("--qbits", type=int, help="single ciphertext modulus bit size")
    ap.add_argument(
        "--qbits-range",
        type=parse_int_range,
        help="inclusive range start:stop:step, or comma-separated values",
    )
    ap.add_argument("--B", type=int, default=15, help="digit-error bound")
    ap.add_argument(
        "--B-range",
        type=parse_int_range,
        help="inclusive B range; degree defaults to 4*B+1 for each B",
    )
    ap.add_argument("--degree", type=int, help="override polynomial degree")
    ap.add_argument("--nbit", type=int, help="override ring degree bit length")
    ap.add_argument("--plain-modulus", type=int, help="override BFV plaintext modulus t")
    ap.add_argument(
        "--alpha-log2",
        type=float,
        help="override normalized torus fresh-noise stddev for TFHEpp-style symmetric presets",
    )
    ap.add_argument(
        "--error-std",
        type=float,
        help="override integer coefficient error stddev before BFV t/q scaling",
    )
    ap.add_argument(
        "--fresh",
        choices=("symmetric", "public"),
        help="override fresh encryption model",
    )
    ap.add_argument(
        "--key-switch",
        choices=("ghs", "ghs-rns", "hybrid", "hybrid-rns"),
        help="override paper key-switch/relinearization noise variant",
    )
    ap.add_argument("--rns-digits", type=int, help="override RNS digit count")
    ap.add_argument("--hybrid-omega", type=int, help="override hybrid decomposition omega")
    ap.add_argument(
        "--scalar-mode",
        choices=(
            "none",
            "unsigned-average",
            "centered-average",
            "poly-average",
            "unsigned-exact",
            "centered-exact",
        ),
        default="unsigned-exact",
        help="coefficient amplification model for nonconstant PolyEval terms",
    )
    ap.add_argument(
        "--circuit-model",
        choices=("dependent", "independent"),
        default="dependent",
        help="dependent matches PolyEval powers of one ciphertext; independent matches 600.pdf Section 4",
    )
    ap.add_argument(
        "--poly-source",
        choices=("auto", "degree", "digit-removal"),
        default="auto",
        help="auto uses TFHEpp's bounded digit-removal polynomial for *-boot presets",
    )
    ap.add_argument(
        "--digit-prime",
        type=int,
        help="prime p for the p^2 digit-removal plaintext ring; defaults to sqrt(t)",
    )
    ap.add_argument(
        "--input-log2-var",
        type=float,
        help="override input ciphertext invariant-noise variance",
    )
    ap.add_argument("--input-degree", type=int, default=1)
    ap.add_argument(
        "--D",
        type=float,
        default=6.0,
        help="Gaussian tail parameter for V <= 1/(8*D^2)",
    )
    ap.add_argument(
        "--failure-log2",
        type=float,
        help="derive D from target total failure probability 2^x over all coefficients",
    )
    ap.add_argument("--quiet", action="store_true")
    return ap.parse_args()


def make_params(args: argparse.Namespace, q_bits: int | None) -> BFVParams:
    factory = PRESETS[args.preset]
    if args.preset == "paper-u3":
        params = factory(
            nbit=args.nbit if args.nbit is not None else 13,
            q_bits=q_bits if q_bits is not None else 149,
            t=args.plain_modulus if args.plain_modulus is not None else 65_537,
            error_std=args.error_std if args.error_std is not None else 3.19,
        )
    elif args.preset == "openfhe-paper":
        params = factory(
            nbit=args.nbit if args.nbit is not None else 13,
            q_bits=q_bits if q_bits is not None else 60,
            t=args.plain_modulus if args.plain_modulus is not None else 65_537,
            error_std=args.error_std if args.error_std is not None else 3.19,
        )
    else:
        kwargs = {}
        if q_bits is not None:
            kwargs["q_bits"] = q_bits
        if args.alpha_log2 is not None:
            kwargs["alpha_log2"] = args.alpha_log2
        params = factory(**kwargs)

    overrides = {}
    if args.nbit is not None:
        overrides["nbit"] = args.nbit
    if args.plain_modulus is not None:
        overrides["t"] = args.plain_modulus
    if args.error_std is not None:
        overrides["error_log2_std"] = math.log2(args.error_std)
    if args.fresh is not None:
        overrides["fresh"] = args.fresh
    if args.key_switch is not None:
        overrides["key_switch"] = args.key_switch
    if args.rns_digits is not None:
        overrides["rns_digits"] = args.rns_digits
    if args.hybrid_omega is not None:
        overrides["hybrid_omega"] = args.hybrid_omega
    return replace(params, **overrides) if overrides else params


def input_state(args: argparse.Namespace, params: BFVParams) -> NoiseState:
    if args.input_log2_var is None:
        return fresh(params)
    return NoiseState(args.input_log2_var, args.input_degree, "input")


def polynomial_coeffs(
    args: argparse.Namespace, params: BFVParams, b: int
) -> tuple[list[int] | None, str]:
    source = args.poly_source
    if source == "auto":
        source = "digit-removal" if args.preset.endswith("-boot") else "degree"
    if source == "degree":
        return None, "degree"
    p = (
        args.digit_prime
        if args.digit_prime is not None
        else prime_from_prime_square(params.t)
    )
    coeffs = lowest_digit_removal_polynomial_over_range(p, b)
    return coeffs, f"digit-removal(p={p})"


def main() -> int:
    args = parse_args()
    if args.qbits is not None:
        q_values: list[int | None] = [args.qbits]
    elif args.qbits_range is not None:
        q_values = args.qbits_range
    else:
        q_values = [None]
    b_values = args.B_range if args.B_range is not None else [args.B]
    if args.degree is not None and args.B_range is not None:
        raise SystemExit("--degree cannot be combined with --B-range")

    first_ok: tuple[int, int, float] | None = None

    if not args.quiet:
        print("BFV average-case invariant-noise estimate")
        print("  source: 600.pdf Sections 4-5 and Section 7 dependent-circuit bounds")
        if args.circuit_model == "dependent":
            print(
                "  dependent mode omits Var((nu*nu')|i), "
                "as in the paper's identical-input examples"
            )
        print("  security: not estimated by this script")
        print(
            f"  preset={args.preset} circuit_model={args.circuit_model} "
            f"scalar_mode={args.scalar_mode} poly_source={args.poly_source}"
        )

    header = (
        "q_bits  B   degree  nbit  log2(t)  poly  bsgs(k,m)  "
        "log2(V)  threshold  margin  log2(p_fail)  status"
    )
    print(header)
    print("-" * len(header))

    for q_bits in q_values:
        for b in b_values:
            params = make_params(args, q_bits)
            if args.degree is not None:
                coeffs = None
                poly_label = "degree"
                degree = args.degree
            else:
                coeffs, poly_label = polynomial_coeffs(args, params, b)
                degree = len(coeffs) - 1 if coeffs is not None else 4 * b + 1
            x = input_state(args, params)
            state = estimate_polyeval_bsgs(
                x,
                degree,
                params,
                scalar_mode=args.scalar_mode,
                circuit_model=args.circuit_model,
                coeffs=coeffs,
            )
            d = (
                correctness_d_for_failure(params.n, args.failure_log2)
                if args.failure_log2 is not None
                else args.D
            )
            threshold = log2_correctness_threshold(d)
            margin = threshold - state.log2_variance
            p_fail = failure_log2_from_variance(params.n, state.log2_variance)
            status = "OK" if margin >= 0 else "FAIL"
            k, m = find_bsgs_params(degree)
            poly_short = (
                "digit" if poly_label.startswith("digit-removal") else "degree"
            )
            print(
                f"{params.q_bits:6d} {b:3d} {degree:8d} {params.nbit:5d} "
                f"{math.log2(params.t):8.2f}  {poly_short:>5} ({k:3d},{m:1d})    "
                f"{format_log2(state.log2_variance):>8} "
                f"{format_log2(threshold):>10} "
                f"{format_log2(margin):>7} "
                f"{format_log2(p_fail):>13}  {status}"
            )
            if status == "OK" and first_ok is None:
                first_ok = (params.q_bits, b, margin)

    if first_ok is not None:
        q_bits, b, margin = first_ok
        print(f"\nfirst OK: q_bits={q_bits}, B={b}, margin={margin:.2f} bits")
    else:
        print("\nfirst OK: none in sweep")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
