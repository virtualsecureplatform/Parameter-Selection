#!/usr/bin/env python3
"""Check the fixed-weight Renyi mean-energy condition for TFHEpp BFV.

The conditional Renyi route in ``sketch/spanningtree.md`` requires the
descriptor-conditioned whitened energy ``F`` to leave a positive first-order
entropy margin.  Before the linear masked-ratio term and every concentration
penalty, this requires

    E[F_quadratic] < 2 log |T_{n,w}|,

where ``T_{n,w}`` is the exact-weight signed ternary secret support.

For one ordinary HalfTRGSW row with scalar gadget weight ``g`` and spherical
coefficient error of standard deviation ``sigma``, the exact fixed-weight
square covariance gives

    E[F_quadratic] = (g / sigma)^2 E[||S^2||_2^2].

The hidden HNF denominator multiplies both the quadratic message and the
terminal error, so it cancels from this ratio.  Consequently this script can
certify a masked-ratio-descriptor-independent obstruction for the uniform
fixed-weight prior, without selecting a DSPR/NTRU descriptor distribution.
It does not claim that an arbitrary leakage-conditioned posterior has the
same moment.  The screened channel consists of the ``l`` independent ordinary
rows created before public Double Decomposition.  Deterministic decomposition
preserves their independence, but failure of an upstream sufficient bound
does not by data processing alone rule out a sharper analysis of the final
encoded rows; transferring the obstruction needs an exact reconstruction or
Renyi-invariance argument for that encoding.

All moment, gadget, and obstruction comparisons are performed with Python
integers or ``Fraction``.  Floating-point logarithms are reporting aids only.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import combinations, product
from typing import Any


@dataclass(frozen=True)
class Parameters:
    """Parameters entering the descriptor-independent quadratic energy."""

    degree: int = 1 << 15
    secret_weight: int = 96
    modulus_bits: int = 640
    error_std_log2: int = 33
    gadget_base_bits: int = 18
    gadget_levels: int = 35
    target_bits: int = 128

    def validate(self) -> None:
        if self.degree <= 0 or self.degree % 2 != 0:
            raise ValueError("degree must be a positive even integer")
        if not 0 < self.secret_weight <= self.degree:
            raise ValueError("secret weight must lie in [1, degree]")
        if self.modulus_bits <= 0:
            raise ValueError("modulus bits must be positive")
        if self.error_std_log2 < 0:
            raise ValueError("this exact checker expects sigma = 2^k, k >= 0")
        if self.gadget_base_bits <= 0 or self.gadget_levels <= 0:
            raise ValueError("gadget base bits and levels must be positive")
        if self.gadget_levels * self.gadget_base_bits > self.modulus_bits:
            raise ValueError("the final gadget weight would have a negative exponent")
        if self.target_bits <= 0:
            raise ValueError("target bits must be positive")


def fixed_weight_square_second_moment(degree: int, weight: int) -> Fraction:
    """Return the exact ``E[||S^2||_2^2]`` for uniform ``T_{n,w}``.

    For even negacyclic degree, the covariance proved in FormalProof4FHE is

      2*n*p2 I + 2*(p - 3*p2) P_even,

    with ``rank(P_even) = n/2``.  The square has mean zero, so taking the
    trace gives the returned second moment.
    """

    if degree <= 0 or degree % 2 != 0:
        raise ValueError("degree must be a positive even integer")
    if not 0 < weight <= degree:
        raise ValueError("weight must lie in [1, degree]")
    p = Fraction(weight, degree)
    p2 = Fraction(weight * (weight - 1), degree * (degree - 1))
    return 2 * degree * degree * p2 + degree * (p - 3 * p2)


def gadget_exponents(params: Parameters) -> list[int]:
    """Return TFHEpp's ordinary HalfTRGSW gadget exponents."""

    return [
        params.modulus_bits - (level + 1) * params.gadget_base_bits
        for level in range(params.gadget_levels)
    ]


def signed_fixed_weight_support_size(degree: int, weight: int) -> int:
    """Cardinality ``2^w * choose(n,w)`` of the proof secret support."""

    return (1 << weight) * math.comb(degree, weight)


def log2_fraction(value: Fraction) -> float:
    if value <= 0:
        raise ValueError("logarithm requires a positive fraction")
    return math.log2(value.numerator) - math.log2(value.denominator)


def binomial_conditioning_log2_penalty(degree: int, weight: int) -> float:
    """Return ``log2(1 / Pr[Bin(n,w/n)=w])`` from equation (15)."""

    p = weight / degree
    log_mass = (
        math.lgamma(degree + 1)
        - math.lgamma(weight + 1)
        - math.lgamma(degree - weight + 1)
        + weight * math.log(p)
        + (degree - weight) * math.log1p(-p)
    )
    return -log_mass / math.log(2)


def analyse(params: Parameters) -> dict[str, Any]:
    params.validate()
    moment = fixed_weight_square_second_moment(
        params.degree, params.secret_weight
    )
    exponents = gadget_exponents(params)
    squared_weight_sum = sum(1 << (2 * exponent) for exponent in exponents)
    sigma_squared = 1 << (2 * params.error_std_log2)
    total_energy = moment * Fraction(squared_weight_sum, sigma_squared)
    top_energy = moment * Fraction(1 << (2 * exponents[0]), sigma_squared)

    support_size = signed_fixed_weight_support_size(
        params.degree, params.secret_weight
    )
    entropy_bits = math.log2(support_size)
    entropy_nats = math.log(support_size)
    twice_entropy_nats = 2 * entropy_nats

    # Exact elementary entropy bound: for b = bit_length(M),
    # log(M) < b because M < 2^b and log(2) < 1.  Hence 2H < 2b.
    entropy_nats_strict_upper = support_size.bit_length()
    twice_entropy_nats_strict_upper = 2 * entropy_nats_strict_upper
    exact_mean_obstruction = top_energy > twice_entropy_nats_strict_upper

    individually_obstructing = [
        level
        for level, exponent in enumerate(exponents)
        if moment * Fraction(1 << (2 * exponent), sigma_squared)
        > twice_entropy_nats_strict_upper
    ]

    total_energy_log2 = log2_fraction(total_energy)
    top_energy_log2 = log2_fraction(top_energy)
    twice_entropy_log2 = math.log2(twice_entropy_nats)
    energy_to_entropy_log2_ratio = total_energy_log2 - twice_entropy_log2

    # The top-row-only equality threshold.  Strictly exceeding it is needed
    # even before all remaining rows and concentration terms are included.
    required_sigma_log2 = exponents[0] + 0.5 * (
        math.log2(float(moment)) - twice_entropy_log2
    )
    maximum_weight_log2 = params.error_std_log2 + 0.5 * (
        twice_entropy_log2 - math.log2(float(moment))
    )

    return {
        "parameters": asdict(params),
        "support": {
            "size": str(support_size),
            "entropy_bits": entropy_bits,
            "entropy_nats": entropy_nats,
            "conditioning_penalty_bits": binomial_conditioning_log2_penalty(
                params.degree, params.secret_weight
            ),
            "entropy_nats_strict_upper_integer": entropy_nats_strict_upper,
        },
        "quadratic_moment": {
            "numerator": str(moment.numerator),
            "denominator": str(moment.denominator),
            "decimal": float(moment),
        },
        "gadget": {
            "exponents": exponents,
            "independent_rows": len(exponents),
            "individually_obstructing_zero_based_levels": individually_obstructing,
        },
        "energy": {
            "top_row_lower_bound_numerator": str(top_energy.numerator),
            "top_row_lower_bound_denominator": str(top_energy.denominator),
            "top_row_lower_bound_log2": top_energy_log2,
            "all_quadratic_rows_log2": total_energy_log2,
            "twice_entropy_nats": twice_entropy_nats,
            "twice_entropy_nats_log2": twice_entropy_log2,
            "log2_energy_to_twice_entropy_ratio": energy_to_entropy_log2_ratio,
            "exact_twice_entropy_strict_upper_integer": (
                twice_entropy_nats_strict_upper
            ),
        },
        "thresholds": {
            "maximum_single_gadget_weight_log2": maximum_weight_log2,
            "required_top_row_error_std_log2": required_sigma_log2,
            "required_normalized_error_std_log2": (
                required_sigma_log2 - params.modulus_bits
            ),
        },
        "certificate": {
            "top_row_exceeds_exact_entropy_upper": exact_mean_obstruction,
            "mean_margin_positive": total_energy < twice_entropy_nats,
            "concentration_can_repair": False if exact_mean_obstruction else None,
            "result": (
                "FAIL_EQUAL_COVARIANCE_RENYI_MEAN_CONDITION"
                if exact_mean_obstruction
                else "MEAN_CONDITION_NOT_RULED_OUT"
            ),
        },
    }


def _negacyclic_square(secret: tuple[int, ...]) -> tuple[int, ...]:
    degree = len(secret)
    square = [0] * degree
    for left, left_value in enumerate(secret):
        for right, right_value in enumerate(secret):
            output = left + right
            if output < degree:
                square[output] += left_value * right_value
            else:
                square[output - degree] -= left_value * right_value
    return tuple(square)


def self_test() -> None:
    """Exhaustively check the closed moment on small exact-weight slices."""

    for degree in (4, 6, 8):
        for weight in range(1, min(3, degree) + 1):
            total = 0
            count = 0
            for support in combinations(range(degree), weight):
                for signs in product((-1, 1), repeat=weight):
                    secret = [0] * degree
                    for coordinate, sign in zip(support, signs, strict=True):
                        secret[coordinate] = sign
                    square = _negacyclic_square(tuple(secret))
                    total += sum(coefficient * coefficient for coefficient in square)
                    count += 1
            observed = Fraction(total, count)
            expected = fixed_weight_square_second_moment(degree, weight)
            if observed != expected:
                raise AssertionError(
                    f"moment mismatch at n={degree}, w={weight}: "
                    f"{observed} != {expected}"
                )


def print_report(result: dict[str, Any]) -> None:
    params = result["parameters"]
    support = result["support"]
    moment = result["quadratic_moment"]
    gadget = result["gadget"]
    energy = result["energy"]
    thresholds = result["thresholds"]
    certificate = result["certificate"]

    print("TFHEpp lvl5boot fixed-weight Renyi certificate")
    print(
        "  n={degree} q=2^{modulus_bits} w={secret_weight} "
        "sigma=2^{error_std_log2}".format(**params)
    )
    print(
        "  ordinary quadratic rows={independent_rows}, gadget exponents=".format(
            **gadget
        )
        + ",".join(str(value) for value in gadget["exponents"])
    )
    print("secret support")
    print(f"  log2 |T_n,w| = {support['entropy_bits']:.12f} bits")
    print(
        "  fixed-weight conditioning penalty = "
        f"{support['conditioning_penalty_bits']:.12f} bits"
    )
    print("exact fixed-weight square moment")
    print(
        "  E[||S^2||_2^2] = "
        f"{moment['numerator']}/{moment['denominator']} "
        f"= {moment['decimal']:.12f}"
    )
    print("descriptor-independent quadratic energy")
    print(
        "  log2(top-row lower bound) = "
        f"{energy['top_row_lower_bound_log2']:.12f}"
    )
    print(
        "  log2(all-row quadratic mean) = "
        f"{energy['all_quadratic_rows_log2']:.12f}"
    )
    print(
        "  2*entropy = "
        f"{energy['twice_entropy_nats']:.12f} nats "
        f"(log2={energy['twice_entropy_nats_log2']:.12f})"
    )
    print(
        "  log2(E[F_quad] / (2*entropy)) = "
        f"{energy['log2_energy_to_twice_entropy_ratio']:.12f}"
    )
    print("exact obstruction check")
    print(
        "  2*entropy < "
        f"{energy['exact_twice_entropy_strict_upper_integer']} and the "
        "top-row rational lower bound exceeds it: "
        f"{certificate['top_row_exceeds_exact_entropy_upper']}"
    )
    print(
        "  individually obstructing rows = "
        f"{len(gadget['individually_obstructing_zero_based_levels'])}/"
        f"{gadget['independent_rows']}"
    )
    print("necessary top-row-only thresholds")
    print(
        "  maximum gadget log2 weight at sigma=2^"
        f"{params['error_std_log2']}: "
        f"{thresholds['maximum_single_gadget_weight_log2']:.6f}"
    )
    print(
        "  required log2 error std for gadget 2^"
        f"{gadget['exponents'][0]}: "
        f"> {thresholds['required_top_row_error_std_log2']:.6f}"
    )
    print(
        "  required normalized log2 error std: "
        f"> {thresholds['required_normalized_error_std_log2']:.6f}"
    )
    print(f"RESULT: {certificate['result']}")
    if certificate["top_row_exceeds_exact_entropy_upper"]:
        print(
            "The mean term is already too large; nonnegative concentration "
            "and descriptor terms cannot repair this Renyi certificate."
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=1 << 15)
    parser.add_argument("--secret-weight", type=int, default=96)
    parser.add_argument("--modulus-bits", type=int, default=640)
    parser.add_argument("--error-std-log2", type=int, default=33)
    parser.add_argument("--gadget-base-bits", type=int, default=18)
    parser.add_argument("--gadget-levels", type=int, default=35)
    parser.add_argument("--target-bits", type=int, default=128)
    parser.add_argument("--json", action="store_true", help="emit JSON")
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="exhaustively validate the square-moment formula on small rings",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
    params = Parameters(
        degree=args.degree,
        secret_weight=args.secret_weight,
        modulus_bits=args.modulus_bits,
        error_std_log2=args.error_std_log2,
        gadget_base_bits=args.gadget_base_bits,
        gadget_levels=args.gadget_levels,
        target_bits=args.target_bits,
    )
    result = analyse(params)
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print_report(result)


if __name__ == "__main__":
    main()
