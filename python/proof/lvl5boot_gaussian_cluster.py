#!/usr/bin/env python3
"""Screen a symmetry-reduced Gaussian cluster for TFHEpp ``lvl5bootparam``.

The candidate cloud contains every signed weight-``w`` secret obtained after
replacing at most two support positions.  For a fixed centre its three orbit
sizes are

    2^w * choose(w, k) * choose(n-w, k),  k = 0, 1, 2.

This is already larger than ``2^128`` at the production parameters, but it is
still represented by three integers rather than by enumerating its members.

For every nonidentity candidate ``t``, one coefficient of

    (t-s) * a + g * (t^2-s^2)

is uniform either on all of ``Z/(2^640)`` or on one parity coset.  The ordinary
top HalfTRGSW row therefore has expected Gaussian kernel overlap below
``2^-604``.  The only analytic input is the elementary theta-sum estimate

    sum_{z in Z} exp(-z^2/(2 sigma^2)) < 4 sigma,  sigma >= 1.

All orbit counts, rational union/Markov bounds, and implementation parameters
are exact.  Floating-point values are presentation-only.  The companion Lean
module kernel-checks the integer/rational certificate and the bounded signed
noise support; it deliberately leaves the displayed theta estimate as an
explicit analytic boundary.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import re
from typing import Any


@dataclass(frozen=True)
class SourceParameters:
    degree: int
    secret_weight: int
    modulus_bits: int
    error_std_log2: int
    gadget_base_bits: int
    gadget_levels: int
    auxiliary_base_bits: int
    auxiliary_levels: int
    signed_noise_bits: int = 64
    target_bits: int = 128

    def validate(self) -> None:
        if self.degree <= 0 or self.degree & (self.degree - 1):
            raise ValueError("degree must be a positive power of two")
        if not 0 < self.secret_weight <= self.degree:
            raise ValueError("secret weight must lie in [1, degree]")
        if self.modulus_bits <= 0 or self.error_std_log2 < 0:
            raise ValueError("modulus and error exponents must be nonnegative")
        if self.error_std_log2 >= self.modulus_bits:
            raise ValueError("the coefficient error must be smaller than q")
        if self.gadget_levels * self.gadget_base_bits > self.modulus_bits:
            raise ValueError("ordinary gadget levels exceed the modulus")
        if self.auxiliary_levels * self.auxiliary_base_bits != self.modulus_bits:
            raise ValueError("auxiliary decomposition does not exactly cover q")
        if self.signed_noise_bits <= 0 or self.target_bits <= 0:
            raise ValueError("noise and target bit counts must be positive")


@dataclass(frozen=True)
class SourceBinding:
    parameters: SourceParameters
    parameter_header: str
    parameter_header_sha256: str
    sampler_header: str
    sampler_header_sha256: str
    modular_gaussian: bool
    saturating_signed_lift: bool
    signed_lift_into_multilimb: bool


def _struct_body(text: str, name: str) -> str:
    match = re.search(rf"\bstruct\s+{re.escape(name)}\s*\{{", text)
    if match is None:
        raise ValueError(f"could not find struct {name}")
    end = text.find("\n};", match.end())
    if end < 0:
        raise ValueError(f"could not find end of struct {name}")
    return text[match.end() : end]


def _decimal_field(body: str, name: str) -> int:
    match = re.search(rf"\b{re.escape(name)}\s*=\s*(\d+)\s*;", body)
    if match is None:
        raise ValueError(f"could not parse {name}")
    return int(match.group(1))


def read_tfhepp_sources(tfhepp_root: Path) -> SourceBinding:
    parameter_header = tfhepp_root / "include" / "params" / "128bit.hpp"
    sampler_header = tfhepp_root / "include" / "utils.hpp"
    parameter_bytes = parameter_header.read_bytes()
    sampler_bytes = sampler_header.read_bytes()
    parameter_text = parameter_bytes.decode("utf-8")
    sampler_text = sampler_bytes.decode("utf-8")
    body = _struct_body(parameter_text, "lvl5bootparam")

    nbit = _decimal_field(body, "nbit")
    limbs_match = re.search(r"using\s+T\s*=\s*MultiLimbUInt<(\d+)>\s*;", body)
    alpha_match = re.search(
        r"\bα\s*=\s*std::pow\(2\.0,\s*(-?\d+)\s*\)\s*;", body
    )
    if limbs_match is None or alpha_match is None:
        raise ValueError("could not parse the multi-limb type or alpha exponent")
    modulus_bits = 64 * int(limbs_match.group(1))
    alpha_log2 = int(alpha_match.group(1))
    parameters = SourceParameters(
        degree=1 << nbit,
        secret_weight=_decimal_field(body, "bfv_key_hamming_weight"),
        modulus_bits=modulus_bits,
        error_std_log2=modulus_bits + alpha_log2,
        gadget_base_bits=_decimal_field(body, "Bgbit"),
        gadget_levels=_decimal_field(body, "l"),
        auxiliary_base_bits=_decimal_field(body, "B̅gbit"),
        auxiliary_levels=_decimal_field(body, "l̅"),
    )
    parameters.validate()

    modular_gaussian = "ErrorDistribution::ModularGaussian" in body
    saturating_signed_lift = (
        "doubleToI64Saturating(distribution(generator) * scale)" in sampler_text
    )
    signed_lift_into_multilimb = (
        "return center + T::from_signed_i64(ival);" in sampler_text
    )
    if not (modular_gaussian and saturating_signed_lift and signed_lift_into_multilimb):
        raise ValueError(
            "TFHEpp source no longer matches the checked modular-Gaussian signed-lift path"
        )

    return SourceBinding(
        parameters=parameters,
        parameter_header=str(parameter_header.resolve()),
        parameter_header_sha256=hashlib.sha256(parameter_bytes).hexdigest(),
        sampler_header=str(sampler_header.resolve()),
        sampler_header_sha256=hashlib.sha256(sampler_bytes).hexdigest(),
        modular_gaussian=modular_gaussian,
        saturating_signed_lift=saturating_signed_lift,
        signed_lift_into_multilimb=signed_lift_into_multilimb,
    )


def orbit_size(degree: int, weight: int, replacements: int) -> int:
    """Number of candidates in the support-replacement orbit ``k``."""

    if not 0 <= replacements <= min(weight, degree - weight):
        return 0
    return (
        (1 << weight)
        * math.comb(weight, replacements)
        * math.comb(degree - weight, replacements)
    )


def radius_orbit_size(degree: int, weight: int, radius: int) -> int:
    return sum(orbit_size(degree, weight, k) for k in range(radius + 1))


def log2_int(value: int) -> float:
    if value <= 0:
        raise ValueError("logarithm requires a positive integer")
    return math.log2(value)


def analyse(binding: SourceBinding) -> dict[str, Any]:
    params = binding.parameters
    params.validate()
    orbit_sizes = [orbit_size(params.degree, params.secret_weight, k) for k in range(3)]
    cloud_size = sum(orbit_sizes)
    nonidentity = cloud_size - 1
    radius_one = sum(orbit_sizes[:2])

    # A parity-coset average has denominator q/2.  Bounding its numerator by
    # the full integer theta sum < 4*sigma gives 8*sigma/q.
    gaussian_pair_overlap_bits = (
        params.modulus_bits - params.error_std_log2 - 3
    )
    expected_excess_denominator_bits = gaussian_pair_overlap_bits
    expected_excess_log2 = (
        log2_int(nonidentity) - expected_excess_denominator_bits
    )

    markov_threshold_bits = params.target_bits
    markov_denominator_bits = (
        gaussian_pair_overlap_bits - markov_threshold_bits
    )
    markov_bad_log2 = log2_int(nonidentity) - markov_denominator_bits

    # The signed-lift support is [-2^(b-1), 2^(b-1)-1], whose difference
    # diameter is below 2^b.  In either a full or parity uniform coset, the
    # chance that a centred shift has magnitude at most 2^b is below
    # 2^(b+2-qbits).
    close_shift_probability_bits = (
        params.modulus_bits - params.signed_noise_bits - 2
    )
    disjointness_bad_log2 = (
        log2_int(nonidentity) - close_shift_probability_bits
    )

    cloud_bits = log2_int(cloud_size)
    min_perfect_r = params.target_bits / (cloud_bits - params.target_bits)
    epsilon = math.ldexp(1.0, -params.target_bits)
    good_effective_security_limit = math.log1p(epsilon) / math.log(2.0)
    orders = []
    for r in (1, 4, 8, 14, 16, 64, 256, 1024):
        exponent = r / (1 + r)
        orders.append(
            {
                "r": r,
                "renyi_order": 1 + r,
                "exponent": exponent,
                "perfect_overlap_security_bits": exponent * cloud_bits,
                "good_mask_security_bits_ceiling": (
                    exponent * good_effective_security_limit
                ),
            }
        )

    return {
        "source_binding": {
            **asdict(binding),
            "parameters": asdict(params),
        },
        "ordinary_top_row": {
            "gadget_exponent": params.modulus_bits - params.gadget_base_bits,
            "coefficient_error_std_exponent": params.error_std_log2,
            "mean_difference": "(t-s)*a + g*(t^2-s^2)",
            "uniform_coordinate_reason": (
                "a coefficient of t-s is odd when supports differ; otherwise "
                "all coefficients are even and (t-s)/2 has a unit coefficient"
            ),
        },
        "symmetry_reduced_cloud": {
            "description": "all signed weight-w secrets replacing at most two support positions",
            "orbit_formula": "2^w * choose(w,k) * choose(n-w,k)",
            "orbit_sizes_k_0_1_2": [str(value) for value in orbit_sizes],
            "radius_one_size": str(radius_one),
            "radius_one_log2": log2_int(radius_one),
            "radius_two_size": str(cloud_size),
            "radius_two_nonidentity": str(nonidentity),
            "radius_two_log2": cloud_bits,
            "radius_one_has_128_bit_capacity": radius_one > (1 << params.target_bits),
            "radius_two_has_128_bit_capacity": cloud_size > (1 << params.target_bits),
            "minimum_r_for_128_bits_if_every_component_overlapped": min_perfect_r,
        },
        "continuous_gaussian_screen": {
            "analytic_boundary": (
                "sum_Z exp(-z^2/(2*sigma^2)) < 4*sigma for sigma >= 1"
            ),
            "per_nonidentity_kernel_overlap_strict_upper": (
                f"2^-{gaussian_pair_overlap_bits}"
            ),
            "expected_raw_effective_size_excess_upper_numerator": str(nonidentity),
            "expected_raw_effective_size_excess_upper_denominator": (
                f"2^{expected_excess_denominator_bits}"
            ),
            "expected_raw_effective_size_excess_upper_log2": expected_excess_log2,
            "markov_threshold": f"2^-{markov_threshold_bits}",
            "markov_bad_probability_upper_numerator": str(nonidentity),
            "markov_bad_probability_upper_denominator": (
                f"2^{markov_denominator_bits}"
            ),
            "markov_bad_probability_upper_log2": markov_bad_log2,
            "coarse_lean_checked_bad_probability": "< 2^-338",
            "good_mask_raw_effective_size_upper": "1 + 2^-128",
            "barycenter_factor_can_increase_effective_size": False,
            "order_optimization": orders,
            "asymptotic_good_mask_security_bits_ceiling": (
                good_effective_security_limit
            ),
            "result": "FAIL_RADIUS_TWO_GAUSSIAN_CLUSTER_CERTIFICATE",
        },
        "implemented_channel": {
            "noise_generation": (
                "normal_distribution<double>, scale 2^33, total saturating "
                "signed-int64 lift, then MultiLimbUInt::from_signed_i64"
            ),
            "signed_lift_support": "[-2^63, 2^63-1]",
            "noise_difference_magnitude_strict_upper": "2^64",
            "per_pair_close_shift_probability_strict_upper": (
                f"2^-{close_shift_probability_bits}"
            ),
            "all_neighbors_close_shift_union_upper_numerator": str(nonidentity),
            "all_neighbors_close_shift_union_upper_denominator": (
                f"2^{close_shift_probability_bits}"
            ),
            "all_neighbors_close_shift_union_upper_log2": disjointness_bad_log2,
            "coarse_lean_checked_bad_probability": "< 2^-436",
            "good_mask_conclusion": (
                "the centre component and every nonidentity radius-two component "
                "have disjoint support in at least one top-row coefficient"
            ),
            "result": "NATURAL_CLUSTER_DISJOINT_ON_GOOD_MASK",
        },
        "certificate": {
            "source_matches_checked_path": True,
            "cloud_large_enough_under_perfect_overlap": (
                cloud_size > (1 << params.target_bits)
            ),
            "continuous_cluster_proves_128_bit_lossiness": False,
            "implemented_channel_improves_overlap": False,
            "result": "NO_128_BIT_CERTIFICATE_FROM_TESTED_RADIUS_TWO_CLUSTER",
            "scope": (
                "This rejects the support-radius-two Gaussian/Gibbs cluster. "
                "It is not an insecurity result and does not exclude nonlocal, "
                "mask-conditioned, or computational reductions."
            ),
        },
    }


def _negacyclic_mul(left: tuple[int, ...], right: tuple[int, ...], q: int) -> tuple[int, ...]:
    degree = len(left)
    output = [0] * degree
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            index = i + j
            sign = 1
            if index >= degree:
                index -= degree
                sign = -1
            output[index] = (output[index] + sign * left_value * right_value) % q
    return tuple(output)


def self_test(binding: SourceBinding) -> None:
    """Check orbit counts and the full/parity coordinate claims on toy rings."""

    if binding.parameters != SourceParameters(
        degree=1 << 15,
        secret_weight=96,
        modulus_bits=640,
        error_std_log2=33,
        gadget_base_bits=18,
        gadget_levels=35,
        auxiliary_base_bits=16,
        auxiliary_levels=40,
    ):
        raise AssertionError("parsed lvl5bootparam no longer matches the proof constants")

    degree, weight = 6, 2
    centre_support = frozenset(range(weight))
    observed = Counter()
    for support_bits in product((0, 1), repeat=degree):
        if sum(support_bits) != weight:
            continue
        support = frozenset(i for i, bit in enumerate(support_bits) if bit)
        replacements = weight - len(centre_support & support)
        for _signs in product((-1, 1), repeat=weight):
            observed[replacements] += 1
    for replacements in range(3):
        expected = orbit_size(degree, weight, replacements)
        if observed[replacements] != expected:
            raise AssertionError(
                f"orbit mismatch for k={replacements}: "
                f"{observed[replacements]} != {expected}"
            )

    q = 8
    centre = (1, 1, 0, 0)
    candidates = ((-1, 1, 0, 0), (0, 1, 1, 0))
    expected_supports = ({0, 2, 4, 6}, set(range(q)))
    for candidate, expected_support in zip(candidates, expected_supports, strict=True):
        difference = tuple(t - s for s, t in zip(centre, candidate, strict=True))
        counts = Counter(
            _negacyclic_mul(difference, mask, q)[0]
            for mask in product(range(q), repeat=len(centre))
        )
        if set(counts) != expected_support or len(set(counts.values())) != 1:
            raise AssertionError("toy mask coordinate is not uniform on its expected coset")


def print_report(result: dict[str, Any]) -> None:
    params = result["source_binding"]["parameters"]
    cloud = result["symmetry_reduced_cloud"]
    gaussian = result["continuous_gaussian_screen"]
    actual = result["implemented_channel"]
    certificate = result["certificate"]

    print("TFHEpp lvl5boot Gaussian-cluster screen")
    print(
        "  source-bound n={degree} q=2^{modulus_bits} w={secret_weight} "
        "sigma=2^{error_std_log2}".format(**params)
    )
    print("symmetry-reduced support-radius-two cloud")
    print("  orbit sizes k=0,1,2: " + ", ".join(cloud["orbit_sizes_k_0_1_2"]))
    print(
        f"  radius one: log2 M={cloud['radius_one_log2']:.12f}; "
        f"128-bit capacity={cloud['radius_one_has_128_bit_capacity']}"
    )
    print(
        f"  radius two: log2 M={cloud['radius_two_log2']:.12f}; "
        f"128-bit capacity={cloud['radius_two_has_128_bit_capacity']}"
    )
    print(
        "  perfect-overlap minimum r for 128 bits: "
        f"{cloud['minimum_r_for_128_bits_if_every_component_overlapped']:.6f}"
    )
    print("continuous Gaussian top-row overlap")
    print(
        "  per nonidentity component < "
        f"{gaussian['per_nonidentity_kernel_overlap_strict_upper']}"
    )
    print(
        "  E[K_raw-1] log2 upper = "
        f"{gaussian['expected_raw_effective_size_excess_upper_log2']:.12f}"
    )
    print(
        "  Pr[K_raw > 1+2^-128] log2 upper = "
        f"{gaussian['markov_bad_probability_upper_log2']:.12f} "
        f"({gaussian['coarse_lean_checked_bad_probability']})"
    )
    print(
        "  best good-mask security-bit ceiling over all r = "
        f"{gaussian['asymptotic_good_mask_security_bits_ceiling']:.6e}"
    )
    print(f"  RESULT: {gaussian['result']}")
    print("actual TFHEpp finite-support channel")
    print(f"  signed error support: {actual['signed_lift_support']}")
    print(
        "  probability some radius-two neighbor remains within one noise "
        "diameter: log2 < "
        f"{actual['all_neighbors_close_shift_union_upper_log2']:.12f} "
        f"({actual['coarse_lean_checked_bad_probability']})"
    )
    print(f"  good-mask conclusion: {actual['good_mask_conclusion']}")
    print(f"  RESULT: {actual['result']}")
    print(f"OVERALL: {certificate['result']}")
    print(certificate["scope"])


def parse_args() -> argparse.Namespace:
    workspace = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tfhepp-root",
        type=Path,
        default=workspace / "TFHEpp",
        help="path to the TFHEpp checkout (default: sibling checkout)",
    )
    parser.add_argument("--json", action="store_true", help="emit JSON")
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="validate source parsing, orbit counts, and toy mask cosets",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    binding = read_tfhepp_sources(args.tfhepp_root)
    if args.self_test:
        self_test(binding)
    result = analyse(binding)
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print_report(result)


if __name__ == "__main__":
    main()
