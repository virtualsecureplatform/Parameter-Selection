#!/usr/bin/env python3
"""Exact theorem-side screen for TFHEpp delayed short preimages.

The current invertible-minor construction first solves over the 16-bit target
ring and lifts the result by the complete 32-to-16 scale.  It therefore cannot
benefit from high-word error scaling.  This script checks that obstruction and
then screens the remaining alternative:

    L A = 2^16 H  (mod 2^32)

with a genuinely short high-modulus row ``L`` which need not be divisible by
``2^16``.  All parameter arithmetic and cardinality comparisons are exact.

For the theorem screen, coefficients are restricted to the canonical nonzero
ternary family: the first nonzero coefficient is positive.  Its exact size is
``(3^m - 1)/2``.  The support-cluster second moment gives the one-target
failure bound

    (Q^d - 1)/N + (2^d - 1) H/N^2,

where ``H = (5^m - 2*3^m + 1)/4``.  All comparisons, including the union bound
over every requested KSK row, use exact integers and rational numbers.  This
proves information-theoretic existence with the reported failure bound.  It
does not construct an efficient public preimage solver; two related gadget
targets expose the separate SIS search barrier.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

from tfhe_subset_joint_screen import read_parameters


def _fraction_json(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": str(float(value)),
    }


def _floor_sqrt_fraction(value: Fraction) -> int:
    if value < 0:
        raise ValueError("square-root input must be nonnegative")
    return math.isqrt(value.numerator // value.denominator)


def _canonical_ternary_count(rows: int) -> int:
    if rows < 0:
        raise ValueError("row count must be nonnegative")
    return (pow(3, rows) - 1) // 2


def _same_support_ordered_pair_count(rows: int) -> int:
    """Exact H for the full-support-budget canonical ternary family."""
    if rows < 0:
        raise ValueError("row count must be nonnegative")
    return (pow(5, rows) - 2 * pow(3, rows) + 1) // 4


def _minimum_canonical_rows_for_count(target: int) -> int:
    """Least m satisfying (3^m - 1)/2 >= target."""
    if target < 0:
        raise ValueError("candidate target must be nonnegative")
    low = 0
    high = max(1, target.bit_length())
    while _canonical_ternary_count(high) < target:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if _canonical_ternary_count(middle) >= target:
            high = middle
        else:
            low = middle + 1
    return low


def _canonical_second_moment_terms(
    rows: int, dimension: int, modulus_bits: int
) -> tuple[Fraction, Fraction]:
    """Return the singleton and same-support terms of the exact failure bound."""
    candidates = _canonical_ternary_count(rows)
    if candidates == 0:
        raise ValueError("the canonical candidate family must be nonempty")
    target_space = 1 << (modulus_bits * dimension)
    clustered_pairs = _same_support_ordered_pair_count(rows)
    singleton_term = Fraction(target_space - 1, candidates)
    cluster_term = Fraction(
        ((1 << dimension) - 1) * clustered_pairs,
        candidates * candidates,
    )
    return singleton_term, cluster_term


def _canonical_second_moment_failure(
    rows: int, dimension: int, modulus_bits: int
) -> Fraction:
    singleton, cluster = _canonical_second_moment_terms(
        rows, dimension, modulus_bits
    )
    return singleton + cluster


def _certified_failure_bits(value: Fraction) -> int:
    """Largest b >= 0 for which value <= 2^-b; return -1 when value > 1."""
    if value < 0:
        raise ValueError("failure bound must be nonnegative")
    if value == 0:
        raise ValueError("zero has no finite negative-log certificate")
    if value > 1:
        return -1
    bits = max(0, value.denominator.bit_length() - value.numerator.bit_length())
    while value.numerator * (1 << bits) > value.denominator:
        bits -= 1
    while value.numerator * (1 << (bits + 1)) <= value.denominator:
        bits += 1
    return bits


def _minimum_rows_for_target_family_failure(
    *,
    dimension: int,
    modulus_bits: int,
    target_count: int,
    security_bits: int,
) -> int:
    """Exact least m whose target-family union bound is at most 2^-security_bits.

    The singleton term supplies a rigorous lower bound on m.  Starting there
    and scanning upward avoids relying on an unproved monotonicity shortcut for
    the complete two-term expression.
    """
    if target_count <= 0:
        raise ValueError("target count must be positive")
    if security_bits < 0:
        raise ValueError("security bits must be nonnegative")
    target_space = 1 << (modulus_bits * dimension)
    necessary_candidates = target_count * (target_space - 1) * (1 << security_bits)
    rows = _minimum_canonical_rows_for_count(necessary_candidates)
    threshold = Fraction(1, 1 << security_bits)
    while target_count * _canonical_second_moment_failure(
        rows, dimension, modulus_bits
    ) > threshold:
        rows += 1
    return rows


def _dimension_screen(
    *,
    name: str,
    dimension: int,
    modulus_bits: int,
    multiplicity_slack_bits: int,
    radius_squared: int,
    output_rows: int,
) -> dict[str, Any]:
    target_space_bits = modulus_bits * dimension
    target_space = 1 << target_space_bits
    one_preimage_rows = _minimum_canonical_rows_for_count(target_space)
    cardinality_slack_rows = _minimum_canonical_rows_for_count(
        target_space << multiplicity_slack_bits
    )
    one_target_rows = _minimum_rows_for_target_family_failure(
        dimension=dimension,
        modulus_bits=modulus_bits,
        target_count=1,
        security_bits=multiplicity_slack_bits,
    )
    target_family_rows = _minimum_rows_for_target_family_failure(
        dimension=dimension,
        modulus_bits=modulus_bits,
        target_count=output_rows,
        security_bits=multiplicity_slack_bits,
    )
    candidate_count = _canonical_ternary_count(target_family_rows)
    previous_failure = output_rows * _canonical_second_moment_failure(
        target_family_rows - 1, dimension, modulus_bits
    )
    singleton_term, cluster_term = _canonical_second_moment_terms(
        target_family_rows, dimension, modulus_bits
    )
    one_target_failure = singleton_term + cluster_term
    target_family_failure = output_rows * one_target_failure
    requested_bound = Fraction(1, 1 << multiplicity_slack_bits)
    worst_case_energy = target_family_rows
    return {
        "name": name,
        "dimension": dimension,
        "target_space_bits": target_space_bits,
        "coefficient_alphabet": [-1, 0, 1],
        "canonical_sign_rule": "first nonzero coefficient is +1",
        "canonical_candidate_count_formula": "(3^rows - 1) / 2",
        "same_support_ordered_pair_count_formula": (
            "(5^rows - 2*3^rows + 1) / 4"
        ),
        "second_moment_failure_formula": (
            "(Q^dimension - 1)/N + (2^dimension - 1)*H/N^2"
        ),
        "minimum_rows_for_target_space_cardinality": one_preimage_rows,
        "multiplicity_slack_bits": multiplicity_slack_bits,
        "minimum_rows_for_requested_candidate_multiplicity": cardinality_slack_rows,
        "minimum_rows_for_one_target_second_moment_failure": one_target_rows,
        "minimum_rows_for_all_ksk_targets_union_bound": target_family_rows,
        "minimum_is_exact": previous_failure > requested_bound,
        "candidate_count_floor_log2": candidate_count.bit_length() - 1,
        "candidate_capacity_target_bits": (
            target_space_bits + multiplicity_slack_bits
        ),
        "one_target_failure_certified_bits": _certified_failure_bits(
            one_target_failure
        ),
        "same_support_term_certified_bits": _certified_failure_bits(cluster_term),
        "all_ksk_targets_union_bound_certified_bits": _certified_failure_bits(
            target_family_failure
        ),
        "all_ksk_targets_meet_requested_failure": (
            target_family_failure <= requested_bound
        ),
        "previous_row_count_fails_requested_failure": (
            previous_failure > requested_bound
        ),
        "worst_case_row_energy": worst_case_energy,
        "available_integer_radius_squared": radius_squared,
        "all_candidates_fit_noise_radius": worst_case_energy <= radius_squared,
        "disjoint_source_rows_for_all_ksk_outputs": (
            output_rows * target_family_rows
        ),
        "disjoint_block_gram_bound": "L L^T <= rows_per_block * I",
        "result": (
            "PASS_INFORMATION_THEORETIC_EXISTENCE_AND_NORM_SOLVER_OPEN"
            if (
                worst_case_energy <= radius_squared
                and target_family_failure <= requested_bound
            )
            else "FAIL_NOISE_RADIUS"
        ),
    }


def analyse(
    root: Path,
    *,
    multiplicity_slack_bits: int = 128,
) -> dict[str, Any]:
    if multiplicity_slack_bits < 0:
        raise ValueError("multiplicity slack must be nonnegative")
    parameters, evidence = read_parameters(root)

    scale_bits = parameters.lvl1_torus_bits - parameters.lvl0_torus_bits
    scale = 1 << scale_bits
    source_sigma = Fraction(
        1 << (parameters.lvl1_torus_bits + parameters.lvl1_alpha_log2),
        1,
    )
    source_variance = source_sigma * source_sigma
    target_sigma = parameters.lvl0_alpha * (1 << parameters.lvl0_torus_bits)
    target_variance = target_sigma * target_sigma

    high_modulus_energy_budget = (
        Fraction(scale * scale, 1) * target_variance / source_variance
    )
    integer_radius = _floor_sqrt_fraction(high_modulus_energy_budget)
    integer_radius_squared = integer_radius * integer_radius

    lifted_minor_minimum_norm = scale
    lifted_minor_minimum_target_ring_energy = 1
    lifted_minor_minimum_derived_variance = source_variance

    dimensions = [
        _dimension_screen(
            name="suffix-only",
            dimension=parameters.suffix_dimension,
            modulus_bits=parameters.lvl1_torus_bits,
            multiplicity_slack_bits=multiplicity_slack_bits,
            radius_squared=integer_radius_squared,
            output_rows=parameters.ksk_rows,
        ),
        _dimension_screen(
            name="full-secret",
            dimension=parameters.lvl1_dimension,
            modulus_bits=parameters.lvl1_torus_bits,
            multiplicity_slack_bits=multiplicity_slack_bits,
            radius_squared=integer_radius_squared,
            output_rows=parameters.ksk_rows,
        ),
    ]

    return {
        "scope": {
            "statement": (
                "Exact canonical-ternary second-moment and noise screen.  It proves "
                "random-matrix existence at the reported union-bound failure but does "
                "not construct an efficient public short-preimage algorithm."
            ),
            "tfhepp_root": str(root.resolve()),
            "source_binding_hashes": evidence["sha256"],
        },
        "parameters": {
            "large_modulus_bits": parameters.lvl1_torus_bits,
            "target_modulus_bits": parameters.lvl0_torus_bits,
            "scale_bits": scale_bits,
            "scale": scale,
            "suffix_dimension": parameters.suffix_dimension,
            "full_secret_dimension": parameters.lvl1_dimension,
            "ksk_rows": parameters.ksk_rows,
            "source_sigma": _fraction_json(source_sigma),
            "source_variance": _fraction_json(source_variance),
            "target_sigma": _fraction_json(target_sigma),
            "target_variance": _fraction_json(target_variance),
        },
        "lifted_invertible_minor": {
            "minimum_nonzero_high_modulus_norm": lifted_minor_minimum_norm,
            "allowed_high_modulus_integer_radius": integer_radius,
            "minimum_nonzero_target_ring_energy": (
                lifted_minor_minimum_target_ring_energy
            ),
            "minimum_derived_variance": _fraction_json(
                lifted_minor_minimum_derived_variance
            ),
            "target_variance": _fraction_json(target_variance),
            "result": (
                "FAIL_SOURCE_VARIANCE_ALREADY_EXCEEDS_TARGET"
                if lifted_minor_minimum_derived_variance > target_variance
                else "PASS"
            ),
        },
        "genuinely_short_high_modulus": {
            "exact_energy_budget": _fraction_json(high_modulus_energy_budget),
            "maximum_integer_euclidean_radius": integer_radius,
            "radius_squared": integer_radius_squared,
            "next_integer_radius_exceeds_budget": (
                (integer_radius + 1) ** 2 > high_modulus_energy_budget
            ),
            "dimension_screens": dimensions,
            "result": (
                "PASSES_EXISTENCE_GRAM_AND_NOISE_SCREEN_SOLVER_SIS_BARRIER"
                if all(
                    item["all_candidates_fit_noise_radius"]
                    and item["all_ksk_targets_meet_requested_failure"]
                    for item in dimensions
                )
                else "FAILS_NECESSARY_SCREEN"
            ),
        },
        "remaining_boundary": {
            "random_matrix_high_probability_existence": "PROVED_BY_SECOND_MOMENT",
            "disjoint_block_joint_gram": "PROVED_INFORMATION_THEORETICALLY",
            "public_polynomial_time_solver": "OPEN_WITH_GEOMETRIC_TARGET_SIS_BARRIER",
            "full_secret_joint_brk_ksk_reduction": "PROVED_CONDITIONALLY_ON_SOLVER_AND_ERRORS",
            "concrete_cpp_error_law_bridge": "OPEN_EXACT_SAMPLER_COMPARISON",
        },
    }


def self_test(root: Path) -> None:
    report = analyse(root)
    parameters = report["parameters"]
    minor = report["lifted_invertible_minor"]
    short = report["genuinely_short_high_modulus"]
    suffix, full = short["dimension_screens"]

    assert parameters["suffix_dimension"] == 394
    assert parameters["full_secret_dimension"] == 1024
    assert parameters["ksk_rows"] == 5516
    assert parameters["scale"] == 1 << 16
    assert parameters["source_sigma"]["numerator"] == "128"
    assert parameters["source_variance"]["numerator"] == "16384"
    assert minor["result"] == "FAIL_SOURCE_VARIANCE_ALREADY_EXCEEDS_TARGET"
    assert short["maximum_integer_euclidean_radius"] == 3104
    assert short["next_integer_radius_exceeds_budget"]

    assert suffix["minimum_rows_for_target_space_cardinality"] == 7956
    assert suffix["minimum_rows_for_requested_candidate_multiplicity"] == 8037
    assert suffix["minimum_rows_for_one_target_second_moment_failure"] == 8037
    assert suffix["minimum_rows_for_all_ksk_targets_union_bound"] == 8044
    assert full["minimum_rows_for_target_space_cardinality"] == 20675
    assert full["minimum_rows_for_requested_candidate_multiplicity"] == 20756
    assert full["minimum_rows_for_one_target_second_moment_failure"] == 20756
    assert full["minimum_rows_for_all_ksk_targets_union_bound"] == 20764
    assert suffix["minimum_is_exact"] and full["minimum_is_exact"]
    assert suffix["all_ksk_targets_meet_requested_failure"]
    assert full["all_ksk_targets_meet_requested_failure"]
    assert suffix["previous_row_count_fails_requested_failure"]
    assert full["previous_row_count_fails_requested_failure"]
    assert suffix["same_support_term_certified_bits"] > 6000
    assert full["same_support_term_certified_bits"] > 16000
    assert suffix["all_candidates_fit_noise_radius"]
    assert full["all_candidates_fit_noise_radius"]


def _print_human(report: dict[str, Any]) -> None:
    parameters = report["parameters"]
    minor = report["lifted_invertible_minor"]
    short = report["genuinely_short_high_modulus"]
    print("TFHEpp delayed short-preimage second-moment screen")
    print(
        "  suffix / full dimension / KSK rows: "
        f"{parameters['suffix_dimension']} / {parameters['full_secret_dimension']} / "
        f"{parameters['ksk_rows']}"
    )
    print(f"  exact high-modulus integer radius: {short['maximum_integer_euclidean_radius']}")
    print(f"  lifted invertible-minor route: {minor['result']}")
    for item in short["dimension_screens"]:
        print(
            f"  {item['name']} ternary candidates: "
            f"one-target rows={item['minimum_rows_for_one_target_second_moment_failure']}, "
            f"all-target rows={item['minimum_rows_for_all_ksk_targets_union_bound']}, "
            f"energy={item['worst_case_row_energy']}, result={item['result']}"
        )
    print(f"  conclusion: {short['result']}")


def main() -> None:
    default_root = Path(__file__).resolve().parents[3] / "TFHEpp"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tfhepp-root", type=Path, default=default_root)
    parser.add_argument("--multiplicity-slack-bits", type=int, default=128)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test(args.tfhepp_root)
    report = analyse(
        args.tfhepp_root,
        multiplicity_slack_bits=args.multiplicity_slack_bits,
    )
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        _print_human(report)
        if args.self_test:
            print("  self-test: PASS")


if __name__ == "__main__":
    main()
