#!/usr/bin/env python3
"""Exact necessary-condition screen for TFHEpp delayed short preimages.

The current invertible-minor construction first solves over the 16-bit target
ring and lifts the result by the complete 32-to-16 scale.  It therefore cannot
benefit from high-word error scaling.  This script checks that obstruction and
then screens the remaining alternative:

    L A = 2^16 H  (mod 2^32)

with a genuinely short high-modulus row ``L`` which need not be divisible by
``2^16``.  All parameter arithmetic and cardinality comparisons are exact.

For the candidate-capacity screen, coefficients are restricted to
``{-1, 0, 1}``.  The script finds the least row count for which the complete
candidate family has at least the target-space cardinality, both without
slack and with a requested multiplicity slack.  Passing this test proves only
that the elementary cardinality and worst-case Euclidean-norm obstructions do
not reject the family.  It is not a lower bound on the probability that a
random matrix has a preimage, and it does not construct an efficient solver.
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


def _minimum_ternary_rows(target_bits: int) -> int:
    """Least m satisfying 3^m - 1 >= 2^target_bits, using integers only."""
    if target_bits < 0:
        raise ValueError("target bits must be nonnegative")
    target = 1 << target_bits
    low = 0
    high = max(1, target_bits)
    while pow(3, high) - 1 < target:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if pow(3, middle) - 1 >= target:
            high = middle
        else:
            low = middle + 1
    return low


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
    one_preimage_rows = _minimum_ternary_rows(target_space_bits)
    slack_rows = _minimum_ternary_rows(
        target_space_bits + multiplicity_slack_bits
    )
    candidate_count = pow(3, slack_rows) - 1
    previous_candidate_count = pow(3, slack_rows - 1) - 1
    capacity_target = 1 << (target_space_bits + multiplicity_slack_bits)
    worst_case_energy = slack_rows
    return {
        "name": name,
        "dimension": dimension,
        "target_space_bits": target_space_bits,
        "coefficient_alphabet": [-1, 0, 1],
        "nonzero_candidate_count_formula": "3^rows - 1",
        "minimum_rows_for_target_space_cardinality": one_preimage_rows,
        "multiplicity_slack_bits": multiplicity_slack_bits,
        "minimum_rows_for_requested_candidate_multiplicity": slack_rows,
        "minimum_is_exact": (
            candidate_count >= capacity_target
            and previous_candidate_count < capacity_target
        ),
        "candidate_count_floor_log2": candidate_count.bit_length() - 1,
        "candidate_capacity_target_bits": (
            target_space_bits + multiplicity_slack_bits
        ),
        "worst_case_row_energy": worst_case_energy,
        "available_integer_radius_squared": radius_squared,
        "all_candidates_fit_noise_radius": worst_case_energy <= radius_squared,
        "disjoint_source_rows_for_all_ksk_outputs": output_rows * slack_rows,
        "result": (
            "PASS_NECESSARY_CARDINALITY_AND_NORM_ONLY"
            if worst_case_energy <= radius_squared
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
                "Exact necessary cardinality/noise screen; neither a random-matrix "
                "existence theorem nor an efficient short-preimage algorithm."
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
                "PASSES_NECESSARY_CAPACITY_AND_NOISE_SCREEN_SOLVER_OPEN"
                if all(item["all_candidates_fit_noise_radius"] for item in dimensions)
                else "FAILS_NECESSARY_SCREEN"
            ),
        },
        "remaining_boundary": {
            "random_matrix_high_probability_existence": "OPEN_RESEARCH",
            "public_polynomial_time_solver": "OPEN_RESEARCH_INHOMOGENEOUS_SIS",
            "full_secret_joint_brk_ksk_reduction": "OPEN_RESEARCH",
            "concrete_cpp_error_law_bridge": "OPEN_IMPLEMENTATION_ANALYSIS",
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

    assert suffix["minimum_rows_for_target_space_cardinality"] == 7955
    assert suffix["minimum_rows_for_requested_candidate_multiplicity"] == 8036
    assert full["minimum_rows_for_target_space_cardinality"] == 20675
    assert full["minimum_rows_for_requested_candidate_multiplicity"] == 20756
    assert suffix["minimum_is_exact"] and full["minimum_is_exact"]
    assert suffix["all_candidates_fit_noise_radius"]
    assert full["all_candidates_fit_noise_radius"]


def _print_human(report: dict[str, Any]) -> None:
    parameters = report["parameters"]
    minor = report["lifted_invertible_minor"]
    short = report["genuinely_short_high_modulus"]
    print("TFHEpp delayed short-preimage necessary-condition screen")
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
            f"rows={item['minimum_rows_for_requested_candidate_multiplicity']}, "
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
