#!/usr/bin/env python3
"""Extract the N=65536 thin-BGV automorphism schedule and frontier bounds.

The formulas mirror Bootstrapping_BGV_BFV's `GetOurHypercubeStructure`,
`GetGeneralBabyGiantParams`, `MatMulGeneralBabyGiantSwitchKeys`, and the
non-cyclic thin-recryption branch.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from typing import Any


def multiplicative_order(value: int, modulus: int) -> int:
    if math.gcd(value, modulus) != 1:
        raise ValueError("value is not a unit")
    current = 1
    for order in range(1, modulus + 1):
        current = current * value % modulus
        if current == 1:
            return order
    raise AssertionError("unreachable")


def index_to_sequence(index: int, sizes: tuple[int, ...]) -> tuple[int, ...]:
    """Magma's one-based mixed-radix IndexToSequence."""
    if index < 1:
        raise ValueError("index must be one based")
    value = index - 1
    result = []
    for size in sizes:
        result.append(value % size)
        value //= size
    if value != 0:
        raise ValueError("index exceeds mixed-radix domain")
    return tuple(result)


def divisors(value: int) -> list[int]:
    return [candidate for candidate in range(1, value + 1) if value % candidate == 0]


def general_baby_giant(sizes: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    total = math.prod(sizes)
    maximum = max(sizes)
    maximum_index = sizes.index(maximum)
    choices = []
    for divisor in divisors(total // maximum):
        baby = math.ceil(math.sqrt(total) / divisor)
        giant = math.ceil(maximum / baby)
        cost = baby * divisor + giant * (total // maximum // divisor)
        choices.append((cost, divisor, baby, giant))
    _, divisor, final_baby, final_giant = min(choices)

    babies: list[int] = []
    giants: list[int] = []
    for index, size in enumerate(sizes):
        if index == maximum_index:
            babies.append(final_baby)
            giants.append(final_giant)
        else:
            factor = math.gcd(size, divisor)
            divisor //= factor
            babies.append(factor)
            giants.append(size // factor)
    return tuple(babies), tuple(giants)


def rotation_exponent(
    generators: tuple[int, ...], rotation: tuple[int, ...], modulus: int
) -> int:
    result = 1
    for generator, exponent in zip(generators, rotation, strict=True):
        result = result * pow(generator, exponent, modulus) % modulus
    return result


def generated_subgroup(generators: set[int], modulus: int) -> set[int]:
    subgroup = {1}
    frontier = [1]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = current * generator % modulus
            if candidate not in subgroup:
                subgroup.add(candidate)
                frontier.append(candidate)
    return subgroup


@dataclass(frozen=True)
class Schedule:
    degree: int
    cyclotomic_index: int
    plaintext_prime: int
    frobenius_order: int
    slot_count: int
    dimension_sizes: tuple[int, ...]
    dimension_generators: tuple[int, ...]
    baby_sizes: tuple[int, ...]
    giant_sizes: tuple[int, ...]
    baby_product: int
    giant_product: int
    switch_key_exponents: frozenset[int]
    generated_group_size: int
    peak_live_ciphertexts: int


@dataclass(frozen=True)
class SwitchUse:
    family: str
    one_based_index: int
    mixed_radix: tuple[int, ...]
    rotation: tuple[int, ...]
    exponent: int


@dataclass(frozen=True)
class GateManifest:
    schedule: Schedule
    baby_switches: tuple[SwitchUse, ...]
    giant_switches: tuple[SwitchUse, ...]
    backward_exponent: int
    trace_exponent: int
    stages: tuple[dict[str, Any], ...]

    def canonical_object(self) -> dict[str, Any]:
        schedule = self.schedule
        return {
            "schema": "tfhepp-compact-bgv-gate-manifest-v1",
            "degree": schedule.degree,
            "cyclotomic_index": schedule.cyclotomic_index,
            "plaintext_prime": schedule.plaintext_prime,
            "hensel_input_exponent": 1,
            "hensel_bootstrap_exponent": 2,
            "frobenius_order": schedule.frobenius_order,
            "slot_count": schedule.slot_count,
            "dimension_sizes": list(schedule.dimension_sizes),
            "dimension_generators": list(schedule.dimension_generators),
            "baby_sizes": list(schedule.baby_sizes),
            "giant_sizes": list(schedule.giant_sizes),
            "baby_product": schedule.baby_product,
            "giant_product": schedule.giant_product,
            "peak_live_ciphertexts": schedule.peak_live_ciphertexts,
            "baby_switches": [
                {
                    "index": use.one_based_index,
                    "mixed_radix": list(use.mixed_radix),
                    "rotation": list(use.rotation),
                    "exponent": use.exponent,
                }
                for use in self.baby_switches
            ],
            "giant_switches": [
                {
                    "index": use.one_based_index,
                    "mixed_radix": list(use.mixed_radix),
                    "rotation": list(use.rotation),
                    "exponent": use.exponent,
                }
                for use in self.giant_switches
            ],
            "backward_exponent": self.backward_exponent,
            "trace_exponent": self.trace_exponent,
            "unique_switch_exponents": sorted(schedule.switch_key_exponents),
            "generated_group_size": schedule.generated_group_size,
            "stages": list(self.stages),
        }

    def canonical_json(self) -> str:
        return json.dumps(
            self.canonical_object(), sort_keys=True, separators=(",", ":")
        )

    def sha256(self) -> str:
        return hashlib.sha256(self.canonical_json().encode("utf-8")).hexdigest()


def thin_power_of_two_schedule(
    degree: int = 65_536, plaintext_prime: int = 65_537
) -> Schedule:
    if degree <= 0 or degree & (degree - 1):
        raise ValueError("degree must be a positive power of two")
    modulus = 2 * degree
    order = multiplicative_order(plaintext_prime % modulus, modulus)
    slot_count = degree // order
    if plaintext_prime % 4 != 1:
        raise ValueError("the targeted non-cyclic branch requires p mod 4 = 1")

    # GetOurHypercubeStructure with one pre-adjustment matrix dimension.
    requested_dimension = modulus // (2 * order)
    dimensions = (requested_dimension // 2, 2)
    generators = (5 % modulus, modulus - 1)

    # The first non-cyclic stage calls MatMul2DBabyGiant with dimensions
    # [GetNbDimensions(), 1], hence order (-1, 5).
    stage_sizes = (dimensions[1], dimensions[0])
    stage_generators = (generators[1], generators[0])
    babies, giants = general_baby_giant(stage_sizes)
    baby_product = math.prod(babies)
    giant_product = math.prod(giants)

    exponents: list[int] = []
    for index in range(2, baby_product + 1):
        exponents.append(
            rotation_exponent(stage_generators, index_to_sequence(index, babies), modulus)
        )
    for index in range(2, giant_product + 1):
        sequence = index_to_sequence(index, giants)
        rotation = tuple(babies[i] * sequence[i] for i in range(len(babies)))
        exponents.append(rotation_exponent(stage_generators, rotation, modulus))

    # Bad-dimension backward key and the one trace key for d=2.
    backward = pow(generators[0], -dimensions[0], modulus)
    trace = pow(plaintext_prime, order // 2, modulus)
    exponents.extend((backward, trace))
    unique_exponents = frozenset(exponents)
    subgroup = generated_subgroup(set(unique_exponents), modulus)

    # MatMulGeneralBadDimensionBabyGiant keeps v and v_prime, plus c,
    # c_prime, tmp, and the output accumulator.
    peak_live = 2 * baby_product + 4
    return Schedule(
        degree=degree,
        cyclotomic_index=modulus,
        plaintext_prime=plaintext_prime,
        frobenius_order=order,
        slot_count=slot_count,
        dimension_sizes=dimensions,
        dimension_generators=generators,
        baby_sizes=babies,
        giant_sizes=giants,
        baby_product=baby_product,
        giant_product=giant_product,
        switch_key_exponents=unique_exponents,
        generated_group_size=len(subgroup),
        peak_live_ciphertexts=peak_live,
    )


def thin_power_of_two_gate_manifest(
    degree: int = 65_536, plaintext_prime: int = 65_537
) -> GateManifest:
    schedule = thin_power_of_two_schedule(degree, plaintext_prime)
    modulus = schedule.cyclotomic_index
    stage_sizes = (schedule.dimension_sizes[1], schedule.dimension_sizes[0])
    stage_generators = (
        schedule.dimension_generators[1],
        schedule.dimension_generators[0],
    )

    baby_switches: list[SwitchUse] = []
    for index in range(2, schedule.baby_product + 1):
        mixed = index_to_sequence(index, schedule.baby_sizes)
        rotation = mixed
        baby_switches.append(
            SwitchUse(
                "baby",
                index,
                mixed,
                rotation,
                rotation_exponent(stage_generators, rotation, modulus),
            )
        )

    giant_switches: list[SwitchUse] = []
    for index in range(2, schedule.giant_product + 1):
        mixed = index_to_sequence(index, schedule.giant_sizes)
        rotation = tuple(
            schedule.baby_sizes[i] * mixed[i] for i in range(len(mixed))
        )
        giant_switches.append(
            SwitchUse(
                "giant",
                index,
                mixed,
                rotation,
                rotation_exponent(stage_generators, rotation, modulus),
            )
        )

    backward = pow(
        schedule.dimension_generators[0],
        -schedule.dimension_sizes[0],
        modulus,
    )
    trace = pow(plaintext_prime, schedule.frobenius_order // 2, modulus)
    stages = (
        {
            "name": "forward_bad_dimension",
            "operation": "matmul_2d_bad_dimension_bsgs",
            "plaintext_hensel_exponent": 1,
            "source_width": 1,
            "peak_width": schedule.peak_live_ciphertexts,
            "target_width": 1,
            "baby_switch_count": len(baby_switches),
            "giant_switch_count": len(giant_switches),
            "uses_backward": True,
        },
        {
            "name": "homomorphic_inner_product",
            "operation": "encrypted_secret_inner_product",
            "plaintext_hensel_exponent": 2,
            "source_width": 1,
            "peak_width": 4,
            "target_width": 1,
        },
        {
            "name": "degree_two_trace",
            "operation": "trace_add_automorphism",
            "plaintext_hensel_exponent": 2,
            "source_width": 1,
            "peak_width": 2,
            "target_width": 1,
            "automorphism": trace,
        },
        {
            "name": "inverse_bad_dimension",
            "operation": "matmul_2d_bad_dimension_bsgs",
            "plaintext_hensel_exponent": 2,
            "source_width": 1,
            "peak_width": schedule.peak_live_ciphertexts,
            "target_width": 1,
            "baby_switch_count": len(baby_switches),
            "giant_switch_count": len(giant_switches),
            "uses_backward": True,
        },
        {
            "name": "bounded_digit_extraction",
            "operation": "lowest_digit_removal_then_exact_division",
            "plaintext_hensel_exponent": 2,
            "source_width": 1,
            "peak_width": 8,
            "target_width": 1,
        },
        {
            "name": "cyclic_return",
            "operation": "frontier_transition_and_rerandomize",
            "plaintext_hensel_exponent": 1,
            "source_width": 1,
            "peak_width": 2,
            "target_width": 1,
        },
    )
    return GateManifest(
        schedule=schedule,
        baby_switches=tuple(baby_switches),
        giant_switches=tuple(giant_switches),
        backward_exponent=backward,
        trace_exponent=trace,
        stages=stages,
    )


def scalar_direct_gate_manifest_object(
    degree: int = 65_536, plaintext_prime: int = 65_537
) -> dict[str, Any]:
    """Canonical manifest for the genuine scalar-only BGV bootstrap."""
    return {
        "schema": "tfhepp-compact-bgv-scalar-bootstrap-manifest-v2",
        "degree": degree,
        "cyclotomic_index": 2 * degree,
        "plaintext_prime": plaintext_prime,
        "plaintext_square": plaintext_prime * plaintext_prime,
        "full_rns_limbs": 23,
        "phase_lift_gadget_digits": 2,
        "trace_gadget_digits": 23,
        "trace_exponents": [
            pow(5, 1 << index, 2 * degree) for index in range(15)
        ]
        + [2 * degree - 1],
        "trace_drop_after": [8, 16],
        "digit_error_bound": 23,
        "digit_polynomial_degree": 93,
        "stages": [
            {
                "name": "modulus_down",
                "operation": "bgv_lsb_drop_to_one_limb",
                "source_limbs": 23,
                "target_limbs": 1,
            },
            {
                "name": "low_to_high_phase_lift",
                "operation": "centered_scale_and_cross_modulus_transition",
                "source_limbs": 1,
                "target_limbs": 23,
                "gadget_digits": 2,
                "key_error_scale": plaintext_prime * plaintext_prime,
            },
            {
                "name": "constant_projection",
                "operation": "normalized_16_stage_galois_trace",
                "source_limbs": 23,
                "target_limbs": 21,
            },
            {
                "name": "bounded_digit_removal",
                "operation": "centered_odd_bsgs_with_level_drops",
                "source_plaintext": plaintext_prime * plaintext_prime,
                "carry_bound": 23,
            },
            {
                "name": "exact_public_division",
                "operation": "multiply_components_by_inverse_plaintext_prime",
                "source_plaintext": plaintext_prime * plaintext_prime,
                "target_plaintext": plaintext_prime,
            },
            {
                "name": "multiplication_closure",
                "operation": "quadratic_hint_multiply_and_level_drop",
                "target_input_limbs": 1,
            },
        ],
    }


def scalar_direct_gate_manifest_json(
    degree: int = 65_536, plaintext_prime: int = 65_537
) -> str:
    return json.dumps(
        scalar_direct_gate_manifest_object(degree, plaintext_prime),
        sort_keys=True,
        separators=(",", ":"),
    )


def scalar_direct_gate_manifest_sha256(
    degree: int = 65_536, plaintext_prime: int = 65_537
) -> str:
    return hashlib.sha256(
        scalar_direct_gate_manifest_json(degree, plaintext_prime).encode()
    ).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=65_536)
    parser.add_argument("--plaintext-prime", type=int, default=65_537)
    parser.add_argument("--rns-limbs", type=int, default=15)
    parser.add_argument(
        "--manifest-json",
        action="store_true",
        help="emit the canonical gate manifest as JSON",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = thin_power_of_two_gate_manifest(args.degree, args.plaintext_prime)
    schedule = manifest.schedule
    if args.manifest_json:
        print(json.dumps(manifest.canonical_object(), sort_keys=True, indent=2))
        print(f"sha256={manifest.sha256()}")
        return 0
    residue_bytes = 8
    full_ciphertext_bytes = (
        2 * schedule.degree * schedule.generated_group_size * residue_bytes * args.rns_limbs
    )
    frontier_bytes = (
        2 * schedule.degree * schedule.peak_live_ciphertexts * residue_bytes * args.rns_limbs
    )
    print("Compact-cover thin-BGV schedule")
    print(
        f"  N={schedule.degree} m={schedule.cyclotomic_index} "
        f"p={schedule.plaintext_prime} d={schedule.frobenius_order} "
        f"slots={schedule.slot_count}"
    )
    print(
        f"  hypercube dimensions={schedule.dimension_sizes} "
        f"generators={schedule.dimension_generators}"
    )
    print(
        f"  baby={schedule.baby_sizes} ({schedule.baby_product}) "
        f"giant={schedule.giant_sizes} ({schedule.giant_product})"
    )
    print(f"  unique switch-key automorphisms={len(schedule.switch_key_exponents)}")
    print(f"  generated subgroup size={schedule.generated_group_size}")
    print(f"  peak live ciphertexts={schedule.peak_live_ciphertexts}")
    print(f"  gate manifest sha256={manifest.sha256()}")
    print(f"  full-cover ciphertext={full_ciphertext_bytes / 2**30:.2f} GiB")
    print(f"  scheduled frontier={frontier_bytes / 2**30:.2f} GiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
