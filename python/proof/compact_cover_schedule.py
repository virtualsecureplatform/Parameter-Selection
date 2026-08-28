#!/usr/bin/env python3
"""Extract the N=65536 thin-BGV automorphism schedule and frontier bounds.

The formulas mirror Bootstrapping_BGV_BFV's `GetOurHypercubeStructure`,
`GetGeneralBabyGiantParams`, `MatMulGeneralBabyGiantSwitchKeys`, and the
non-cyclic thin-recryption branch.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=65_536)
    parser.add_argument("--plaintext-prime", type=int, default=65_537)
    parser.add_argument("--rns-limbs", type=int, default=15)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    schedule = thin_power_of_two_schedule(args.degree, args.plaintext_prime)
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
    print(f"  full-cover ciphertext={full_ciphertext_bytes / 2**30:.2f} GiB")
    print(f"  scheduled frontier={frontier_bytes / 2**30:.2f} GiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
