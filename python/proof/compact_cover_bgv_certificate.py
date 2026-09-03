#!/usr/bin/env python3
"""Deterministic certificate for the N=65536 scalar Binary-NTT BGV cycle."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import cache
import hashlib
import json
import math
from pathlib import Path
import sys

_HERE = Path(__file__).resolve().parent
_PYTHON = _HERE.parent
sys.path[:0] = [str(_HERE), str(_PYTHON)]

from compact_cover_schedule import scalar_direct_gate_manifest_sha256  # noqa: E402
from noiseestimation.bfv_polynomial import (  # noqa: E402
    lowest_digit_removal_polynomial_over_range,
)


DEGREE = 65_536
PLAINTEXT_PRIME = 65_537
PLAINTEXT_SQUARE = PLAINTEXT_PRIME**2
SECRET_WEIGHT = 32
CBD_ETA = 20
PHASE_LIFT_DIGITS = 2
TRACE_DIGITS = 23
TRACE_KEY_COUNT = 16
ACTIVE_RNS_LIMBS = 20
TRACE_EXPONENTS = (
    5, 25, 625, 128_481, 28_609, 61_313, 7_937, 81_409,
    31_745, 63_489, 126_977, 122_881, 114_689, 98_305, 65_537, 131_071,
)
DIGIT_ERROR_BOUND = 23
SECURITY_TARGET_BITS = 128.0
EVALUATION_SECURITY_PROXY_BITS = 162.94
CONTEXT_SECURITY_PROXY_BITS = 155.93

RNS_PRIMES_AND_ROOTS = (
    (2_301_972_608_560_791_553, 5),
    (2_295_217_002_959_732_737, 5),
    (2_291_839_200_159_203_329, 7),
    (2_280_016_890_357_350_401, 7),
    (2_274_950_186_156_556_289, 7),
    (2_271_009_416_222_605_313, 3),
    (2_265_942_712_021_811_201, 3),
    (2_252_994_467_953_115_137, 7),
    (2_244_549_960_951_791_617, 7),
    (2_230_475_782_616_252_417, 3),
    (2_227_097_979_815_723_009, 3),
    (2_217_527_538_547_556_353, 5),
    (2_203_453_360_212_017_153, 3),
    (2_179_808_740_608_311_297, 5),
    (2_156_164_121_004_605_441, 3),
    (2_152_786_318_204_076_033, 3),
    (2_124_637_961_532_997_633, 11),
    (2_114_504_553_131_409_409, 14),
    (2_109_437_848_930_615_297, 5),
    (2_102_682_243_329_556_481, 13),
    (2_078_474_656_592_429_057, 3),
    (2_065_526_412_523_732_993, 11),
    (2_057_081_905_522_409_473, 5),
)


@dataclass(frozen=True)
class ErrorState:
    limbs: int
    bound: int


def modulus(limbs: int) -> int:
    return math.prod(value for value, _ in RNS_PRIMES_AND_ROOTS[:limbs])


def modulus_drop(state: ErrorState, target_limbs: int) -> ErrorState:
    if not 0 < target_limbs <= state.limbs:
        raise ValueError("invalid target limb count")
    error = state.bound
    for active in range(state.limbs, target_limbs, -1):
        dropped = RNS_PRIMES_AND_ROOTS[active - 1][0]
        # Exact BGV LSB division plus one body and h signed mask roundings.
        error = (error + dropped - 1) // dropped + (SECRET_WEIGHT + 2) // 2
    return ErrorState(target_limbs, error)


def key_switch_error(limbs: int, rows: int, gadget_bits: int | None = None) -> int:
    if gadget_bits is None:
        gadget_bits = math.ceil(modulus(limbs).bit_length() / rows)
    digit_bound = 1 << (gadget_bits - 1)
    return rows * DEGREE * digit_bound * CBD_ETA


def add(left: ErrorState, right: ErrorState) -> ErrorState:
    limbs = min(left.limbs, right.limbs)
    return ErrorState(
        limbs,
        modulus_drop(left, limbs).bound + modulus_drop(right, limbs).bound,
    )


def multiply_and_drop(left: ErrorState, right: ErrorState) -> ErrorState:
    limbs = min(left.limbs, right.limbs)
    if limbs < 2:
        raise ValueError("multiplication exhausted RNS levels")
    left_error = modulus_drop(left, limbs).bound
    right_error = modulus_drop(right, limbs).bound
    message_bound = PLAINTEXT_SQUARE // 2
    raw = (
        message_bound * (left_error + right_error)
        + PLAINTEXT_SQUARE * left_error * right_error
        + (message_bound * message_bound + PLAINTEXT_SQUARE - 1)
        // PLAINTEXT_SQUARE
        + 1
    )
    return modulus_drop(ErrorState(limbs, raw), limbs - 1)


def scale(state: ErrorState, scalar: int) -> ErrorState:
    scalar %= PLAINTEXT_SQUARE
    if scalar > PLAINTEXT_SQUARE // 2:
        scalar -= PLAINTEXT_SQUARE
    return ErrorState(state.limbs, abs(scalar) * state.bound)


def generic_bsgs_error(coefficients: list[int], value: ErrorState) -> ErrorState:
    baby_count = 3
    giant_count = 6
    baby = [ErrorState(value.limbs, 0), value]
    baby.append(multiply_and_drop(baby[1], baby[1]))
    baby.append(multiply_and_drop(baby[2], baby[1]))
    giant = [baby[baby_count]]
    for _ in range(1, giant_count):
        giant.append(multiply_and_drop(giant[-1], giant[-1]))

    def evaluate(offset: int, length: int, level: int) -> ErrorState:
        accumulator = ErrorState(value.limbs, 0)
        if length <= 1:
            return accumulator
        if level == 0:
            for power in range(1, min(length - 1, baby_count) + 1):
                coefficient = coefficients[offset + power]
                if coefficient:
                    accumulator = add(accumulator, scale(baby[power], coefficient))
            return accumulator
        split = min(length, baby_count * (1 << (level - 1)))
        lower = evaluate(offset, split, level - 1)
        if split == length:
            return lower
        upper = evaluate(offset + split, length - split, level - 1)
        return add(lower, multiply_and_drop(upper, giant[level - 1]))

    return evaluate(0, len(coefficients), giant_count)


def digit_polynomial_error(coefficients: list[int], value: ErrorState) -> ErrorState:
    odd = coefficients[0] == 0 and all(
        coefficients[index] == 0 for index in range(2, len(coefficients), 2)
    )
    if not odd:
        return generic_bsgs_error(coefficients, value)
    square = multiply_and_drop(value, value)
    inner = generic_bsgs_error(coefficients[1::2], square)
    return multiply_and_drop(inner, value)


@dataclass(frozen=True)
class Certificate:
    schema: str
    gate_manifest_sha256: str
    degree: int
    plaintext_prime: int
    plaintext_square: int
    secret_weight: int
    operational_secret_distribution: str
    error_distribution: str
    rns_primes: tuple[int, ...]
    primitive_roots: tuple[int, ...]
    modulus_bits: int
    modulus_log2: float
    low_modulus_bits: int
    phase_lift_digits: int
    trace_digits: int
    trace_key_count: int
    trace_exponents: tuple[int, ...]
    trace_drop_after: tuple[int, ...]
    digit_error_bound: int
    digit_polynomial_degree: int
    digit_polynomial_coefficients: tuple[int, ...]
    accepted_input_error_bound: int
    accepted_input_error_log2: float
    projected_error_bound: int
    projected_error_log2_bound: float
    output_limbs: int
    output_error_bound: int
    output_error_log2_bound: float
    output_capacity: int
    output_capacity_log2: float
    multiplication_input_limbs: int
    addition_output_error_bound: int
    multiplication_output_error_bound: int
    multiplication_output_error_log2_bound: float
    multiplication_capacity: int
    multiplication_capacity_log2: float
    cycle_contraction_bits: float
    bootstrap_key_bytes: int
    evaluation_security_proxy_bits: float
    context_security_proxy_bits: float
    unit_pivot_success_probability: float
    unit_pivot_failure_log2: float
    conditioned_context_security_proxy_bits: float
    combined_security_proxy_bits: float
    correctness_failure_log2_bound: float
    correctness_certified: bool
    estimated_security_meets_target: bool

    def canonical_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))

    def sha256(self) -> str:
        return hashlib.sha256(self.canonical_json().encode()).hexdigest()


@cache
def build_certificate() -> Certificate:
    selected = RNS_PRIMES_AND_ROOTS[:ACTIVE_RNS_LIMBS]
    primes = tuple(value for value, _ in selected)
    roots = tuple(root for _, root in selected)
    congruence = math.lcm(2 * DEGREE, PLAINTEXT_SQUARE)
    if any((prime - 1) % congruence for prime in primes):
        raise AssertionError("invalid BGV/NTT prime congruence")
    full_modulus = modulus(len(primes))
    low_modulus = primes[0]

    # The first term of the BGV carry is bounded by (h+1)/2.  Reserve six
    # additional units for the scaled old error, leaving strict room below 23.
    accepted_input_error = 6 * low_modulus // PLAINTEXT_SQUARE
    carry_bound_twice = SECRET_WEIGHT + 1 + 12
    if carry_bound_twice >= 2 * DIGIT_ERROR_BOUND:
        raise AssertionError("carry interval is not strict")

    low_gadget_bits = math.ceil(low_modulus.bit_length() / PHASE_LIFT_DIGITS)
    quotient_factor = (low_modulus - 1) // PLAINTEXT_SQUARE
    phase_error = (
        accepted_input_error
        + quotient_factor * DIGIT_ERROR_BOUND
        + key_switch_error(len(primes), PHASE_LIFT_DIGITS, low_gadget_bits)
    )
    state = ErrorState(len(primes), phase_error)
    for stage in range(TRACE_KEY_COUNT):
        state = ErrorState(
            state.limbs,
            2 * state.bound
            + key_switch_error(state.limbs, TRACE_DIGITS),
        )
        if stage in (7, 15):
            state = modulus_drop(state, state.limbs - 1)
    # 65536^(-1) mod 65537^2 has centered representative -65538.
    state = scale(state, PLAINTEXT_SQUARE - 65_538)
    projected_error = state.bound

    polynomial = lowest_digit_removal_polynomial_over_range(
        PLAINTEXT_PRIME, DIGIT_ERROR_BOUND
    )
    output = digit_polynomial_error(polynomial, state)
    output_capacity = modulus(output.limbs) // (2 * PLAINTEXT_PRIME)

    multiplication_input = modulus_drop(output, 2)
    multiplied = multiply_and_drop(multiplication_input, multiplication_input)
    multiplication_capacity = low_modulus // (2 * PLAINTEXT_PRIME)
    addition_input = modulus_drop(output, 1)
    added = add(addition_input, addition_input)

    unit_pivot_log_success = sum(
        DEGREE * math.log1p(-1.0 / prime) for prime in primes
    )
    unit_pivot_success = math.exp(unit_pivot_log_success)
    unit_pivot_failure = -math.expm1(unit_pivot_log_success)
    conditioned_context_security = (
        CONTEXT_SECURITY_PROXY_BITS + math.log2(unit_pivot_success)
    )
    # The P0--P2 full-source hybrid contains the evaluation table and one
    # context row; the P2--P1 auxiliary-leakage triangle charges the context
    # row once more.  Keep both context terms explicit rather than silently
    # treating the heterogeneous full problem as the evaluation-row problem.
    combined_security = -math.log2(
        2.0 ** (-EVALUATION_SECURITY_PROXY_BITS)
        + 2.0 * 2.0 ** (-conditioned_context_security)
    )
    cycle_contraction = math.log2(accepted_input_error) - math.log2(multiplied.bound)
    ciphertext_bytes = 2 * DEGREE * 8
    phase_key_bytes = PHASE_LIFT_DIGITS * len(primes) * ciphertext_bytes + 48
    trace_key_bytes = sum(
        TRACE_DIGITS * (len(primes) if index < 8 else len(primes) - 1)
        * ciphertext_bytes + 48
        for index in range(TRACE_KEY_COUNT)
    )
    hint_bytes = 4 * DEGREE * len(primes) * 8
    key_bytes = phase_key_bytes + trace_key_bytes + hint_bytes

    correctness_certified = (
        output.bound < output_capacity
        and multiplied.bound < multiplication_capacity
        and added.bound <= accepted_input_error
        and cycle_contraction > 0
        and len(polynomial) - 1 == 93
    )
    estimated_security_meets_target = combined_security >= SECURITY_TARGET_BITS
    return Certificate(
        schema="tfhepp-compact-bgv-scalar-certificate-v6",
        gate_manifest_sha256=scalar_direct_gate_manifest_sha256(),
        degree=DEGREE,
        plaintext_prime=PLAINTEXT_PRIME,
        plaintext_square=PLAINTEXT_SQUARE,
        secret_weight=SECRET_WEIGHT,
        operational_secret_distribution=f"fixed-weight-{SECRET_WEIGHT}-signed-ternary",
        error_distribution=f"CBD({CBD_ETA})",
        rns_primes=primes,
        primitive_roots=roots,
        modulus_bits=full_modulus.bit_length(),
        modulus_log2=math.log2(full_modulus),
        low_modulus_bits=low_modulus.bit_length(),
        phase_lift_digits=PHASE_LIFT_DIGITS,
        trace_digits=TRACE_DIGITS,
        trace_key_count=TRACE_KEY_COUNT,
        trace_exponents=TRACE_EXPONENTS,
        trace_drop_after=(8, 16),
        digit_error_bound=DIGIT_ERROR_BOUND,
        digit_polynomial_degree=len(polynomial) - 1,
        digit_polynomial_coefficients=tuple(polynomial),
        accepted_input_error_bound=accepted_input_error,
        accepted_input_error_log2=math.log2(accepted_input_error),
        projected_error_bound=projected_error,
        projected_error_log2_bound=math.log2(projected_error),
        output_limbs=output.limbs,
        output_error_bound=output.bound,
        output_error_log2_bound=math.log2(output.bound),
        output_capacity=output_capacity,
        output_capacity_log2=math.log2(output_capacity),
        multiplication_input_limbs=2,
        addition_output_error_bound=added.bound,
        multiplication_output_error_bound=multiplied.bound,
        multiplication_output_error_log2_bound=math.log2(multiplied.bound),
        multiplication_capacity=multiplication_capacity,
        multiplication_capacity_log2=math.log2(multiplication_capacity),
        cycle_contraction_bits=cycle_contraction,
        bootstrap_key_bytes=key_bytes,
        evaluation_security_proxy_bits=EVALUATION_SECURITY_PROXY_BITS,
        context_security_proxy_bits=CONTEXT_SECURITY_PROXY_BITS,
        unit_pivot_success_probability=unit_pivot_success,
        unit_pivot_failure_log2=math.log2(unit_pivot_failure),
        conditioned_context_security_proxy_bits=conditioned_context_security,
        combined_security_proxy_bits=combined_security,
        correctness_failure_log2_bound=-1_000_000.0,
        correctness_certified=correctness_certified,
        estimated_security_meets_target=estimated_security_meets_target,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.json:
        print(json.dumps(asdict(certificate), sort_keys=True, indent=2))
    else:
        print("Binary-NTT scalar BGV cycle certificate")
        print(
            f"  Q={certificate.modulus_bits} bits, q0={certificate.low_modulus_bits} bits"
        )
        print(
            f"  projected error <2^{certificate.projected_error_log2_bound:.2f}"
        )
        print(
            f"  output: {certificate.output_limbs} limbs, error "
            f"<2^{certificate.output_error_log2_bound:.2f}, capacity "
            f"2^{certificate.output_capacity_log2:.2f}"
        )
        print(
            f"  multiply/drop error <2^{certificate.multiplication_output_error_log2_bound:.2f}, "
            f"capacity 2^{certificate.multiplication_capacity_log2:.2f}"
        )
        print(
            f"  cycle contraction={certificate.cycle_contraction_bits:.2f} bits, "
            f"combined security proxy={certificate.combined_security_proxy_bits:.2f} bits"
        )
        print(f"  evaluation key={certificate.bootstrap_key_bytes / 2**30:.2f} GiB")
        print(f"  manifest={certificate.gate_manifest_sha256}")
        print(f"  sha256={certificate.sha256()}")
        print(
            "  correctness="
            f"{'CERTIFIED' if certificate.correctness_certified else 'FAIL'}"
        )
        print(
            "  estimated security="
            f"{'MEETS TARGET' if certificate.estimated_security_meets_target else 'FAIL'}"
        )
    return 0 if (certificate.correctness_certified and
                  certificate.estimated_security_meets_target) else 1


if __name__ == "__main__":
    raise SystemExit(main())
