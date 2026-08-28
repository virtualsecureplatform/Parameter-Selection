#!/usr/bin/env python3
"""Structural parameter screen for the regular-cover BGV construction.

This script evaluates proof-induced costs that are absent from a native BGV
noise estimator: pivot unit density, cover expansion, the dense-binary
ring-Regev masking term, ordinary-RLWE sample count, and the union-bound loss
across cover components.  Native bootstrap noise is supplied as an external
per-component bound because the regular cover does not change that recurrence.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ScreenInput:
    degree: int
    group_size: int
    q_bits: int
    target_bits: float = 128.0
    pivot_attempts: int = 128
    evaluation_rows: int = 0
    public_key_rows: int | None = None
    native_failure_bits: float = math.inf
    native_input_noise_bits: float | None = None
    native_fresh_noise_bits: float | None = None


@dataclass(frozen=True)
class ScreenResult:
    cover_log2_cardinality: int
    minimum_public_key_rows: int
    public_key_rows: int
    leftover_log2_bound: float
    unit_probability: float
    pivot_failure_log2: float
    base_rlwe_samples: int
    ciphertext_storage_bits: int
    cover_failure_bits: float
    pivot_condition_ok: bool
    leftover_ok: bool
    failure_ok: bool
    refresh_contract_ok: bool | None


def minimum_public_key_rows(cover_log2_cardinality: int, target_bits: float) -> int:
    """Rows needed for sqrt(|Cover|^2 / 2^m)/2 <= 2^-target."""
    return max(0, math.ceil(2.0 * (cover_log2_cardinality + target_bits - 1.0)))


def screen(candidate: ScreenInput) -> ScreenResult:
    if candidate.degree <= 0 or candidate.group_size <= 0 or candidate.q_bits <= 0:
        raise ValueError("degree, group_size, and q_bits must be positive")
    if candidate.pivot_attempts <= 0 or candidate.evaluation_rows < 0:
        raise ValueError("pivot_attempts must be positive and evaluation_rows nonnegative")

    cover_dimension = candidate.degree * candidate.group_size
    cover_log_card = cover_dimension * candidate.q_bits
    min_rows = minimum_public_key_rows(cover_log_card, candidate.target_bits)
    pk_rows = candidate.public_key_rows if candidate.public_key_rows is not None else min_rows
    if pk_rows < 0:
        raise ValueError("public_key_rows must be nonnegative")

    leftover_log2 = -1.0 + cover_log_card - pk_rows / 2.0

    # Use logs so very large q and dimensions do not underflow prematurely.
    q = 2**candidate.q_bits
    log_unit_probability = cover_dimension * math.log1p(-1.0 / q)
    unit_probability = math.exp(log_unit_probability)
    failure_probability = -math.expm1(log_unit_probability)
    pivot_failure_log2 = (
        -math.inf
        if failure_probability == 0.0
        else candidate.pivot_attempts * math.log2(failure_probability)
    )

    cover_failure_bits = (
        math.inf
        if math.isinf(candidate.native_failure_bits)
        else candidate.native_failure_bits - math.log2(candidate.group_size)
    )

    refresh_contract_ok: bool | None
    if candidate.native_input_noise_bits is None or candidate.native_fresh_noise_bits is None:
        refresh_contract_ok = None
    else:
        refresh_contract_ok = (
            candidate.native_fresh_noise_bits < candidate.native_input_noise_bits
        )

    total_cover_rows = pk_rows + candidate.evaluation_rows + candidate.pivot_attempts
    return ScreenResult(
        cover_log2_cardinality=cover_log_card,
        minimum_public_key_rows=min_rows,
        public_key_rows=pk_rows,
        leftover_log2_bound=leftover_log2,
        unit_probability=unit_probability,
        pivot_failure_log2=pivot_failure_log2,
        base_rlwe_samples=candidate.group_size * total_cover_rows,
        ciphertext_storage_bits=2 * cover_log_card,
        cover_failure_bits=cover_failure_bits,
        pivot_condition_ok=(2 * cover_dimension <= q),
        leftover_ok=(leftover_log2 <= -candidate.target_bits),
        failure_ok=(
            not math.isinf(candidate.native_failure_bits)
            and cover_failure_bits >= candidate.target_bits
        ),
        refresh_contract_ok=refresh_contract_ok,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=512)
    parser.add_argument("--group-size", type=int, default=512)
    parser.add_argument("--qbits", type=int, default=900)
    parser.add_argument("--target-bits", type=float, default=128.0)
    parser.add_argument("--pivot-attempts", type=int, default=128)
    parser.add_argument("--evaluation-rows", type=int, default=0)
    parser.add_argument("--public-key-rows", type=int)
    parser.add_argument("--native-failure-bits", type=float, default=math.inf)
    parser.add_argument("--native-input-noise-bits", type=float)
    parser.add_argument("--native-fresh-noise-bits", type=float)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = screen(
        ScreenInput(
            degree=args.degree,
            group_size=args.group_size,
            q_bits=args.qbits,
            target_bits=args.target_bits,
            pivot_attempts=args.pivot_attempts,
            evaluation_rows=args.evaluation_rows,
            public_key_rows=args.public_key_rows,
            native_failure_bits=args.native_failure_bits,
            native_input_noise_bits=args.native_input_noise_bits,
            native_fresh_noise_bits=args.native_fresh_noise_bits,
        )
    )

    print("Regular-cover BGV proof screen")
    print(f"  log2|Cover|={result.cover_log2_cardinality}")
    print(
        f"  dense-binary public-key rows={result.public_key_rows} "
        f"(minimum {result.minimum_public_key_rows})"
    )
    print(
        f"  leftover log2(distance)<={result.leftover_log2_bound:.2f} "
        f"{'PASS' if result.leftover_ok else 'FAIL'}"
    )
    print(f"  pivot unit probability={result.unit_probability:.12g}")
    print(
        f"  pivot log2(failure)<={result.pivot_failure_log2:.2f} "
        f"{'PASS' if result.pivot_condition_ok else 'UNPROVED'}"
    )
    print(f"  ordinary base-RLWE samples={result.base_rlwe_samples}")
    print(f"  ciphertext storage={result.ciphertext_storage_bits} bits")
    print(
        f"  cover failure security={result.cover_failure_bits:.2f} bits "
        f"{'PASS' if result.failure_ok else 'FAIL/UNSPECIFIED'}"
    )
    if result.refresh_contract_ok is not None:
        print(f"  refresh contraction={'PASS' if result.refresh_contract_ok else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
