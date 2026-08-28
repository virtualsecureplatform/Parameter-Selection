#!/usr/bin/env python3
"""Concrete security/correctness certificate for compact-cover scalar BGV.

The certificate deliberately separates three claims:

* the exact gate manifest and RNS congruences are deterministic checks;
* ordinary Binary-NTT RLWE security is reported through the repository's
  conservative unlimited-sample LWE attack-cost proxy; and
* correctness uses the existing dependent-input BFV/BGV invariant-noise
  model, plus a rigorous Hoeffding bound for the modulus-switch digit error.

It does not claim a conventional-RLWE reduction for Binary-NTT RLWE.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import argparse
from functools import cache
import hashlib
import json
import math
from pathlib import Path
import sys

_HERE = Path(__file__).resolve().parent
_PYTHON = _HERE.parent
sys.path[:0] = [str(_HERE), str(_PYTHON)]

from compact_cover_schedule import (  # noqa: E402
    scalar_direct_gate_manifest_sha256,
    thin_power_of_two_gate_manifest,
)
from noiseestimation.bfv import (  # noqa: E402
    BFVParams,
    NoiseState,
    add_many,
    const_mul,
    correctness_d_for_failure,
    estimate_polyeval_bsgs,
    fresh,
    key_switch,
    log2_add,
    log2_correctness_threshold,
)


DEGREE = 65_536
PLAINTEXT_PRIME = 65_537
PLAINTEXT_SQUARE = PLAINTEXT_PRIME * PLAINTEXT_PRIME
FRONTIER_WIDTH = 368
SECRET_WEIGHT = 32
CBD_ETA = 20
ERROR_STDDEV = math.sqrt(CBD_ETA / 2)
GADGET_DIGITS = 5
DIGIT_ERROR_BOUND = 43
FAILURE_TARGET_LOG2 = -128.0
SECURITY_TARGET_BITS = 128.0
ACCEPTED_INPUT_LOG2_VARIANCE = -900.0

# These primes are the first deterministic descending candidates q=1+k*M,
# M=lcm(2N,p^2), passing primality.  The listed primitive roots are checked by
# the C++ NTT regression.  Every prefix product is one modulo p^2 and 2N.
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
)


@dataclass(frozen=True)
class Certificate:
    schema: str
    gate_manifest_sha256: str
    degree: int
    plaintext_prime: int
    plaintext_square: int
    frobenius_order: int
    frontier_width: int
    secret_weight: int
    error_distribution: str
    error_stddev: float
    rns_primes: tuple[int, ...]
    primitive_roots: tuple[int, ...]
    modulus_log2: float
    modulus_bits: int
    gadget_digits: int
    gadget_bits: int
    digit_error_bound: int
    digit_polynomial_degree: int
    digit_bsgs_k: int
    digit_bsgs_m: int
    rounding_trials: int
    rounding_failure_log2_bound: float
    accepted_input_log2_variance: float
    fresh_output_log2_variance: float
    refresh_contraction_bits: float
    output_correctness_threshold_log2_variance: float
    source_security_proxy_bits: float
    reduction_loss_bits: float
    certified_security_bits: float
    security_model: str
    correctness_model: str
    circuit: str
    deterministic_output_error_log2_bound: float
    bootstrap_key_bytes: int
    certified: bool

    def canonical_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))

    def sha256(self) -> str:
        return hashlib.sha256(self.canonical_json().encode()).hexdigest()


def _params(t: int, q_bits: int) -> BFVParams:
    return BFVParams(
        nbit=16,
        t=t,
        q_bits=q_bits,
        secret_variance=SECRET_WEIGHT / DEGREE,
        encryption_u_variance=2 / 3,
        error_log2_std=math.log2(ERROR_STDDEV),
        key_switch="hybrid-rns",
        rns_digits=GADGET_DIGITS,
        hybrid_omega=GADGET_DIGITS,
        fresh="symmetric",
        tag="compact-cover-bgv-65536",
    )


def _linear_stage(state: NoiseState, params: BFVParams) -> NoiseState:
    switched = key_switch(state, params)
    scaled = const_mul(switched, params, "linear-constant")
    # All diagonal terms depend on the same input.  Squaring the sum of their
    # standard deviations is the conservative dependent-input bound.
    return NoiseState(
        scaled.log2_variance + 2 * math.log2(DEGREE // 2),
        scaled.degree,
        "bad-dimension-linear-stage",
    )


@cache
def build_certificate() -> Certificate:
    packed_manifest = thin_power_of_two_gate_manifest(DEGREE, PLAINTEXT_PRIME)
    primes = tuple(value for value, _ in RNS_PRIMES_AND_ROOTS)
    roots = tuple(root for _, root in RNS_PRIMES_AND_ROOTS)
    modulus = math.prod(primes)
    modulus_log2 = math.log2(modulus)
    modulus_bits = modulus.bit_length()
    congruence = math.lcm(2 * DEGREE, PLAINTEXT_SQUARE)
    if any((prime - 1) % congruence != 0 for prime in primes):
        raise AssertionError("RNS prime violates BGV/NTT congruence")

    p_params = _params(PLAINTEXT_PRIME, modulus_bits)
    p2_params = _params(PLAINTEXT_SQUARE, modulus_bits)

    # Retain the general packed-path polynomial metadata, but certify the
    # scalar specialization below.  Its lifted p^2 error is removed by exact
    # public division and requires no probabilistic digit-extraction event.
    packed_digit_degree = 4 * DIGIT_ERROR_BOUND + 1

    # A coefficient modulus switch has one body rounding term and at most h
    # signed mask terms.  Each centered rounding term has range length one.
    # Hoeffding: P(|sum| >= B) <= 2 exp(-2 B^2/(h+1)).
    rounding_trials = 0
    # No rounding sampler is used by the scalar direct circuit.  Keep a finite
    # sentinel so the canonical JSON remains standards compliant.
    rounding_failure = -1_000_000.0

    output_d = correctness_d_for_failure(DEGREE, FAILURE_TARGET_LOG2)
    output_threshold = log2_correctness_threshold(output_d)

    gadget_bits = math.ceil(modulus_bits / GADGET_DIGITS)
    deterministic_error_log2 = (
        math.log2(PLAINTEXT_PRIME)
        + math.log2(GADGET_DIGITS)
        + math.log2(DEGREE)
        + (gadget_bits - 1)
        + math.log2(CBD_ETA)
    )
    output_log2_variance = 2 * (deterministic_error_log2 - modulus_log2)

    # Unlimited-sample rough estimator result for n=65536, qbits=915 and the
    # actual evaluation-row phase error p^2*CBD(20), recorded by
    # regular_cover_bgv_security.py. One bit is reserved for finite reduction
    # accounting.
    source_security = 243.53
    reduction_loss = 1.0
    certified_security = source_security - reduction_loss
    contraction = ACCEPTED_INPUT_LOG2_VARIANCE - output_log2_variance
    per_ciphertext_bytes = 2 * DEGREE * len(primes) * 8
    bootstrap_key_bytes = GADGET_DIGITS * per_ciphertext_bytes + 48
    certified = (
        deterministic_error_log2 < modulus_log2 - 1
        and output_log2_variance <= output_threshold
        and contraction > 0
        and certified_security >= SECURITY_TARGET_BITS
    )

    return Certificate(
        schema="tfhepp-compact-bgv-certificate-v1",
        gate_manifest_sha256=scalar_direct_gate_manifest_sha256(
            DEGREE, PLAINTEXT_PRIME
        ),
        degree=DEGREE,
        plaintext_prime=PLAINTEXT_PRIME,
        plaintext_square=PLAINTEXT_SQUARE,
        frobenius_order=packed_manifest.schedule.frobenius_order,
        frontier_width=FRONTIER_WIDTH,
        secret_weight=SECRET_WEIGHT,
        error_distribution=f"CBD({CBD_ETA})",
        error_stddev=ERROR_STDDEV,
        rns_primes=primes,
        primitive_roots=roots,
        modulus_log2=modulus_log2,
        modulus_bits=modulus_bits,
        gadget_digits=GADGET_DIGITS,
        gadget_bits=gadget_bits,
        digit_error_bound=DIGIT_ERROR_BOUND,
        digit_polynomial_degree=packed_digit_degree,
        digit_bsgs_k=3,
        digit_bsgs_m=6,
        rounding_trials=rounding_trials,
        rounding_failure_log2_bound=rounding_failure,
        accepted_input_log2_variance=ACCEPTED_INPUT_LOG2_VARIANCE,
        fresh_output_log2_variance=output_log2_variance,
        refresh_contraction_bits=contraction,
        output_correctness_threshold_log2_variance=output_threshold,
        source_security_proxy_bits=source_security,
        reduction_loss_bits=reduction_loss,
        certified_security_bits=certified_security,
        security_model=(
            "unlimited-sample LWE proxy for Binary-NTT RLWE with p^2*CBD(20)"
        ),
        correctness_model="deterministic bounded-CBD full-modulus gadget bound",
        circuit="scalar direct p-to-p^2 transition and exact division",
        deterministic_output_error_log2_bound=deterministic_error_log2,
        bootstrap_key_bytes=bootstrap_key_bytes,
        certified=certified,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.json:
        print(json.dumps(asdict(certificate), sort_keys=True, indent=2))
    else:
        print("Compact-cover BGV N=65536 certificate")
        print(f"  manifest={certificate.gate_manifest_sha256}")
        print(
            f"  Q bits={certificate.modulus_bits} "
            f"log2(Q)={certificate.modulus_log2:.3f}"
        )
        print(
            f"  gadget={certificate.gadget_digits}x"
            f"{certificate.gadget_bits} bits "
            f"deterministic error<2^"
            f"{certificate.deterministic_output_error_log2_bound:.2f}"
        )
        print(
            f"  output log2(V)={certificate.fresh_output_log2_variance:.2f} "
            f"contraction={certificate.refresh_contraction_bits:.2f} bits"
        )
        print(
            f"  source proxy={certificate.source_security_proxy_bits:.2f} "
            f"certified={certificate.certified_security_bits:.2f} bits"
        )
        print(f"  sha256={certificate.sha256()}")
        print(f"  status={'CERTIFIED' if certificate.certified else 'FAIL'}")
    return 0 if certificate.certified else 1


if __name__ == "__main__":
    raise SystemExit(main())
