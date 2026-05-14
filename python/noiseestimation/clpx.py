#!/usr/bin/env python3
"""Noise estimates for TFHEpp's CLPX scheme-switching helpers.

The estimates compose the existing TFHE noise formulas in keyvariation.py along
the operation sequences implemented in TFHEpp/include/bfv-clpx.hpp.  They are
screening estimates: programmable bootstrapping is treated as refreshing the
output encryption noise, while the largest internal PBS-input variance is
reported separately because it is the main semantic-correctness risk.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from scipy.special import erfc

try:
    from python.noiseestimation.keyvariation import (  # noqa: E402
        annihilatecalc,
        brnoisecalc,
        extpnoisecalc,
        iksnoisecalc,
    )
    from python.noiseestimation.params import clpx as default_params  # noqa: E402
except ModuleNotFoundError as exc:
    if exc.name != "python":
        raise
    from noiseestimation.keyvariation import (  # type: ignore[no-redef] # noqa: E402
        annihilatecalc,
        brnoisecalc,
        extpnoisecalc,
        iksnoisecalc,
    )
    from noiseestimation.params import clpx as default_params  # type: ignore[no-redef] # noqa: E402


NEG_INF = float("-inf")


@dataclass(frozen=True)
class TLWES2CLPXEstimate:
    input_variance: float
    iks_variance: float
    pbs_variance: float
    max_packed_input_variance: float
    packed_variance: float
    clpx_value_variance: float
    temp_tlwe_count: int
    max_terms_per_temp: int
    log2_failure: float


@dataclass(frozen=True)
class CLPX2TLWESEstimate:
    input_variance: float
    pbs02_variance: float
    pbs01_variance: float
    sumpra_variance: float
    max_homdecomp_sum_variance: float
    max_internal_pbs_input_variance: float
    max_final_pbs_input_variance: float
    output_variance: float
    produced_tlwes: int
    log2_failure: float


@dataclass(frozen=True)
class CLPXMultiplicationStep:
    multiplication_count: int
    coefficient_variance: float
    digit_value_variance: float
    margin_bits: float
    log2_failure: float
    status: str


@dataclass(frozen=True)
class CLPXMultiplicationEstimate:
    initial_coefficient_variance: float
    initial_digit_value_variance: float
    relin_coefficient_variance: float
    fft_round_coefficient_variance: float
    message_bound: float
    chain: str
    d: float
    steps: list[CLPXMultiplicationStep]
    supported_multiplications: int


@dataclass(frozen=True)
class HomDecompEstimate:
    sums: list[float]
    max_mid_to_low_input_variance: float
    max_sub_pbs_input_variance: float


def q_bits(P) -> int:
    return int(P.q).bit_length() - 1


def log2_variance(variance: float) -> float:
    if variance <= 0:
        return NEG_INF
    return math.log2(float(variance))


def format_log2(variance: float) -> str:
    value = log2_variance(variance)
    if value == NEG_INF:
        return "-inf"
    return f"{value:.2f}"


def _log2_probability(probability: float) -> float:
    if probability <= 0:
        return NEG_INF
    return math.log2(min(1.0, probability))


def _failure_log2(threshold: float, variance: float, count: int = 1) -> float:
    if variance <= 0:
        return NEG_INF
    one = float(erfc(threshold / math.sqrt(2.0 * variance)))
    return _log2_probability(count * one)


def tlwe_failure_log2(P, variance: float, count: int = 1) -> float:
    threshold = P.q / (2.0 * P.plain_modulus)
    return _failure_log2(threshold, variance, count)


def clpx_digit_failure_log2(P, coefficient_variance: float, count: int = 1) -> tuple[float, float]:
    value_variance = (1.0 + float(P.plain_modulus) ** 2) * coefficient_variance
    threshold = P.q / 4.0
    return value_variance, _failure_log2(threshold, value_variance, count)


def clpx_value_margin_log2(P, value_variance: float, d: float = 6.0) -> float:
    if value_variance <= 0:
        return float("inf")
    threshold_log2 = 2.0 * math.log2(P.q / 4.0) - 1.0 - 2.0 * math.log2(d)
    return threshold_log2 - math.log2(value_variance)


def pbs_input_margin_log2(brP, input_variance: float, num_out: int = 1) -> float:
    if input_variance <= 0:
        return float("inf")
    bitwidth = max(0, int(num_out - 1).bit_length())
    half_bin_log2 = q_bits(brP.domainP) - 2 - brP.targetP.nbit + bitwidth
    return 2.0 * half_bin_log2 - math.log2(input_variance)


def scale_variance_between(variance: float, domainP, targetP) -> float:
    if variance <= 0:
        return 0.0
    shift = q_bits(targetP) - q_bits(domainP)
    return math.ldexp(float(variance), 2 * shift)


def shift_variance(variance: float, bits: int) -> float:
    if variance <= 0:
        return 0.0
    return math.ldexp(float(variance), 2 * bits)


def identity_key_switch_variance(ksP, input_variance: float = 0.0) -> float:
    propagated = scale_variance_between(input_variance, ksP.domainP, ksP.targetP)
    return propagated + float(iksnoisecalc(ksP))


def pbs_variance(brP) -> float:
    return float(brnoisecalc(brP))


def _key_square_coefficient_variance(P) -> float:
    second_moment = P.variance_key_coefficient + P.expectation_key_coefficient**2
    return float(P.n * second_moment * second_moment)


def clpx_relinearization_variance(P=default_params.lvl2param) -> float:
    """Approximate coefficient variance added by CLPX Relinearization().

    CLPXMult relinearizes the c2 term through the standard TFHEpp
    relinKeySwitch path.  This treats that operation as an external product
    with a HalfTRGSW encryption of s^2.  It is a screening model, not a sampled
    validation of FFT/rescale error.
    """
    return float(
        extpnoisecalc(
            P,
            P.σ,
            0.0,
            0.0,
            _key_square_coefficient_variance(P),
        )
    )


def hom_decomp_variances(
    input_variance: float,
    *,
    high2midP=default_params.lvl21param,
    mid2lowP=default_params.lvl1hparam,
    brP=default_params.lvlh1param,
    basebit: int = 2,
    numdigit: int = 6,
) -> HomDecompEstimate:
    domain_bits = q_bits(high2midP.domainP)
    subtlwe_variance = 0.0
    sums: list[float] = []
    max_mid_input = 0.0
    max_pbs_input = 0.0

    for digit in range(1, numdigit + 1):
        shifted_input = shift_variance(input_variance, domain_bits - basebit * digit)
        cres = identity_key_switch_variance(high2midP, shifted_input)
        if digit != 1:
            cres += subtlwe_variance
        sums.append(cres)

        max_mid_input = max(max_mid_input, cres)
        tlwelvlhalf = identity_key_switch_variance(mid2lowP, cres)
        max_pbs_input = max(max_pbs_input, tlwelvlhalf)
        if digit != numdigit:
            subtlwe_variance = pbs_variance(brP)

    return HomDecompEstimate(
        sums=sums,
        max_mid_to_low_input_variance=max_mid_input,
        max_sub_pbs_input_variance=max_pbs_input,
    )


def estimate_tlwes_to_clpx(
    *,
    validbit: int = 8,
    num_multi: int = 4,
    shift: int = 0,
    iksP=default_params.lvl10param,
    bkP=default_params.lvl02param,
    sskP=default_params.lvl22param,
    input_variance: float | None = None,
) -> TLWES2CLPXEstimate:
    if validbit <= 0:
        raise ValueError("validbit must be positive")
    if num_multi <= 0:
        raise ValueError("num_multi must be positive")
    if shift < 0:
        raise ValueError("shift must be non-negative")
    if sskP.domainP is not bkP.targetP:
        raise ValueError("sskP.domainP must match bkP.targetP")

    input_var = float(iksP.domainP.σ if input_variance is None else input_variance)
    iks_var = identity_key_switch_variance(iksP, input_var)
    pbs_var = pbs_variance(bkP)

    target_bits = q_bits(bkP.targetP)
    temp_count = validbit + target_bits - 1
    temp_vars = [0.0 for _ in range(temp_count)]
    term_counts = [0 for _ in range(temp_count)]
    step = num_multi * (shift + 1)

    for input_index in range(validbit):
        for j in range(target_bits - 1, -1, -step):
            for lane in range(num_multi):
                for small_shift in range(shift + 1):
                    bit_index = j - lane * (shift + 1) - small_shift
                    if 0 <= bit_index < target_bits:
                        out_index = input_index + bit_index
                        temp_vars[out_index] += shift_variance(pbs_var, small_shift)
                        term_counts[out_index] += 1

    max_temp = max(temp_vars) if temp_vars else 0.0
    packed_var = float(annihilatecalc(bkP.targetP, max_temp))
    value_var, log2_fail = clpx_digit_failure_log2(
        bkP.targetP, packed_var, count=validbit
    )
    return TLWES2CLPXEstimate(
        input_variance=input_var,
        iks_variance=iks_var,
        pbs_variance=pbs_var,
        max_packed_input_variance=max_temp,
        packed_variance=packed_var,
        clpx_value_variance=value_var,
        temp_tlwe_count=temp_count,
        max_terms_per_temp=max(term_counts) if term_counts else 0,
        log2_failure=log2_fail,
    )


def estimate_clpx_to_tlwes(
    *,
    validbit: int = 8,
    batch_size: int = 16,
    numdigit: int = 4,
    basebit: int = 2,
    iksP10=default_params.lvl1hparam,
    iksP21=default_params.lvl21param,
    bkP01=default_params.lvlh1param,
    bkP02=default_params.lvlh2param,
    iksP20=default_params.lvl2hparam,
    input_variance: float | None = None,
) -> CLPX2TLWESEstimate:
    if validbit <= 0:
        raise ValueError("validbit must be positive")
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if numdigit < 2:
        raise ValueError("numdigit must be at least 2")
    if basebit <= 0:
        raise ValueError("basebit must be positive")

    input_var = float(iksP20.domainP.σ if input_variance is None else input_variance)
    pbs02_var = pbs_variance(bkP02)
    pbs01_var = pbs_variance(bkP01)
    max_internal_input = 0.0
    max_final_input = 0.0
    max_sum_var = 0.0
    output_vars: list[float] = []

    remaining = validbit
    global_index = 0
    epoch_count = (validbit + batch_size - 1) // batch_size
    sumpra_var = 0.0
    prev_temps1_var = 0.0

    for epoch in range(epoch_count):
        current_batch = min(batch_size, remaining)
        for _ in range(current_batch):
            temp = identity_key_switch_variance(iksP20, input_var)
            max_internal_input = max(max_internal_input, temp)

            temp31 = pbs02_var
            temp = identity_key_switch_variance(iksP20, temp31)
            max_internal_input = max(max_internal_input, temp)

            temps0 = pbs02_var
            temps1 = pbs02_var
            diff = temps0 + (prev_temps1_var if global_index > 0 else 0.0)
            temp10 = identity_key_switch_variance(iksP20, diff)

            shifted = shift_variance(temp10, 2)
            max_internal_input = max(max_internal_input, shifted)
            temp32_0 = pbs02_var
            temp32_1 = pbs02_var
            temp10 += identity_key_switch_variance(iksP20, temp32_1)

            max_internal_input = max(
                max_internal_input,
                shift_variance(temp10, 1),
                temp10,
            )
            sumpra_var += pbs02_var + pbs02_var + temp32_0
            prev_temps1_var = temps1
            global_index += 1

        hom = hom_decomp_variances(
            sumpra_var,
            high2midP=iksP21,
            mid2lowP=iksP10,
            brP=bkP01,
            basebit=basebit,
            numdigit=numdigit + 2,
        )
        max_sum_var = max(max_sum_var, max(hom.sums))
        max_internal_input = max(max_internal_input, hom.max_sub_pbs_input_variance)

        for digit in range(1, numdigit):
            temp = identity_key_switch_variance(iksP10, hom.sums[digit])
            for k in range(basebit - 1):
                final_input = shift_variance(temp, basebit - k - 1)
                max_final_input = max(max_final_input, final_input)
                output_vars.append(pbs01_var)
            max_final_input = max(max_final_input, temp)
            output_vars.append(pbs01_var)

        remaining -= current_batch
        if epoch < epoch_count - 1:
            carry_inputs = [
                identity_key_switch_variance(iksP10, hom.sums[numdigit + i])
                for i in range(2)
            ]
            max_internal_input = max(max_internal_input, *carry_inputs)
            max_internal_input = max(
                max_internal_input,
                shift_variance(carry_inputs[0], 1),
            )
            sumpra_var = 3.0 * pbs02_var

    output_var = max(output_vars) if output_vars else 0.0
    log2_fail = tlwe_failure_log2(bkP01.targetP, output_var, count=len(output_vars))
    return CLPX2TLWESEstimate(
        input_variance=input_var,
        pbs02_variance=pbs02_var,
        pbs01_variance=pbs01_var,
        sumpra_variance=sumpra_var,
        max_homdecomp_sum_variance=max_sum_var,
        max_internal_pbs_input_variance=max_internal_input,
        max_final_pbs_input_variance=max_final_input,
        output_variance=output_var,
        produced_tlwes=len(output_vars),
        log2_failure=log2_fail,
    )


def clpx_multiply_value_variance(
    lhs_value_variance: float,
    rhs_value_variance: float,
    *,
    P=default_params.lvl2param,
    message_bound: float | None = None,
    relin_coefficient_variance: float | None = None,
    fft_round_coefficient_variance: float = 1.0 / 12.0,
) -> float:
    """Approximate digit-value variance after one CLPXMult.

    The plaintext digit values of both operands are assumed bounded in absolute
    value by `message_bound`.  TFHEpp's CLPX multiplication rescales by
    `P::plain_modulus`, so the first-order message-noise terms are divided by
    that value squared.  The caller must choose a bound matching the application
    semantics; growing integer products are not modeled automatically.
    """
    if message_bound is None:
        message_bound = float(P.plain_modulus - 1)
    if relin_coefficient_variance is None:
        relin_coefficient_variance = clpx_relinearization_variance(P)

    plain = float(P.plain_modulus)
    # clpxSymIntDecrypt divides the digit-value expression by q/2.  Do the
    # multiplication recurrence in normalized digit units, then scale back to
    # torus-value units for reporting and failure estimates.
    norm_scale = (P.q / 2.0) ** 2
    lhs_normalized = lhs_value_variance / norm_scale
    rhs_normalized = rhs_value_variance / norm_scale
    message_term = (
        message_bound * message_bound * (lhs_normalized + rhs_normalized)
        + lhs_normalized * rhs_normalized
    ) / (plain * plain)
    coefficient_to_value = 1.0 + plain * plain
    relin_term = coefficient_to_value * (
        relin_coefficient_variance + fft_round_coefficient_variance
    ) / norm_scale
    return (message_term + relin_term) * norm_scale


def estimate_clpx_multiplication_depth(
    *,
    initial_coefficient_variance: float | None = None,
    validbit: int = 8,
    max_multiplications: int = 16,
    chain: str = "fresh",
    P=default_params.lvl2param,
    message_bound: float | None = None,
    d: float = 6.0,
) -> CLPXMultiplicationEstimate:
    if max_multiplications < 0:
        raise ValueError("max_multiplications must be non-negative")
    if chain not in {"fresh", "square"}:
        raise ValueError("chain must be 'fresh' or 'square'")
    if d <= 0:
        raise ValueError("d must be positive")

    if initial_coefficient_variance is None:
        initial_coefficient_variance = estimate_tlwes_to_clpx(
            validbit=validbit
        ).packed_variance
    if message_bound is None:
        message_bound = float(P.plain_modulus - 1)

    relin_var = clpx_relinearization_variance(P)
    fft_round_var = 1.0 / 12.0
    coeff_to_value = 1.0 + float(P.plain_modulus) ** 2
    initial_value_var = coeff_to_value * initial_coefficient_variance
    current_value_var = initial_value_var
    fresh_value_var = initial_value_var
    steps: list[CLPXMultiplicationStep] = []
    supported = 0

    for count in range(0, max_multiplications + 1):
        coeff_var = current_value_var / coeff_to_value
        margin = clpx_value_margin_log2(P, current_value_var, d)
        log2_fail = _failure_log2(P.q / 4.0, current_value_var, count=1)
        status = "OK" if margin >= 0 else "FAIL"
        if status == "OK":
            supported = count
        steps.append(
            CLPXMultiplicationStep(
                multiplication_count=count,
                coefficient_variance=coeff_var,
                digit_value_variance=current_value_var,
                margin_bits=margin,
                log2_failure=log2_fail,
                status=status,
            )
        )
        if count == max_multiplications:
            break
        rhs = fresh_value_var if chain == "fresh" else current_value_var
        current_value_var = clpx_multiply_value_variance(
            current_value_var,
            rhs,
            P=P,
            message_bound=message_bound,
            relin_coefficient_variance=relin_var,
            fft_round_coefficient_variance=fft_round_var,
        )

    return CLPXMultiplicationEstimate(
        initial_coefficient_variance=initial_coefficient_variance,
        initial_digit_value_variance=initial_value_var,
        relin_coefficient_variance=relin_var,
        fft_round_coefficient_variance=fft_round_var,
        message_bound=message_bound,
        chain=chain,
        d=d,
        steps=steps,
        supported_multiplications=supported,
    )
