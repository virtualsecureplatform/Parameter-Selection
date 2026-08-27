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

from scipy.special import erfcx

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
    approximation_variance: float
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
    rounded_digit_margin_bits: float
    homdecomp_margin_bits: float
    final_extraction_margin_bits: float
    carry_margin_bits: float
    output_variance: float
    produced_tlwes: int
    fid_round_pbs_count: int
    homdecomp_pbs_count: int
    bit_extraction_pbs_count: int
    carry_pbs_count: int
    total_pbs_count: int
    semantic_log2_failure: float
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


def _failure_log2(threshold: float, variance: float, count: int = 1) -> float:
    if variance <= 0:
        return NEG_INF
    if threshold <= 0:
        return 0.0
    # erfc underflows for precisely the small failure probabilities for which
    # this estimator is useful.  erfc(x) = exp(-x^2) * erfcx(x) keeps the
    # calculation in the log domain.
    x = threshold / math.sqrt(2.0 * variance)
    log_probability = math.log(float(erfcx(x))) - x * x + math.log(count)
    return min(0.0, log_probability / math.log(2.0))


def _union_log2(log_probabilities: list[float]) -> float:
    finite = [value for value in log_probabilities if value != NEG_INF]
    if not finite:
        return NEG_INF
    largest = max(finite)
    if largest >= 0:
        return 0.0
    total = sum(2.0 ** (value - largest) for value in finite)
    return min(0.0, largest + math.log2(total))


def _margin_failure_log2(margin_bits: float, count: int) -> float:
    if count <= 0 or math.isinf(margin_bits):
        return NEG_INF
    return _failure_log2(2.0 ** (margin_bits / 2.0), 1.0, count)


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


def pbs_input_margin_log2(
    brP,
    input_variance: float,
    num_out: int = 1,
    input_plain_modulus: int | None = None,
) -> float:
    """Variance-bit margin to the nearest semantic PBS decision boundary.

    ``num_out`` is retained for API compatibility, but ManyLUT's coefficient
    layout does not shrink the input message interval.  The previous model
    compared TLWE noise with one modulus-switch rounding bin; ordinary TFHE
    PBS noise is not required to fit inside such a bin.
    """
    if input_variance <= 0:
        return float("inf")
    del num_out
    plain_modulus = int(
        brP.domainP.plain_modulus
        if input_plain_modulus is None
        else input_plain_modulus
    )
    half_interval = float(brP.domainP.q) / (2.0 * plain_modulus)
    return 2.0 * math.log2(half_interval) - math.log2(input_variance)


def decision_margin_log2(threshold: float, variance: float) -> float:
    """Return log2(threshold^2 / variance) for a custom PBS decision."""
    if variance <= 0:
        return float("inf")
    if threshold <= 0:
        return NEG_INF
    return 2.0 * math.log2(threshold) - math.log2(variance)


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


def reverse_lvl2_pbs_variance(
    bgbit: int,
    *,
    levels: int = 4,
    domainP=default_params.lvlhalfparam,
    targetP=default_params.lvl2param,
) -> float:
    """Evaluate a reverse-switch lvl2 PBS gadget decomposition candidate."""
    if bgbit <= 0 or levels <= 0:
        raise ValueError("bgbit and levels must be positive")
    base = 1 << bgbit
    candidate_target = type(
        f"ReverseLvl2_l{levels}_b{bgbit}",
        (targetP,),
        {
            "l": levels,
            "la": levels,
            "Bbit": bgbit,
            "Babit": bgbit,
            "B": base,
            "Ba": base,
        },
    )
    candidate_bk = type(
        f"ReverseLvlh2_l{levels}_b{bgbit}",
        (),
        {"domainP": domainP, "targetP": candidate_target},
    )
    return pbs_variance(candidate_bk)


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
    source_bits: int | None = None,
) -> HomDecompEstimate:
    domain_bits = q_bits(high2midP.domainP)
    if source_bits is None:
        source_bits = domain_bits
    if source_bits < basebit * numdigit or source_bits > domain_bits:
        raise ValueError("invalid HomDecomp source window")
    subtlwe_variance = 0.0
    sums: list[float] = []
    max_mid_input = 0.0
    max_pbs_input = 0.0

    for digit in range(1, numdigit + 1):
        shifted_input = shift_variance(input_variance, source_bits - basebit * digit)
        cres = identity_key_switch_variance(high2midP, shifted_input)
        if digit != 1:
            cres += subtlwe_variance
        sums.append(cres)

        max_mid_input = max(max_mid_input, cres)
        tlwelvlhalf = identity_key_switch_variance(mid2lowP, cres)
        if digit != numdigit:
            max_pbs_input = max(max_pbs_input, tlwelvlhalf)
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
    w: int | None = None,
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
    truncation_bits = target_bits if w is None else int(w)
    if truncation_bits <= 0:
        raise ValueError("w must be positive")
    if truncation_bits > target_bits:
        raise ValueError("w cannot exceed the target torus bit width")

    temp_count = validbit + truncation_bits - 1
    temp_vars = [0.0 for _ in range(temp_count)]
    term_counts = [0 for _ in range(temp_count)]
    step = num_multi * (shift + 1)

    for input_index in range(validbit):
        for j in range(truncation_bits - 1, -1, -step):
            for lane in range(num_multi):
                for small_shift in range(shift + 1):
                    bit_index = j - lane * (shift + 1) - small_shift
                    if 0 <= bit_index < truncation_bits:
                        out_index = input_index + bit_index
                        temp_vars[out_index] += shift_variance(pbs_var, small_shift)
                        term_counts[out_index] += 1

    max_temp = max(temp_vars) if temp_vars else 0.0
    # Nagai et al., T2CPVar (Eq. 43 in the checked-in draft): approximating
    # Delta_b to w bits contributes q^2 / 2^(2w) before annihilate packing.
    # This term is what creates the observed correctness transition around
    # w=16--20; omitting it incorrectly makes very small w look preferable.
    approximation_var = math.ldexp(
        1.0, 2 * (target_bits - truncation_bits)
    )
    packed_var = float(
        annihilatecalc(bkP.targetP, max_temp + approximation_var)
    )
    value_var, log2_fail = clpx_digit_failure_log2(
        bkP.targetP, packed_var, count=validbit
    )
    return TLWES2CLPXEstimate(
        input_variance=input_var,
        iks_variance=iks_var,
        pbs_variance=pbs_var,
        approximation_variance=approximation_var,
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
    batch_size: int | None = None,
    numdigit: int = 9,
    basebit: int = 2,
    carry_mode: str = "legacy",
    iksP10=default_params.lvl1hparam,
    iksP21=default_params.lvl21param,
    bkP01=default_params.lvlh1param,
    bkP02=default_params.lvlh2param,
    iksP20=default_params.lvl2hparam,
    input_variance: float | None = None,
) -> CLPX2TLWESEstimate:
    if validbit <= 0:
        raise ValueError("validbit must be positive")
    if numdigit < 2:
        raise ValueError("numdigit must be at least 2")
    if basebit <= 0:
        raise ValueError("basebit must be positive")
    if carry_mode not in {"legacy", "single"}:
        raise ValueError("carry_mode must be 'legacy' or 'single'")
    if carry_mode == "legacy" and basebit != 2:
        raise ValueError("the legacy two-digit carry is modeled only for basebit=2")
    if carry_mode == "single" and basebit < 4:
        raise ValueError("the single-PBS carry requires basebit >= 4")
    implementation_block_size = (numdigit - 1) * basebit
    if batch_size is None:
        batch_size = implementation_block_size
    elif batch_size != implementation_block_size:
        raise ValueError(
            "batch_size must equal (numdigit - 1) * basebit to match TFHEpp"
        )

    input_var = float(iksP20.domainP.σ if input_variance is None else input_variance)
    pbs02_var = pbs_variance(bkP02)
    pbs01_var = pbs_variance(bkP01)
    max_internal_input = 0.0
    max_final_input = 0.0
    max_sum_var = 0.0
    rounded_margin = float("inf")
    hom_margin = float("inf")
    final_margin = float("inf")
    carry_margin = float("inf")
    output_vars: list[float] = []
    fid_round_pbs_count = 0
    homdecomp_pbs_count = 0
    bit_extraction_pbs_count = 0
    carry_pbs_count = 0

    remaining = validbit
    global_index = 0
    epoch_count = (validbit + batch_size - 1) // batch_size
    sumpra_var = 0.0
    previous_scaled_var = 0.0

    for epoch in range(epoch_count):
        current_batch = min(batch_size, remaining)
        for _ in range(current_batch):
            # Fid_1 refreshes the extracted coefficient. Fid_2 refreshes its
            # scaled representative before x-b is formed. The current Fid_2
            # ciphertext is shifted once, whereas the previous coefficient is
            # used without a shift.
            temp = identity_key_switch_variance(iksP20, input_var)
            max_internal_input = max(max_internal_input, temp)

            fid_stage1 = pbs02_var
            temp = identity_key_switch_variance(iksP20, fid_stage1)
            max_internal_input = max(max_internal_input, temp)

            fid_stage2 = pbs02_var
            x_minus_b = shift_variance(fid_stage2, 1)
            if global_index > 0:
                x_minus_b += previous_scaled_var
            rounded_input = identity_key_switch_variance(iksP20, x_minus_b)
            max_internal_input = max(max_internal_input, rounded_input)
            rounded_threshold = float(iksP20.targetP.q) / 32.0
            rounded_margin = min(
                rounded_margin,
                decision_margin_log2(rounded_threshold, rounded_input),
            )

            # The combined rounded-digit/weight LUT replaces the previous
            # three-PBS chain with one refreshed lvl2 ciphertext.
            sumpra_var += pbs02_var
            previous_scaled_var = fid_stage2
            global_index += 1
            fid_round_pbs_count += 3

        carry_digits = 2 if carry_mode == "legacy" else 1
        hom_numdigit = numdigit + carry_digits
        hom = hom_decomp_variances(
            sumpra_var,
            high2midP=iksP21,
            mid2lowP=iksP10,
            brP=bkP01,
            basebit=basebit,
            numdigit=hom_numdigit,
            source_bits=basebit * hom_numdigit,
        )
        homdecomp_pbs_count += hom_numdigit - 1
        max_sum_var = max(max_sum_var, max(hom.sums))
        max_internal_input = max(max_internal_input, hom.max_sub_pbs_input_variance)
        hom_threshold = float(bkP01.domainP.q) / (1 << (basebit + 1))
        hom_margin = min(
            hom_margin,
            decision_margin_log2(
                hom_threshold, hom.max_sub_pbs_input_variance
            ),
        )

        for digit in range(1, numdigit):
            temp = identity_key_switch_variance(iksP10, hom.sums[digit])
            active_bits = min(
                basebit,
                max(0, current_batch - (digit - 1) * basebit),
            )
            if active_bits == 0:
                continue

            for k in range(active_bits):
                shift = basebit - k - 1
                final_input = shift_variance(temp, shift)
                threshold = float(iksP10.targetP.q) / (
                    1 << (basebit + 1 - shift)
                )
                max_final_input = max(max_final_input, final_input)
                final_margin = min(
                    final_margin,
                    decision_margin_log2(threshold, final_input),
                )
                output_vars.append(pbs01_var)
            bit_extraction_pbs_count += active_bits

        remaining -= current_batch
        if epoch < epoch_count - 1:
            if carry_mode == "legacy":
                carry_inputs = [
                    identity_key_switch_variance(iksP10, hom.sums[numdigit + i])
                    for i in range(2)
                ]
                shifted_carry = shift_variance(carry_inputs[0], 1)
                max_internal_input = max(
                    max_internal_input, *carry_inputs, shifted_carry
                )
                unshifted_threshold = float(iksP10.targetP.q) / 8.0
                shifted_threshold = float(iksP10.targetP.q) / 4.0
                carry_margin = min(
                    carry_margin,
                    *(decision_margin_log2(unshifted_threshold, variance)
                      for variance in carry_inputs),
                    decision_margin_log2(shifted_threshold, shifted_carry),
                )
                sumpra_var = 3.0 * pbs02_var
                carry_pbs_count += 3
            else:
                carry_input = identity_key_switch_variance(
                    iksP10, hom.sums[numdigit]
                )
                max_internal_input = max(max_internal_input, carry_input)
                carry_threshold = float(iksP10.targetP.q) / (
                    1 << (basebit + 1)
                )
                carry_margin = min(
                    carry_margin,
                    decision_margin_log2(carry_threshold, carry_input),
                )
                sumpra_var = pbs02_var
                carry_pbs_count += 1

    output_var = max(output_vars) if output_vars else 0.0
    log2_fail = tlwe_failure_log2(bkP01.targetP, output_var, count=len(output_vars))
    total_pbs_count = (
        fid_round_pbs_count
        + homdecomp_pbs_count
        + bit_extraction_pbs_count
        + carry_pbs_count
    )
    semantic_log2_failure = _union_log2(
        [
            _margin_failure_log2(rounded_margin, fid_round_pbs_count // 3),
            _margin_failure_log2(hom_margin, homdecomp_pbs_count),
            _margin_failure_log2(final_margin, len(output_vars)),
            _margin_failure_log2(carry_margin, carry_pbs_count),
        ]
    )
    return CLPX2TLWESEstimate(
        input_variance=input_var,
        pbs02_variance=pbs02_var,
        pbs01_variance=pbs01_var,
        sumpra_variance=sumpra_var,
        max_homdecomp_sum_variance=max_sum_var,
        max_internal_pbs_input_variance=max_internal_input,
        max_final_pbs_input_variance=max_final_input,
        rounded_digit_margin_bits=rounded_margin,
        homdecomp_margin_bits=hom_margin,
        final_extraction_margin_bits=final_margin,
        carry_margin_bits=carry_margin,
        output_variance=output_var,
        produced_tlwes=len(output_vars),
        fid_round_pbs_count=fid_round_pbs_count,
        homdecomp_pbs_count=homdecomp_pbs_count,
        bit_extraction_pbs_count=bit_extraction_pbs_count,
        carry_pbs_count=carry_pbs_count,
        total_pbs_count=total_pbs_count,
        semantic_log2_failure=semantic_log2_failure,
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


def paper_clpx_multiplication_variance(
    lhs_coefficient_variance: float,
    rhs_coefficient_variance: float,
    lhs_nonzero_terms: int,
    rhs_nonzero_terms: int,
    *,
    P=default_params.SS2CLPXlvl2param,
) -> float:
    """Equation (44) from Nagai et al. for CLPX multiplication noise."""
    b = float(P.plain_modulus)
    common = 6.0 * P.n + math.sqrt(3.0 * P.n)
    lhs_factor = ((b + 1.0) * lhs_nonzero_terms + common) ** 2
    rhs_factor = ((b + 1.0) * rhs_nonzero_terms + common) ** 2
    return (
        lhs_factor * rhs_coefficient_variance
        + rhs_factor * lhs_coefficient_variance
    )


def paper_clpx_fixed_noise_threshold(
    *,
    P=default_params.SS2CLPXlvl2param,
    relinearization_base: int | None = None,
) -> float:
    """Equations (45)-(47) from Nagai et al.

    The paper's ``w_rlin`` is the integer base used to decompose the
    relinearized polynomial, not the unrelated Delta_b truncation width ``w``.
    TFHEpp has ``P.l`` gadget rows whose signed digits use base ``P.B``.
    """
    b = float(P.plain_modulus)
    w_rlin = int(P.B if relinearization_base is None else relinearization_base)
    if w_rlin <= 0:
        raise ValueError("relinearization_base must be positive")
    a1 = (b + 1.0) * math.sqrt(3.0 * P.n) * (
        1.0 + 2.0 * math.sqrt(3.0 * P.n) + 12.0 * P.n
    )
    a2 = (
        6.0
        * math.sqrt(3.0)
        * (b + 1.0)
        * P.n
        * P.alpha
        * P.l
        * w_rlin
    )
    return P.q / 2.0 - a1 - a2


def estimate_clpx_multiplication_depth(
    *,
    initial_coefficient_variance: float | None = None,
    validbit: int = 8,
    num_multi: int = 4,
    shift: int = 0,
    w: int | None = None,
    max_multiplications: int = 16,
    chain: str = "fresh",
    P=default_params.lvl2param,
    iksP=default_params.lvl10param,
    bkP=default_params.lvl02param,
    sskP=default_params.lvl22param,
    message_bound: float | None = None,
    d: float = 6.0,
    model: str = "tfhepp",
) -> CLPXMultiplicationEstimate:
    if max_multiplications < 0:
        raise ValueError("max_multiplications must be non-negative")
    if chain not in {"fresh", "square"}:
        raise ValueError("chain must be 'fresh' or 'square'")
    if d <= 0:
        raise ValueError("d must be positive")
    if model not in {"tfhepp", "paper"}:
        raise ValueError("model must be 'tfhepp' or 'paper'")

    if initial_coefficient_variance is None:
        initial_coefficient_variance = estimate_tlwes_to_clpx(
            validbit=validbit,
            num_multi=num_multi,
            shift=shift,
            w=w,
            iksP=iksP,
            bkP=bkP,
            sskP=sskP,
        ).packed_variance
    if message_bound is None:
        message_bound = float(P.plain_modulus - 1)

    relin_var = clpx_relinearization_variance(P)
    fft_round_var = 1.0 / 12.0
    coeff_to_value = 1.0 + float(P.plain_modulus) ** 2
    initial_value_var = coeff_to_value * initial_coefficient_variance

    if model == "paper":
        threshold = paper_clpx_fixed_noise_threshold(P=P)
        current_coeff_var = initial_coefficient_variance
        fresh_coeff_var = initial_coefficient_variance
        current_nonzero_terms = validbit
        fresh_nonzero_terms = validbit
        steps: list[CLPXMultiplicationStep] = []
        supported = 0
        for count in range(0, max_multiplications + 1):
            output_variance = (
                (float(P.plain_modulus) + 1.0) ** 2 * current_coeff_var
            )
            margin = (
                2.0 * math.log2(threshold) - math.log2(output_variance)
                if threshold > 0 and output_variance > 0
                else float("-inf")
            )
            log2_fail = (
                _failure_log2(threshold, output_variance, count=1)
                if threshold > 0
                else 0.0
            )
            status = "OK" if log2_fail <= -64.0 else "FAIL"
            if status == "OK":
                supported = count
            steps.append(
                CLPXMultiplicationStep(
                    multiplication_count=count,
                    coefficient_variance=current_coeff_var,
                    digit_value_variance=output_variance,
                    margin_bits=margin,
                    log2_failure=log2_fail,
                    status=status,
                )
            )
            if count == max_multiplications:
                break
            rhs_coeff_var = (
                fresh_coeff_var if chain == "fresh" else current_coeff_var
            )
            rhs_nonzero_terms = (
                fresh_nonzero_terms if chain == "fresh" else current_nonzero_terms
            )
            current_coeff_var = paper_clpx_multiplication_variance(
                current_coeff_var,
                rhs_coeff_var,
                current_nonzero_terms,
                rhs_nonzero_terms,
                P=P,
            )
            current_nonzero_terms += rhs_nonzero_terms

        return CLPXMultiplicationEstimate(
            initial_coefficient_variance=initial_coefficient_variance,
            initial_digit_value_variance=(
                (float(P.plain_modulus) + 1.0) ** 2
                * initial_coefficient_variance
            ),
            relin_coefficient_variance=relin_var,
            fft_round_coefficient_variance=fft_round_var,
            message_bound=message_bound,
            chain=chain,
            d=d,
            steps=steps,
            supported_multiplications=supported,
        )

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
