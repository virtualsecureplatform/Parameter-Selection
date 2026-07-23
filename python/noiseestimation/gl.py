#!/usr/bin/env python3
"""Noise screening model for TFHEpp's Gentry--Lee SHIP bootstrap.

The GL and GL-SHIP papers prove algebraic correctness but deliberately collect
the implementation error in an abstract ``epsilon_HE`` term.  This module
expands that term for the coefficient-domain Double Decomposition (DD) path in
TFHEpp.  It is an average-case variance model, not a proof of correctness or
security.

The model tracks two domains:

* coefficient phase variance through grouped slots-to-coefficients (StC) and
  dense-to-sparse switching; and
* decoded (X,W)-slot variance through masked-column selection, X-only HMux,
  the balanced product tree, conjugation, and Y reassembly.

All variances are stored in log2-space.  This is required for the 214--861 bit
``P*Q`` profiles and also makes parameter sweeps numerically stable.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Iterable, Mapping


NEG_INF = float("-inf")


def log2_add(*values: float) -> float:
    """Return ``log2(sum(2**value))`` without overflowing."""
    finite = [value for value in values if value != NEG_INF]
    if not finite:
        return NEG_INF
    maximum = max(finite)
    return maximum + math.log2(
        sum(2.0 ** (value - maximum) for value in finite)
    )


def log2_sum(values: Iterable[float]) -> float:
    return log2_add(*tuple(values))


def log2_sum_stddevs(values: Iterable[float]) -> float:
    """Return the variance of a worst-aligned sum from component variances."""
    finite = [value for value in values if value != NEG_INF]
    if not finite:
        return NEG_INF
    return 2.0 * log2_sum(0.5 * value for value in finite)


def combine_variances(values: Iterable[float], model: str) -> float:
    values = tuple(values)
    if model == "independent":
        return log2_sum(values)
    if model == "correlated":
        return log2_sum_stddevs(values)
    raise ValueError(f"unknown variance model: {model}")


def repeated_variance(log2_variance: float, count: float, model: str) -> float:
    if count <= 0.0 or log2_variance == NEG_INF:
        return NEG_INF
    multiplier = math.log2(count)
    if model == "correlated":
        multiplier *= 2.0
    return log2_variance + multiplier


def _sum_gain(count: int, model: str) -> float:
    """Variance gain for a sum of ``count`` equal-scale values."""
    if count <= 0:
        raise ValueError("sum length must be positive")
    return float(count if model == "independent" else count * count)


def _ceil_div(value: int, divisor: int) -> int:
    return (value + divisor - 1) // divisor


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


@dataclass(frozen=True)
class GLNoiseParams:
    """Parameters needed by the TFHEpp GL/DD noise model.

    ``q0_bits`` and the per-level scale sizes are reconstructed defaults.  The
    GL-SHIP paper publishes total ``log Q``, ``log P``, the grouped StC limb,
    and the gap, but not every individual RNS prime.  The CLI exposes all of
    these values so a concrete implementation schedule can override them.
    """

    tag: str
    n: int
    p: int
    log_q: int
    log_p: int
    q0_bits: int
    gap_bits: int
    stc_bits: int
    x_scale_bits: int
    w_scale_bits: int
    tree_scale_bits: int
    outside_multiplicative_depth: int
    outside_scale_bits: int
    dnum: int
    theta: int
    window_width: int
    masked_column_count: int
    primary_bit: int
    bbar_bit: int
    storage_bits: int
    security_limit_log_pq: int
    paper_precision_bits: float
    error_stddev: float = 3.2
    sparse_hamming_weight: int = 31
    hmux_radix: int = 4
    w_baby_step: int = 4
    dense_secret_variance: float = 2.0 / 3.0
    message_bound: float = 1.0

    @property
    def phi(self) -> int:
        return self.p - 1

    @property
    def ring_degree(self) -> int:
        return 2 * self.n * self.phi

    @property
    def full_ring_degree(self) -> int:
        return self.n * self.ring_degree

    @property
    def factor_count(self) -> int:
        return self.sparse_hamming_weight + 1

    @property
    def tree_depth(self) -> int:
        return (self.factor_count - 1).bit_length()

    @property
    def input_log_q(self) -> int:
        return self.q0_bits + self.stc_bits

    @property
    def input_scale_bits(self) -> int:
        return self.q0_bits - self.gap_bits

    @property
    def output_log_q(self) -> int:
        return self.log_q - self.tree_depth * self.tree_scale_bits

    @property
    def required_output_log_q(self) -> int:
        return (
            self.tree_scale_bits
            + self.outside_multiplicative_depth * self.outside_scale_bits
        )

    @property
    def output_depth_margin_bits(self) -> int:
        return self.output_log_q - self.required_output_log_q

    @property
    def tree_scale_headroom_bits(self) -> int:
        # Raising every tree scale by one bit consumes tree_depth bits from
        # output_log_q and raises the required output scale by one more bit.
        return max(0, self.output_depth_margin_bits // (self.tree_depth + 1))

    @property
    def gamma(self) -> int:
        return 1 << self.gap_bits

    @property
    def w_giant_steps(self) -> int:
        return _ceil_div(self.phi, self.w_baby_step)

    @property
    def hmux_stages(self) -> int:
        choices = self.n // self.theta
        stages = 0
        covered = 1
        while covered < choices:
            covered *= self.hmux_radix
            stages += 1
        return stages

    @property
    def average_candidates_per_sparse_term(self) -> float:
        return self.masked_column_count / self.sparse_hamming_weight

    @property
    def log_pq(self) -> int:
        return self.log_q + self.log_p

    @property
    def storage_margin_bits(self) -> int:
        return self.storage_bits - self.log_pq

    @property
    def security_margin_bits(self) -> int:
        return self.security_limit_log_pq - self.log_pq

    @property
    def bbar_cover_bits(self) -> int:
        return _ceil_div(self.log_pq, self.bbar_bit) * self.bbar_bit

    @property
    def primary_cover_bits(self) -> int:
        return _ceil_div(self.log_q, self.primary_bit) * self.primary_bit

    def validate(self) -> None:
        if self.n <= 0 or self.n & (self.n - 1):
            raise ValueError("n must be a positive power of two")
        if self.p < 3 or not _is_prime(self.p):
            raise ValueError("p must be an odd prime for the modeled GL ring")
        if self.ring_degree <= 0:
            raise ValueError("invalid GL ring degree")
        if self.log_q <= 0 or self.log_p <= 0 or self.q0_bits <= 0:
            raise ValueError("modulus sizes must be positive")
        if self.q0_bits >= 63:
            raise ValueError("TFHEpp's phase-root table requires q0 below 2^63")
        if self.gap_bits < 0:
            raise ValueError("encoding gap bits cannot be negative")
        if self.q0_bits <= self.gap_bits:
            raise ValueError("q0 must be wider than the encoding gap")
        if self.x_scale_bits <= 0 or self.w_scale_bits <= 0:
            raise ValueError("StC transform scales must be positive")
        if self.x_scale_bits + self.w_scale_bits != self.stc_bits:
            raise ValueError("X and W transform scales must sum to stc_bits")
        if self.tree_scale_bits <= 0:
            raise ValueError("tree scale must be positive")
        if self.sparse_hamming_weight <= 0:
            raise ValueError("sparse Hamming weight must be positive")
        if self.factor_count & (self.factor_count - 1):
            raise ValueError("sparse_hamming_weight + 1 must be a power of two")
        if self.theta <= 0 or self.n % self.theta:
            raise ValueError("theta must be a positive divisor of n")
        if self.hmux_radix < 2:
            raise ValueError("HMux radix must be at least two")
        if not 0 < self.w_baby_step <= self.phi:
            raise ValueError("invalid W BSGS baby step")
        if self.window_width <= 0 or self.dnum <= 0:
            raise ValueError("window width and paper dnum must be positive")
        if self.outside_multiplicative_depth < 0 or self.outside_scale_bits <= 0:
            raise ValueError("invalid outside-depth schedule")
        if self.primary_bit <= 0 or self.bbar_bit <= 0:
            raise ValueError("DD base bits must be positive")
        if 2 * self.bbar_bit + 2 >= 63:
            raise ValueError("DD auxiliary digit products must fit int64_t")
        if not math.isfinite(self.error_stddev) or self.error_stddev <= 0.0:
            raise ValueError("error standard deviation must be positive")
        if (
            not math.isfinite(self.dense_secret_variance)
            or self.dense_secret_variance < 0.0
        ):
            raise ValueError("secret variance cannot be negative")
        if not math.isfinite(self.message_bound) or self.message_bound <= 0.0:
            raise ValueError("message bound must be positive")
        if self.message_bound >= self.gamma / 2.0:
            raise ValueError("message bound leaves no centered q0 phase headroom")
        if self.masked_column_count <= 0:
            raise ValueError("masked-column count must be positive")
        if self.output_log_q <= self.tree_scale_bits:
            raise ValueError("product tree leaves an unusable output modulus")
        if self.output_log_q < self.required_output_log_q:
            raise ValueError("product tree leaves insufficient outside depth")
        if self.input_log_q + self.log_p > self.storage_bits:
            raise ValueError("StC DD key modulus exceeds torus storage")
        if self.log_q + self.log_p > self.storage_bits:
            raise ValueError("half-bootstrap DD key modulus exceeds torus storage")
        if self.primary_cover_bits > self.storage_bits:
            raise ValueError("DD primary digits do not fit torus storage")
        if self.bbar_cover_bits > self.storage_bits:
            raise ValueError("DD auxiliary limbs do not fit torus storage")

    def with_overrides(self, **overrides: object) -> "GLNoiseParams":
        result = replace(self, **overrides)
        result.validate()
        return result


@dataclass(frozen=True)
class DDKeySwitchNoise:
    log2_eval_key_variance: float
    log2_moddown_rounding_variance: float
    log2_total_variance: float
    primary_rows: int


@dataclass(frozen=True)
class ProductState:
    log2_variance: float
    value_bound: float


@dataclass(frozen=True)
class GLNoiseEstimate:
    params: GLNoiseParams
    model: str
    arithmetic_mode: str
    masked_moddown: str
    tail_bound: float
    stages: Mapping[str, float]
    fresh_phase_log2_variance: float
    stc_phase_log2_variance: float
    half_sparse_phase_log2_variance: float
    full_sparse_phase_log2_variance: float
    half_phase_wrap_margin_bits: float
    full_phase_wrap_margin_bits: float
    half_he_log2_variance: float
    full_he_log2_variance: float
    half_he_error_bound: float
    full_he_error_bound: float
    half_sine_error_bound: float
    full_sine_error_bound: float
    half_quantization_error_bound: float
    full_quantization_error_bound: float
    half_total_error_bound: float
    full_total_error_bound: float
    half_precision_bits: float
    full_precision_bits: float

    @property
    def precision_delta_from_paper(self) -> float:
        return self.full_precision_bits - self.params.paper_precision_bits

    @property
    def sparse_phase_log2_variance(self) -> float:
        """Full-bootstrap sparse phase, kept as the conservative headline."""
        return self.full_sparse_phase_log2_variance

    @property
    def phase_wrap_margin_bits(self) -> float:
        return min(
            self.half_phase_wrap_margin_bits,
            self.full_phase_wrap_margin_bits,
        )


def base_ring_convolution_gain(params: GLNoiseParams) -> int:
    """Maximum coefficient variance gain for one (I,X,W) product.

    The Gaussian/X part contributes ``2*n`` pairs.  Reduction by
    ``Phi_p(W)=1+...+W^(p-1)`` contributes at most ``2*phi(p)-1`` pairs to an
    output W coefficient.
    """
    return 2 * params.n * (2 * params.phi - 1)


def full_ring_convolution_gain(params: GLNoiseParams) -> int:
    """Conservative gain when the Y dimension also participates."""
    return params.n * base_ring_convolution_gain(params)


def dense_secret_phase_gain(params: GLNoiseParams) -> float:
    return base_ring_convolution_gain(params) * params.dense_secret_variance


def sparse_secret_phase_gain(params: GLNoiseParams) -> float:
    # Multiplication by one W monomial has average squared row norm at most
    # (2*phi-1)/phi because one coefficient can fan out through Phi_p.
    monomial_gain = (2.0 * params.phi - 1.0) / params.phi
    return params.sparse_hamming_weight * monomial_gain


def encoded_real_coefficient_variance(params: GLNoiseParams) -> float:
    """Average real-coefficient variance for uniform GL input slots.

    The paper samples each real and imaginary slot component uniformly from
    ``[-message_bound, message_bound]``.  The X and Y embedding matrices are
    unitary up to their factors of ``n``.  For prime ``p``, the degree
    ``p-2`` inverse-W representative has squared row norm ``2/p``.  Hence a
    real GL coefficient has variance ``(B^2/3) * 2/(p*n^2)``.
    """
    return (
        params.message_bound**2
        / 3.0
        * 2.0
        / (params.p * params.n * params.n)
    )


def _phase_rounding_log2_variance(secret_phase_gain: float) -> float:
    # Each rounded ciphertext component has coefficient variance 1/12.  Phase
    # reconstruction is r0 + r1*s.
    return math.log2(1.0 + secret_phase_gain) - math.log2(12.0)


def _tensor_rounding_log2_variance(params: GLNoiseParams) -> float:
    # Component-wise rescaling of (c0*d0, c0*d1, c1*d0, c1*d1) gives
    # r00 + (r01+r10)s + r11*s^2.  The second-moment proxy is
    # (1 + C*E[s^2])^2 / 12.
    phase_gain = dense_secret_phase_gain(params)
    return 2.0 * math.log2(1.0 + phase_gain) - math.log2(12.0)


def dd_key_switch_noise(
    params: GLNoiseParams,
    *,
    log_q: int,
    destination: str,
    full_ring: bool,
) -> DDKeySwitchNoise:
    """Added coefficient phase variance of one TFHEpp DD switch.

    TFHEpp's active primary decomposition and auxiliary Bbar decomposition both
    cover all active bits, so there is no gadget truncation term.  Evaluation
    key error is divided by ``2^log_p``.  The remaining floor is coefficient
    rounding when the two ciphertext components are divided by that modulus.
    """
    if log_q <= 0:
        raise ValueError("key-switch input modulus must be positive")
    if log_q + params.log_p > params.storage_bits:
        raise ValueError("key-switch modulus exceeds torus storage")

    primary_rows = _ceil_div(log_q, params.primary_bit)
    if primary_rows * params.primary_bit > params.storage_bits:
        raise ValueError("key-switch primary digits exceed torus storage")
    digit_log2_variance = log2_add(
        2.0 * params.primary_bit,
        1.0,
    ) - math.log2(12.0)
    convolution_gain = (
        full_ring_convolution_gain(params)
        if full_ring
        else base_ring_convolution_gain(params)
    )
    eval_key_log2 = (
        math.log2(primary_rows)
        + math.log2(convolution_gain)
        + digit_log2_variance
        + 2.0 * math.log2(params.error_stddev)
        - 2.0 * params.log_p
    )

    if destination == "dense":
        secret_gain = dense_secret_phase_gain(params)
    elif destination == "sparse":
        secret_gain = sparse_secret_phase_gain(params)
    else:
        raise ValueError(f"unknown destination secret: {destination}")
    rounding_log2 = _phase_rounding_log2_variance(secret_gain)
    return DDKeySwitchNoise(
        log2_eval_key_variance=eval_key_log2,
        log2_moddown_rounding_variance=rounding_log2,
        log2_total_variance=log2_add(eval_key_log2, rounding_log2),
        primary_rows=primary_rows,
    )


def _embedding_gain(params: GLNoiseParams, model: str) -> float:
    return _sum_gain(params.ring_degree, model)


def _coefficient_to_slot_variance(
    params: GLNoiseParams,
    log2_coefficient_variance: float,
    scale_bits: float,
    model: str,
) -> float:
    if log2_coefficient_variance == NEG_INF:
        return NEG_INF
    return (
        log2_coefficient_variance
        + math.log2(_embedding_gain(params, model))
        - 2.0 * scale_bits
    )


def _stc_phase_noise(
    params: GLNoiseParams, model: str
) -> tuple[float, dict[str, float]]:
    input_log_q = params.input_log_q
    fresh = 2.0 * math.log2(params.error_stddev)
    big_switch = dd_key_switch_noise(
        params, log_q=input_log_q, destination="dense", full_ring=True
    ).log2_total_variance
    small_switch = dd_key_switch_noise(
        params, log_q=input_log_q, destination="dense", full_ring=False
    ).log2_total_variance

    x_gain = _sum_gain(params.n, model)
    w_gain = _sum_gain(params.phi, model)
    input_and_first_switch = log2_add(fresh, big_switch)
    transformed_input = (
        input_and_first_switch + math.log2(x_gain) + math.log2(w_gain)
    )

    second_conjugate = (
        big_switch
        - 2.0 * params.x_scale_bits
        + math.log2(w_gain)
    )
    baby_reuses = (params.w_baby_step - 1) * params.w_giant_steps
    baby_rotations = repeated_variance(
        small_switch - 2.0 * params.x_scale_bits,
        baby_reuses,
        model,
    )
    giant_rotations = repeated_variance(
        small_switch - 2.0 * params.stc_bits,
        params.w_giant_steps - 1,
        model,
    )

    # Transform-plaintext coefficient rounding.  The paper's experiment
    # samples slots uniformly from [-B,B]; it does not place an independent
    # magnitude-B value in every polynomial coefficient.  Parseval for the
    # X/Y transforms and the prime-W embedding gives the coefficient variance
    # below.  Both plaintext-rounding paths then traverse one X and one W sum.
    input_coefficient_variance = encoded_real_coefficient_variance(params)
    transform_rounding = (
        2.0 * params.input_scale_bits
        + math.log2(input_coefficient_variance / 12.0)
        + math.log2(x_gain)
        + math.log2(w_gain)
    )
    x_encoding = (
        transform_rounding - 2.0 * params.x_scale_bits
    )
    w_encoding = (
        transform_rounding - 2.0 * params.w_scale_bits
    )
    final_rescale = _phase_rounding_log2_variance(
        dense_secret_phase_gain(params)
    )

    contributions = {
        "stc/input+first-big-switch": transformed_input,
        "stc/second-big-switch": second_conjugate,
        "stc/W-baby-switches": baby_rotations,
        "stc/W-giant-switches": giant_rotations,
        "stc/X-plaintext-encoding": x_encoding,
        "stc/W-plaintext-encoding": w_encoding,
        "stc/final-rescale": final_rescale,
    }
    return log2_sum(contributions.values()), contributions


def _masked_factor_noise(
    params: GLNoiseParams,
    model: str,
    arithmetic_mode: str,
    masked_moddown: str,
) -> tuple[float, dict[str, float]]:
    scale_bits = params.tree_scale_bits
    embedding_log2 = math.log2(_embedding_gain(params, model))
    plaintext_encoding = (
        embedding_log2 - math.log2(12.0) - 2.0 * scale_bits
    )

    candidates = params.average_candidates_per_sparse_term
    selector_coefficient_variance = params.error_stddev**2 + 1.0 / 12.0
    selector = (
        math.log2(candidates)
        + embedding_log2
        + math.log2(selector_coefficient_variance)
        - 2.0 * params.log_p
    )

    if masked_moddown == "per-candidate":
        moddown_count = candidates
    elif masked_moddown == "fused":
        moddown_count = 1.0
    else:
        raise ValueError(f"unknown masked-column modulus-down mode: {masked_moddown}")
    moddown_round = repeated_variance(
        _coefficient_to_slot_variance(
            params,
            _phase_rounding_log2_variance(dense_secret_phase_gain(params)),
            scale_bits,
            model,
        ),
        moddown_count,
        model,
    )

    hmux_switch = dd_key_switch_noise(
        params, log_q=params.log_q, destination="dense", full_ring=False
    )
    hmux_switch_count = 2 * params.hmux_radix * params.hmux_stages
    if arithmetic_mode == "legacy":
        one_hmux_switch = _coefficient_to_slot_variance(
            params, hmux_switch.log2_total_variance, scale_bits, model
        )
        hmux_eval = repeated_variance(
            _coefficient_to_slot_variance(
                params,
                hmux_switch.log2_eval_key_variance,
                scale_bits,
                model,
            ),
            hmux_switch_count,
            model,
        )
        hmux_moddown = repeated_variance(
            _coefficient_to_slot_variance(
                params,
                hmux_switch.log2_moddown_rounding_variance,
                scale_bits,
                model,
            ),
            hmux_switch_count,
            model,
        )
        hmux = repeated_variance(one_hmux_switch, hmux_switch_count, model)
    elif arithmetic_mode == "fused-dd":
        hmux_eval = repeated_variance(
            _coefficient_to_slot_variance(
                params,
                hmux_switch.log2_eval_key_variance,
                scale_bits,
                model,
            ),
            hmux_switch_count,
            model,
        )
        # Body, mask, and all radix branches remain under P*Q.  There is one
        # ModDown after their sum at each HMux stage.
        hmux_moddown = repeated_variance(
            _coefficient_to_slot_variance(
                params,
                hmux_switch.log2_moddown_rounding_variance,
                scale_bits,
                model,
            ),
            params.hmux_stages,
            model,
        )
        hmux = log2_add(hmux_eval, hmux_moddown)
    else:
        raise ValueError(f"unknown arithmetic mode: {arithmetic_mode}")

    contributions = {
        "half/candidate-encoding": plaintext_encoding,
        "half/encrypted-selectors": selector,
        "half/masked-column-moddown": moddown_round,
        "half/HMux-eval-key": hmux_eval,
        "half/HMux-moddown": hmux_moddown,
        "half/HMux-switches": hmux,
    }
    return log2_sum(contributions.values()), contributions


def _product_tree_noise(
    params: GLNoiseParams,
    model: str,
    arithmetic_mode: str,
    base_factor_log2_variance: float,
    masked_factor_log2_variance: float,
) -> tuple[ProductState, dict[str, float]]:
    gamma_factor = params.gamma / (4.0 * math.pi)
    nodes = [ProductState(base_factor_log2_variance, gamma_factor)]
    nodes.extend(
        ProductState(masked_factor_log2_variance, 1.0)
        for _ in range(params.sparse_hamming_weight)
    )
    if len(nodes) != params.factor_count:
        raise AssertionError("invalid product-tree factor count")

    contributions: dict[str, float] = {}
    for level in range(params.tree_depth):
        input_log_q = params.log_q - level * params.tree_scale_bits
        output_log_q = params.log_q - (level + 1) * params.tree_scale_bits
        if arithmetic_mode == "legacy":
            relin_log_q = output_log_q
            rescale_rounding = _tensor_rounding_log2_variance(params)
            rescale_label = "tensor-rescale"
        elif arithmetic_mode == "fused-dd":
            # Standard CKKS order: multiply modulo Q, relinearize the s^2
            # component at that same Q, and only then rescale the resulting
            # two-component ciphertext.  The final floor is r0+r1*s, not an
            # independently rounded four-component tensor.
            relin_log_q = input_log_q
            rescale_rounding = _phase_rounding_log2_variance(
                dense_secret_phase_gain(params)
            )
            rescale_label = "post-relin-rescale"
        else:
            raise ValueError(f"unknown arithmetic mode: {arithmetic_mode}")

        relin = dd_key_switch_noise(
            params,
            log_q=relin_log_q,
            destination="dense",
            full_ring=False,
        ).log2_total_variance
        rescale_slot = _coefficient_to_slot_variance(
            params, rescale_rounding, params.tree_scale_bits, model
        )
        relin_slot = _coefficient_to_slot_variance(
            params,
            relin,
            params.tree_scale_bits,
            model,
        )
        arithmetic_slot = log2_add(rescale_slot, relin_slot)
        contributions[f"tree/level-{level + 1}-{rescale_label}"] = rescale_slot
        # Stable alias used by the text and older JSON consumers.
        contributions[f"tree/level-{level + 1}-rescale"] = rescale_slot
        contributions[f"tree/level-{level + 1}-relinearization"] = relin_slot
        contributions[f"tree/level-{level + 1}-arithmetic"] = arithmetic_slot

        next_nodes: list[ProductState] = []
        for index in range(0, len(nodes), 2):
            lhs = nodes[index]
            rhs = nodes[index + 1]
            propagated = [
                lhs.log2_variance + 2.0 * math.log2(rhs.value_bound),
                rhs.log2_variance + 2.0 * math.log2(lhs.value_bound),
                lhs.log2_variance + rhs.log2_variance,
                arithmetic_slot,
            ]
            next_nodes.append(
                ProductState(
                    combine_variances(propagated, model),
                    lhs.value_bound * rhs.value_bound,
                )
            )
        nodes = next_nodes

    if len(nodes) != 1:
        raise AssertionError("product tree did not reduce to one node")
    contributions["tree/output"] = nodes[0].log2_variance
    return nodes[0], contributions


def _phase_wrap_margin_bits(
    params: GLNoiseParams,
    sparse_phase_log2_variance: float,
    tail_bound: float,
) -> float:
    # The centered phase interval is [-q0/2, q0/2).  Reserve room for the
    # encoded message q0*mu/gamma before assigning the remainder to noise.
    relative_headroom = 0.5 - params.message_bound / params.gamma
    if relative_headroom <= 0.0:
        return NEG_INF
    headroom_log2 = params.q0_bits + math.log2(relative_headroom)
    noise_bound_log2 = (
        0.5 * sparse_phase_log2_variance + math.log2(tail_bound)
    )
    return headroom_log2 - noise_bound_log2


def estimate_gl_ship(
    params: GLNoiseParams,
    *,
    model: str = "independent",
    arithmetic_mode: str = "fused-dd",
    masked_moddown: str | None = None,
    tail_bound: float = 6.0,
) -> GLNoiseEstimate:
    """Estimate one full complex GL-SHIP bootstrap.

    ``fused-dd`` keeps the hybrid-RNS operation boundaries while replacing
    decomposition/recomposition with coefficient-domain DD.  ``legacy``
    reproduces TFHEpp's former per-contribution ModDown and
    rescale-before-relinearization path.
    """
    params.validate()
    if model not in {"independent", "correlated"}:
        raise ValueError("model must be independent or correlated")
    if arithmetic_mode not in {"legacy", "fused-dd"}:
        raise ValueError("arithmetic_mode must be legacy or fused-dd")
    if masked_moddown is None:
        masked_moddown = (
            "per-candidate" if arithmetic_mode == "legacy" else "fused"
        )
    if masked_moddown not in {"per-candidate", "fused"}:
        raise ValueError("masked_moddown must be per-candidate or fused")
    if not math.isfinite(tail_bound) or tail_bound <= 0.0:
        raise ValueError("tail_bound must be positive")

    stages: dict[str, float] = {}
    stc_phase, stc_stages = _stc_phase_noise(params, model)
    stages.update(stc_stages)
    stages["stc/output-phase"] = stc_phase

    fresh_phase = 2.0 * math.log2(params.error_stddev)
    stages["half/input-fresh-phase"] = fresh_phase

    dense_to_sparse = dd_key_switch_noise(
        params,
        log_q=params.q0_bits,
        destination="sparse",
        full_ring=False,
    ).log2_total_variance
    stages["half/dense-to-sparse-switch"] = dense_to_sparse
    half_sparse_phase = log2_add(fresh_phase, dense_to_sparse)
    full_sparse_phase = log2_add(stc_phase, dense_to_sparse)
    stages["half/sparse-phase"] = half_sparse_phase
    stages["full/sparse-phase-after-StC"] = full_sparse_phase

    half_phase_to_message = half_sparse_phase + 2.0 * (
        params.gap_bits - params.q0_bits
    )
    full_phase_to_message = full_sparse_phase + 2.0 * (
        params.gap_bits - params.q0_bits
    )
    stages["half/input-phase-to-message"] = half_phase_to_message
    stages["full/input-phase-to-message"] = full_phase_to_message

    masked_factor, masked_stages = _masked_factor_noise(
        params, model, arithmetic_mode, masked_moddown
    )
    stages.update(masked_stages)
    stages["half/masked-factor"] = masked_factor

    embedding_log2 = math.log2(_embedding_gain(params, model))
    base_factor = (
        embedding_log2
        - math.log2(12.0)
        - 2.0 * params.tree_scale_bits
    )
    stages["half/base-factor-encoding"] = base_factor

    product, product_stages = _product_tree_noise(
        params, model, arithmetic_mode, base_factor, masked_factor
    )
    stages.update(product_stages)

    output_switch = dd_key_switch_noise(
        params,
        log_q=params.output_log_q,
        destination="dense",
        full_ring=False,
    ).log2_total_variance
    output_switch_slot = _coefficient_to_slot_variance(
        params,
        output_switch,
        params.tree_scale_bits,
        model,
    )
    channel = combine_variances(
        [product.log2_variance, product.log2_variance, output_switch_slot],
        model,
    )
    stages["half/output-conjugation-switch"] = output_switch_slot
    stages["half/one-channel"] = channel

    # Two Gaussian channels are combined as C0 + I*C1.  A direct half
    # bootstrap starts with a freshly encrypted coefficient ciphertext.  The
    # full bootstrap instead feeds each slice with grouped-StC output noise.
    half_he = combine_variances(
        [channel, channel, half_phase_to_message, half_phase_to_message], model
    )
    full_slice_he = combine_variances(
        [channel, channel, full_phase_to_message, full_phase_to_message], model
    )
    y_gain = _sum_gain(params.n, model)
    full_he = full_slice_he + math.log2(y_gain)
    stages["half/complex-output"] = half_he
    stages["full/per-slice-complex-output"] = full_slice_he
    stages["full/Y-reassembly"] = full_he

    half_he_bound = tail_bound * 2.0 ** (0.5 * half_he)
    full_he_bound = tail_bound * 2.0 ** (0.5 * full_he)

    # SHIP, Equation (1)/(17): |gamma/(2*pi) sin(2*pi*m/q0)-mu|.
    # A direct half bootstrap accepts an arbitrary bounded coefficient-domain
    # message, so retain the deterministic cubic bound.  In the paper's full
    # experiment, however, uniform [-B,B] slots have first been transformed to
    # Y coefficients.  Their per-channel variance is B^2/(3n); using the sixth
    # Gaussian moment gives the matching average-case cubic-error variance.
    real_sine_bound = (
        (2.0 * math.pi) ** 2
        * params.message_bound**3
        / (6.0 * params.gamma**2)
    )
    real_quantization_bound = params.gamma / (2.0 * (2.0**params.q0_bits))
    complex_factor = math.sqrt(2.0)
    half_sine = complex_factor * real_sine_bound
    half_quantization = complex_factor * real_quantization_bound

    cubic_coefficient = (2.0 * math.pi) ** 2 / (6.0 * params.gamma**2)
    full_phase_message_variance = params.message_bound**2 / (3.0 * params.n)
    real_sine_variance = (
        15.0 * cubic_coefficient**2 * full_phase_message_variance**3
    )
    full_sine_variance = 2.0 * real_sine_variance * y_gain
    full_sine = tail_bound * math.sqrt(full_sine_variance)

    # GLEncode rounds each real polynomial coefficient.  Reconstructing one
    # complex matrix slot traverses the full real GL dimension.
    full_quantization_variance = (
        params.full_ring_degree
        * params.gamma**2
        / (12.0 * (2.0 ** (2 * params.q0_bits)))
    )
    full_quantization = tail_bound * math.sqrt(full_quantization_variance)

    half_total = half_he_bound + half_sine + half_quantization
    full_total = full_he_bound + full_sine + full_quantization
    half_precision = -math.log2(half_total)
    full_precision = -math.log2(full_total)

    half_phase_wrap_margin = _phase_wrap_margin_bits(
        params, half_sparse_phase, tail_bound
    )
    full_phase_wrap_margin = _phase_wrap_margin_bits(
        params, full_sparse_phase, tail_bound
    )

    return GLNoiseEstimate(
        params=params,
        model=model,
        arithmetic_mode=arithmetic_mode,
        masked_moddown=masked_moddown,
        tail_bound=tail_bound,
        stages=stages,
        fresh_phase_log2_variance=fresh_phase,
        stc_phase_log2_variance=stc_phase,
        half_sparse_phase_log2_variance=half_sparse_phase,
        full_sparse_phase_log2_variance=full_sparse_phase,
        half_phase_wrap_margin_bits=half_phase_wrap_margin,
        full_phase_wrap_margin_bits=full_phase_wrap_margin,
        half_he_log2_variance=half_he,
        full_he_log2_variance=full_he,
        half_he_error_bound=half_he_bound,
        full_he_error_bound=full_he_bound,
        half_sine_error_bound=half_sine,
        full_sine_error_bound=full_sine,
        half_quantization_error_bound=half_quantization,
        full_quantization_error_bound=full_quantization,
        half_total_error_bound=half_total,
        full_total_error_bound=full_total,
        half_precision_bits=half_precision,
        full_precision_bits=full_precision,
    )


def format_log2_variance(log2_variance: float) -> str:
    if log2_variance == NEG_INF:
        return "-inf"
    return f"{log2_variance:.2f}"


def estimate_as_dict(estimate: GLNoiseEstimate) -> dict[str, object]:
    params = estimate.params
    return {
        "profile": params.tag,
        "model": estimate.model,
        "arithmetic_mode": estimate.arithmetic_mode,
        "masked_moddown": estimate.masked_moddown,
        "tail_bound": estimate.tail_bound,
        "parameters": {
            "N": params.ring_degree,
            "n": params.n,
            "p": params.p,
            "phi_p": params.phi,
            "log_q": params.log_q,
            "log_p": params.log_p,
            "log_pq": params.log_pq,
            "q0_bits": params.q0_bits,
            "gap_bits": params.gap_bits,
            "gamma": params.gamma,
            "stc_bits": params.stc_bits,
            "x_scale_bits": params.x_scale_bits,
            "w_scale_bits": params.w_scale_bits,
            "input_log_q": params.input_log_q,
            "input_scale_bits": params.input_scale_bits,
            "tree_scale_bits": params.tree_scale_bits,
            "tree_depth": params.tree_depth,
            "output_log_q": params.output_log_q,
            "required_output_log_q": params.required_output_log_q,
            "output_depth_margin_bits": params.output_depth_margin_bits,
            "tree_scale_headroom_bits": params.tree_scale_headroom_bits,
            "outside_multiplicative_depth": params.outside_multiplicative_depth,
            "outside_scale_bits": params.outside_scale_bits,
            "paper_dnum": params.dnum,
            "theta": params.theta,
            "window_width": params.window_width,
            "masked_column_count": params.masked_column_count,
            "average_candidates_per_sparse_term": (
                params.average_candidates_per_sparse_term
            ),
            "sparse_hamming_weight": params.sparse_hamming_weight,
            "hmux_radix": params.hmux_radix,
            "hmux_stages": params.hmux_stages,
            "w_baby_step": params.w_baby_step,
            "w_giant_steps": params.w_giant_steps,
            "primary_bit": params.primary_bit,
            "primary_cover_bits": params.primary_cover_bits,
            "bbar_bit": params.bbar_bit,
            "bbar_cover_bits": params.bbar_cover_bits,
            "storage_bits": params.storage_bits,
            "storage_margin_bits": params.storage_margin_bits,
            "security_limit_log_pq": params.security_limit_log_pq,
            "security_margin_bits": params.security_margin_bits,
            "error_stddev": params.error_stddev,
            "dense_secret_variance": params.dense_secret_variance,
            "message_bound": params.message_bound,
        },
        "noise": {
            "fresh_phase_log2_variance": estimate.fresh_phase_log2_variance,
            "stc_phase_log2_variance": estimate.stc_phase_log2_variance,
            "half_sparse_phase_log2_variance": (
                estimate.half_sparse_phase_log2_variance
            ),
            "full_sparse_phase_log2_variance": (
                estimate.full_sparse_phase_log2_variance
            ),
            "half_phase_wrap_margin_bits": (
                estimate.half_phase_wrap_margin_bits
            ),
            "full_phase_wrap_margin_bits": (
                estimate.full_phase_wrap_margin_bits
            ),
            "phase_wrap_margin_bits": estimate.phase_wrap_margin_bits,
            "half_he_log2_variance": estimate.half_he_log2_variance,
            "full_he_log2_variance": estimate.full_he_log2_variance,
            "half_he_error_bound": estimate.half_he_error_bound,
            "full_he_error_bound": estimate.full_he_error_bound,
            "half_sine_error_bound": estimate.half_sine_error_bound,
            "full_sine_error_bound": estimate.full_sine_error_bound,
            "half_quantization_error_bound": estimate.half_quantization_error_bound,
            "full_quantization_error_bound": estimate.full_quantization_error_bound,
            "half_total_error_bound": estimate.half_total_error_bound,
            "full_total_error_bound": estimate.full_total_error_bound,
            "half_precision_bits": estimate.half_precision_bits,
            "full_precision_bits": estimate.full_precision_bits,
            "paper_precision_bits": params.paper_precision_bits,
            "precision_delta_from_paper": estimate.precision_delta_from_paper,
        },
        "stages_log2_variance": dict(estimate.stages),
    }
