#!/usr/bin/env python3
"""
Average-case BFV invariant-noise estimator.

This module implements the formulas from `600.pdf`:

    Beatrice Biasioli, Chiara Marcolla, Marco Calderini, Johannes Mono,
    "Improving and Automating BFV Parameters Selection: An Average-Case
    Approach", 2024.

The estimator works in log2-space so large candidate moduli such as 2^900 do
not overflow Python floats.  The formulas are the independent-ciphertext model
from Sections 4 and 5 of the paper.  Polynomial evaluation of one ciphertext is
not independent; use the output as an optimistic screening estimate until a
dependent-circuit model is added.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Iterable

from scipy.special import erfc, erfcinv


NEG_INF = float("-inf")


def log2_add(*xs: float) -> float:
    """Return log2(sum(2**x for x in xs))."""
    finite = [x for x in xs if x != NEG_INF]
    if not finite:
        return NEG_INF
    m = max(finite)
    return m + math.log2(sum(2.0 ** (x - m) for x in finite))


def log2_sum(xs: Iterable[float]) -> float:
    return log2_add(*list(xs))


@dataclass(frozen=True)
class CorrectionFit:
    """Heuristic 1 fit f(i) = -exp(alpha - beta*i - gamma*i^2) + delta."""

    alpha: float
    beta: float
    gamma: float
    delta: float

    def value(self, i: int) -> float:
        if i <= 1:
            return 1.0
        v = self.delta - math.exp(self.alpha - self.beta * i - self.gamma * i * i)
        return max(1.0, v)

    def log2_value(self, i: int) -> float:
        return math.log2(self.value(i))


# Appendix B / Table 1 in 600.pdf, chi_s = U3 (ternary {-1,0,1}).
U3_CORRECTION_FITS = {
    12: CorrectionFit(alpha=2.8732, beta=0.0160, gamma=0.0049, delta=19.1895),
    13: CorrectionFit(alpha=2.9644, beta=0.0196, gamma=0.0046, delta=20.4747),
    14: CorrectionFit(alpha=2.9578, beta=0.0386, gamma=0.0032, delta=19.5755),
    15: CorrectionFit(alpha=2.9765, beta=0.0197, gamma=0.0043, delta=20.7760),
}


def correction_fit_for_nbit(nbit: int) -> CorrectionFit:
    if nbit in U3_CORRECTION_FITS:
        return U3_CORRECTION_FITS[nbit]
    # Conservative fallback: use the largest asymptote available.  This keeps
    # exploratory sweeps working while making unsupported nbits visible.
    return U3_CORRECTION_FITS[max(U3_CORRECTION_FITS)]


@dataclass(frozen=True)
class BFVParams:
    nbit: int
    t: int
    q_bits: int
    secret_variance: float = 2.0 / 3.0
    encryption_u_variance: float = 2.0 / 3.0
    error_log2_std: float = math.log2(3.19)
    key_switch: str = "hybrid-rns"
    rns_digits: int = 4
    hybrid_omega: int = 4
    fresh: str = "public"
    correction_fit: CorrectionFit | None = None
    tag: str = "bfv"

    @property
    def n(self) -> int:
        return 1 << self.nbit

    @property
    def log2_error_variance(self) -> float:
        return 2.0 * self.error_log2_std

    @property
    def fit(self) -> CorrectionFit:
        return self.correction_fit or correction_fit_for_nbit(self.nbit)

    def with_q_bits(self, q_bits: int) -> "BFVParams":
        return replace(self, q_bits=q_bits)


@dataclass(frozen=True)
class NoiseState:
    """Log2 variance of one invariant-noise coefficient and degree in s."""

    log2_variance: float
    degree: int
    label: str = ""

    def add(self, other: "NoiseState", label: str = "add") -> "NoiseState":
        return NoiseState(
            log2_add(self.log2_variance, other.log2_variance),
            max(self.degree, other.degree),
            label,
        )


def fresh(params: BFVParams) -> NoiseState:
    """Fresh encryption variance, Proposition 3 or symmetric BFV variant."""
    log_t2_over_q2 = 2.0 * math.log2(params.t) - 2.0 * params.q_bits
    log_ve = params.log2_error_variance

    if params.fresh == "symmetric":
        inside = log_ve
    elif params.fresh == "public":
        # 1/12 + n*Ve*Vu + Ve + n*Ve*Vs
        inside = log2_add(
            -math.log2(12.0),
            math.log2(params.n) + log_ve + math.log2(params.encryption_u_variance),
            log_ve,
            math.log2(params.n) + log_ve + math.log2(params.secret_variance),
        )
    else:
        raise ValueError(f"unknown fresh encryption mode: {params.fresh}")

    return NoiseState(log_t2_over_q2 + inside, degree=1, label="fresh")


def const_mul(state: NoiseState, params: BFVParams, label: str = "const") -> NoiseState:
    """Average polynomial constant multiplication, Proposition 4."""
    factor = ((params.t * params.t - 1.0) * params.n) / 12.0
    return NoiseState(state.log2_variance + math.log2(factor), state.degree, label)


def scalar_mul_unsigned_average(
    state: NoiseState, params: BFVParams, label: str = "scalar-unsigned"
) -> NoiseState:
    """Average scalar multiplication by coefficients sampled from [0, t)."""
    factor_log2 = (
        math.log2(params.t - 1.0) + math.log2(2.0 * params.t - 1.0) - math.log2(6.0)
    )
    return NoiseState(state.log2_variance + factor_log2, state.degree, label)


def scalar_mul_centered_average(
    state: NoiseState, params: BFVParams, label: str = "scalar-centered"
) -> NoiseState:
    """Average scalar multiplication by centered coefficients modulo t."""
    factor_log2 = math.log2(params.t * params.t - 1.0) - math.log2(12.0)
    return NoiseState(state.log2_variance + factor_log2, state.degree, label)


def mod_switch(state: NoiseState, params: BFVParams, q_prime_bits: int | None = None) -> NoiseState:
    """Modulo switch variance, Proposition 5."""
    if q_prime_bits is None:
        q_prime_bits = params.q_bits
    add_log = (
        2.0 * math.log2(params.t)
        + math.log2(1.0 + params.n * params.secret_variance)
        - math.log2(12.0)
        - 2.0 * q_prime_bits
    )
    return NoiseState(log2_add(state.log2_variance, add_log), state.degree, "mod-switch")


def key_switch_variance_log2(params: BFVParams) -> float:
    """Invariant coefficient variance added by relinearization/key switching."""
    n = params.n
    vs = params.secret_variance
    log_ve = params.log2_error_variance

    if params.key_switch == "ghs":
        a0 = log2_add(math.log2(n) + log_ve, 0.0)
    elif params.key_switch == "ghs-rns":
        a0 = log2_add(math.log2((params.rns_digits + 2) * n) + log_ve, 0.0)
    elif params.key_switch == "hybrid":
        a0 = log2_add(math.log2(params.hybrid_omega * n) + log_ve, 0.0)
    elif params.key_switch == "hybrid-rns":
        a0 = log2_add(
            math.log2((params.rns_digits + 2 * params.hybrid_omega) * n) + log_ve,
            0.0,
        )
    else:
        raise ValueError(f"unsupported key switch variant: {params.key_switch}")

    # a0 term plus a1*s term.  a1's coefficient from the paper is 1.
    invariant_factor = log2_add(a0, math.log2(n * vs))
    return (
        2.0 * math.log2(params.t)
        - math.log2(12.0)
        - 2.0 * params.q_bits
        + invariant_factor
    )


def key_switch(state: NoiseState, params: BFVParams) -> NoiseState:
    return NoiseState(
        log2_add(state.log2_variance, key_switch_variance_log2(params)),
        state.degree,
        "key-switch",
    )


def mul(lhs: NoiseState, rhs: NoiseState, params: BFVParams, relin: bool = True) -> NoiseState:
    """Independent-ciphertext multiplication estimate, Theorem 1."""
    c = (
        2.0 * math.log2(params.t)
        + 2.0 * math.log2(params.n)
        + math.log2(params.secret_variance)
        - math.log2(12.0)
    )
    inner = log2_add(
        lhs.log2_variance + params.fit.log2_value(lhs.degree + 1),
        rhs.log2_variance + params.fit.log2_value(rhs.degree + 1),
    )
    out = NoiseState(c + inner, lhs.degree + rhs.degree, "mul")
    return key_switch(out, params) if relin else out


def floor_pow2(x: int) -> int:
    p = 1
    while 2 * p <= x:
        p *= 2
    return p


def find_bsgs_params(degree: int) -> tuple[int, int]:
    if degree <= 0:
        return 1, 0
    best_k, best_m, best_cost = 1, 0, degree + 1
    for k in range(1, degree + 1):
        m = 0
        reach = k
        while reach < degree:
            reach *= 2
            m += 1
        cost = (k - 1) + (m - 1 if m > 0 else 0) + m
        if cost < best_cost:
            best_k, best_m, best_cost = k, m, cost
        if k > 2 * best_k:
            break
    return best_k, best_m


def estimate_polyeval_bsgs(
    x: NoiseState,
    degree: int,
    params: BFVParams,
    *,
    scalar_mode: str = "none",
) -> NoiseState:
    """
    Estimate TFHEpp's baby-step/giant-step PolyEval shape.

    `scalar_mode=none` ignores cleartext coefficient amplification.  This is
    useful for first-pass depth/noise-budget screening of high-degree digit
    extraction polynomials.

    `scalar_mode=unsigned-average` models TFHEpp's current unsigned scalar
    coefficient multiplication.  `centered-average` models the same operation if
    polynomial coefficients are centered before scalar multiplication.
    `poly-average` applies Proposition 4, modeling multiplication by a random
    plaintext polynomial rather than a scalar coefficient.
    """
    if degree <= 0:
        return NoiseState(NEG_INF, 0, "constant")
    if scalar_mode == "average":
        scalar_mode = "unsigned-average"
    if scalar_mode not in {"none", "unsigned-average", "centered-average", "poly-average"}:
        raise ValueError(
            "scalar_mode must be 'none', 'unsigned-average', 'centered-average', or 'poly-average'"
        )

    k, m = find_bsgs_params(degree)

    baby: list[NoiseState | None] = [None for _ in range(k + 1)]
    baby[1] = x

    def baby_state(i: int) -> NoiseState:
        state = baby[i]
        if state is None:
            raise RuntimeError(f"missing baby-step x^{i}")
        return state

    for i in range(2, k + 1):
        a = floor_pow2(i - 1)
        b = i - a
        baby[i] = mul(baby_state(a), baby_state(b), params)

    giant: list[NoiseState] = []
    if m > 0:
        giant.append(baby_state(k))
        for j in range(1, m):
            giant.append(mul(giant[j - 1], giant[j - 1], params))

    def scale_term(term: NoiseState) -> NoiseState:
        if scalar_mode == "unsigned-average":
            return scalar_mul_unsigned_average(term, params)
        if scalar_mode == "centered-average":
            return scalar_mul_centered_average(term, params)
        if scalar_mode == "poly-average":
            return const_mul(term, params)
        return term

    def eval_recursive(length: int, level: int) -> NoiseState:
        if length <= 1:
            return NoiseState(NEG_INF, 0, "constant")
        if level == 0:
            terms = [scale_term(baby_state(i)) for i in range(1, min(length, k + 1))]
            if not terms:
                return NoiseState(NEG_INF, 0, "constant")
            return NoiseState(
                log2_sum(t.log2_variance for t in terms),
                max(t.degree for t in terms),
                "baby-sum",
            )

        split = k
        for _ in range(level - 1):
            split *= 2
        split = min(split, length)

        lower = eval_recursive(split, level - 1)
        if split >= length:
            return lower
        upper = eval_recursive(length - split, level - 1)
        prod = mul(upper, giant[level - 1], params)
        return lower.add(prod, "recursive-sum")

    return eval_recursive(degree + 1, m)


def correctness_d_for_failure(n: int, failure_log2: float) -> float:
    """Return D such that n*erfc(D) is at most 2**failure_log2."""
    target = 2.0 ** failure_log2 / n
    if target <= 0:
        raise ValueError("failure target underflowed")
    return float(erfcinv(target))


def log2_correctness_threshold(d: float) -> float:
    """Correct decryption threshold V <= 1/(8*D^2)."""
    return -math.log2(8.0 * d * d)


def failure_log2_from_variance(n: int, log2_variance: float) -> float:
    if log2_variance == NEG_INF:
        return NEG_INF
    if log2_variance > 2048.0:
        return math.log2(n)
    if log2_variance < -2048.0:
        return NEG_INF
    d = 1.0 / (2.0 * math.sqrt(2.0) * (2.0 ** (0.5 * log2_variance)))
    p = n * erfc(d)
    if p == 0.0:
        return NEG_INF
    return math.log2(p)


def format_log2(x: float) -> str:
    if x == NEG_INF:
        return "-inf"
    return f"{x:.2f}"
