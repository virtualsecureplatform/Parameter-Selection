"""Modulus budget and average-case noise screen for TFHEpp's full-slot TUTTI.

TUTTI (ePrint 2026/2848) evaluates the decryption equation on the unit
circle: one accumulator per block of a block-binary secret, each a sum of
public root encodings multiplied by encrypted indicators (AccMul), followed
by a product tree and a conjugation. TFHEpp runs it on power-of-two moduli:

* AccMul rescales by `plain_bits + indicator_scale_bits - scale_bits`, the
  role of the paper's auxiliary prime. A power-of-two modulus lets the
  indicators carry a larger scale than the tree: each extra bit halves their
  noise and costs one bit of modulus, where a larger tree scale would cost
  one bit per level;
* the tree relinearizes each tensor before rescaling, so its `scale_bits`
  suppress the double-decomposition switch as a special modulus would;
* dense-key automorphisms (S2C, conjugation) are lifted by `auxiliary_bits`;
* the switch into the sparse key is a plain gadget switch at the base
  modulus, so that key is never published above `base_bits`.

Variances are log2, relative to a unit slot value. FFT error is excluded and
the screen is not a correctness or security certificate.
"""

from dataclasses import dataclass, asdict, replace
import functools
import math

from .gl import ProductState, combine_variances, log2_add, NEG_INF
from .narrow_ckks import NarrowParams, s2c_variance, switch_slot_variance

# Largest integer log q at which the cheapest attack of the paper's model
# still costs 128 bits (their Tables 1 and 2). The sparse caps come from the
# rotated primal hybrid of Hou, Jiang and Ogilvie and are external inputs:
# this repository does not reproduce that attack.
PAPER_DENSE_CAP = {4096: 106, 8192: 214, 16384: 430, 32768: 868, 65536: 1746}
PAPER_DENSE_CAP_WEIGHT_512 = {8192: 212, 16384: 428, 32768: 862}
PAPER_SPARSE_CAP = {(2048, 64): 27, (4096, 64): 51, (8192, 64): 102,
                    (16384, 64): 209, (8192, 32): 27, (16384, 32): 59}


def truncated_normal_rms(sigma):
    """Integer RMS of a continuous normal sample truncated toward zero.

    TFHEpp converts `sigma * N(0, 1)` to an integer by truncation, so the
    integer error is smaller than the configured deviation.
    """
    if not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be positive and finite")
    tail = lambda x: math.erfc(x / (sigma * math.sqrt(2)))  # P(|X| >= x)
    limit = max(64, math.ceil(40 * sigma))
    return math.sqrt(sum((2 * k - 1) * tail(k) for k in range(1, limit)))


@functools.lru_cache(maxsize=None)
def minimax_sine_scale(gap, message_bound):
    """Constant c for which c*g*sin(v/g) - v equioscillates on |v| <= bound.

    g = gap / (2 pi). Returns (c, largest absolute error).
    """
    g = gap / (2 * math.pi)
    x = message_bound / g
    if not 0 < x < math.pi / 2:
        raise ValueError("message bound must stay inside the monotone sine range")
    interior = lambda c: g * (math.sqrt(c * c - 1) - math.acos(1 / c))
    endpoint = lambda c: message_bound - c * g * math.sin(x)
    low, high = 1.0, x / math.sin(x)
    for _ in range(200):
        mid = (low + high) / 2
        if interior(mid) < endpoint(mid):
            low = mid
        else:
            high = mid
    return high, max(interior(high), endpoint(high))


@dataclass(frozen=True)
class TuttiParams:
    n: int = 8192
    blocks: int = 64
    stride: int = 4
    base_bits: int = 24
    input_scale_bits: int = 20
    s2c_bits: int = 17
    scale_bits: int = 22
    indicator_scale_bits: int = 22
    plain_bits: int = 18
    output_bits: int = 61
    output_scale_bits: int = 20
    auxiliary_bits: int = 20
    primary_bits: int = 12
    bbar_bits: int = 8
    sparse_digit_bits: int = 1
    sigma: float = 3.2
    dense_weight: int = 512
    message_bound: float = 1.0
    minimax_sine: bool = True
    # Boolean batch refresh of TFHE TLWEs (TFHEpp ckks/tutti_boolean.hpp): the
    # input is packed into coefficients under the dense key at 2^48, so there
    # is no opening S2C; a closing S2C at the output modulus returns the
    # slots to coefficients at 2^32 for the switch back to TFHE.
    boolean_refresh: bool = False

    @property
    def depth(self):
        return self.blocks.bit_length() - 1

    @property
    def block_length(self):
        return self.n // self.blocks

    @property
    def candidates(self):
        return self.block_length // self.stride

    @property
    def indicators(self):
        return (self.blocks - 1) * self.candidates

    @property
    def tree_output_bits(self):
        return self.output_bits + self.scale_bits - self.output_scale_bits

    @property
    def factor_bits(self):
        return self.tree_output_bits + self.depth * self.scale_bits

    @property
    def accmul_drop_bits(self):
        return self.plain_bits + self.indicator_scale_bits - self.scale_bits

    @property
    def selector_bits(self):
        return self.factor_bits + self.accmul_drop_bits

    @property
    def gap(self):
        return 2.0 ** (self.base_bits - self.input_scale_bits)

    @property
    def dense_norm2(self):
        return self.dense_weight if self.dense_weight else self.n * 2 / 3

    @property
    def sparse_degree(self):
        return self.n // self.stride

    def sine(self):
        """(scale constant, deterministic sine error bound)."""
        if self.minimax_sine:
            return minimax_sine_scale(self.gap, self.message_bound)
        ratio = 2 * math.pi * self.message_bound / self.gap
        return 1.0, self.gap / (2 * math.pi) * (ratio - math.sin(ratio))

    def validate(self):
        power = lambda x: x > 0 and not x & (x - 1)
        if not power(self.n) or self.n < 16:
            raise ValueError("N must be a power of two of at least 16")
        if not power(self.blocks) or not 2 <= self.blocks <= self.n:
            raise ValueError("blocks must be a power of two between 2 and N")
        if not power(self.stride) or self.stride > self.block_length:
            raise ValueError("stride must be a power of two dividing the block length")
        if min(self.scale_bits, self.plain_bits, self.s2c_bits, self.primary_bits,
               self.bbar_bits, self.sparse_digit_bits, self.output_scale_bits) <= 0:
            raise ValueError("scales and digit widths must be positive")
        if self.auxiliary_bits < 0:
            raise ValueError("auxiliary bits must be nonnegative")
        if not 0 < self.input_scale_bits < self.base_bits:
            raise ValueError("input scale must lie below the base modulus")
        if not self.output_scale_bits <= self.scale_bits:
            raise ValueError("output scale cannot exceed the tree scale")
        if self.indicator_scale_bits < self.scale_bits:
            raise ValueError("indicator scale cannot be below the tree scale")
        if not math.isfinite(self.sigma) or self.sigma < 0:
            raise ValueError("sigma must be finite and nonnegative")
        if self.dense_weight < 0 or self.dense_weight > self.n:
            raise ValueError("invalid dense Hamming weight")
        if not 0 < self.message_bound < self.gap / 4:
            raise ValueError("message bound must stay inside the monotone sine range")
        amplitude = math.ceil(math.log2(self.sine()[0] * self.gap / (4 * math.pi) + 1))
        if self.output_bits < self.output_scale_bits + max(amplitude, 0) + 2:
            raise ValueError("output modulus cannot hold the unit-circle product")
        if self.output_bits <= self.base_bits + self.s2c_bits:
            raise ValueError("output modulus does not permit a repeatable refresh")


def budget(p):
    """Moduli at which each secret is exposed, and what the output retains."""
    p.validate()
    input_bits = p.base_bits + p.s2c_bits
    application = p.output_bits - input_bits
    mult_drop = 2 * p.output_scale_bits - p.input_scale_bits
    dense = {
        "indicators": p.selector_bits,
        "tree_relinearization": p.factor_bits,
        "output_conjugation": p.output_bits + p.auxiliary_bits,
        "s2c_rotations": input_bits + p.auxiliary_bits,
        "application_relinearization": p.output_bits,
    }
    if p.boolean_refresh:
        if p.output_bits != 32 + p.s2c_bits:
            raise ValueError("the closing S2C must end at the 32-bit TLWE modulus")
        del dense["application_relinearization"]
        dense["s2c_rotations"] = p.output_bits + p.auxiliary_bits
        dense["tlwe_packing"] = 48
        input_bits, application, mult_drop = p.base_bits, 0, 0
    return {
        "input_bits": input_bits, "output_bits": p.output_bits,
        "tree_output_bits": p.tree_output_bits, "factor_bits": p.factor_bits,
        "selector_bits": p.selector_bits,
        "dense_key_exposures": dense, "maximum_dense_key_bits": max(dense.values()),
        "sparse_key_exposure_bits": p.base_bits,
        "sparse_degree": p.sparse_degree, "sparse_blocks": p.blocks,
        "application_bits": application,
        "multiplication_drop_bits": mult_drop,
        "multiplications_returning_to_input": application // mult_drop if mult_drop > 0 else 0,
        "one_multiplication_returns_to_input": mult_drop > 0 and application == mult_drop,
        "other_published_keys": ["dense CKKS key under the TFHE ring key at 2^32 (subring switch)"]
                                if p.boolean_refresh else [],
    }


def published_caps(p):
    """Compare the exposures with the paper's published 128-bit caps."""
    exposure = budget(p)
    dense_table = PAPER_DENSE_CAP_WEIGHT_512 if p.dense_weight == 512 else PAPER_DENSE_CAP
    dense_cap = dense_table.get(p.n) if p.dense_weight in (0, 512) else None
    sparse_cap = PAPER_SPARSE_CAP.get((p.sparse_degree, p.blocks))
    return {
        "source": "ePrint 2026/2848 Tables 1 and 2; sigma 3.2, lattice-estimator word operations",
        "dense_cap_bits": dense_cap,
        "dense_within_cap": None if dense_cap is None else
            exposure["maximum_dense_key_bits"] <= dense_cap,
        "sparse_row": [p.sparse_degree, p.blocks],
        "sparse_cap_bits": sparse_cap,
        "sparse_within_cap": None if sparse_cap is None else p.base_bits <= sparse_cap,
        "integer_noise_rms_matches_model": p.sigma >= 3.2,
    }


def linear_variance(p, model="independent"):
    """Coefficient error after the opening S2C and the sparse switch.

    The S2C pipeline is lifted by the auxiliary scale and divided under the
    dense key; the sparse switch then runs at the base modulus.
    """
    aux = p.auxiliary_bits
    lifted = NarrowParams(n=p.n, sigma=p.sigma, source_bits=p.base_bits + aux,
                          input_scale_bits=p.input_scale_bits + aux, s2c_bits=p.s2c_bits,
                          bbar_bits=p.bbar_bits, message_bound=p.message_bound)
    _, parts = s2c_variance(lifted, model)
    del parts["encapsulation"]
    n = p.n
    baby = 1 << ((n.bit_length() - 2) // 2)
    giants = (n // 2 + baby - 1) // baby
    scale = 2 * p.input_scale_bits
    parts["fresh"] = NEG_INF if p.sigma == 0 else math.log2(n * p.sigma**2) - scale
    parts["rescale"] = (math.log2(2 * giants * (1 + p.dense_norm2) / 12)
                        - scale - 2 * aux)
    if aux:
        parts["auxiliary_division"] = math.log2((1 + p.dense_norm2) / 12) - scale
    rows = (p.base_bits + p.sparse_digit_bits - 1) // p.sparse_digit_bits
    parts["sparse_switch"] = NEG_INF if p.sigma == 0 else (
        math.log2(rows * n) + log2_add(2 * p.sparse_digit_bits, 1) - math.log2(12)
        + 2 * math.log2(p.sigma) - scale)
    return combine_variances(parts.values(), model), parts


def estimate(p, model="independent"):
    p.validate()
    n, dense = p.n, p.dense_norm2
    scale = 2 * p.scale_bits
    sine_scale, sine_error = p.sine()
    beta = sine_scale * p.gap / (4 * math.pi)
    delta = beta ** (1 / p.blocks)
    rounding = math.log2(n * (1 + dense) / 12) - scale
    fresh = NEG_INF if p.sigma == 0 else (
        math.log2(p.candidates * n) + 2 * math.log2(p.sigma * delta)
        - 2 * p.indicator_scale_bits)
    encoding = math.log2(n * delta**2 / 12) - 2 * p.plain_bits
    factor_v = combine_variances((fresh, encoding, rounding), model)
    public_v = math.log2(n / 12) - scale
    nodes = [ProductState(public_v, delta)] + [ProductState(factor_v, delta)] * (p.blocks - 1)
    levels = []
    for level in range(p.depth):
        q = p.factor_bits - level * p.scale_bits
        relin = switch_slot_variance(replace(_switch(p), scale_bits=scale), q, p.primary_bits)
        arithmetic = combine_variances((rounding, relin), model)
        merged = []
        for a, b in zip(nodes[::2], nodes[1::2]):
            variance = combine_variances((
                a.log2_variance + 2 * math.log2(b.value_bound),
                b.log2_variance + 2 * math.log2(a.value_bound),
                a.log2_variance + b.log2_variance, arithmetic), model)
            merged.append(ProductState(variance, a.value_bound * b.value_bound))
        nodes = merged
        levels.append({"input_bits": q, "output_bits": q - p.scale_bits,
                       "relin_log2_variance": None if relin == NEG_INF else relin,
                       "node_log2_variance": nodes[-1].log2_variance})
    tree_v = nodes[0].log2_variance
    out_scale = 2 * p.output_scale_bits
    out_rounding = math.log2(n * (1 + dense) / 12) - out_scale
    drop_v = out_rounding if p.output_scale_bits < p.scale_bits else NEG_INF
    conjugation = switch_slot_variance(
        replace(_switch(p), scale_bits=p.output_scale_bits + p.auxiliary_bits),
        p.output_bits + p.auxiliary_bits, p.bbar_bits)
    if p.auxiliary_bits:
        conjugation = log2_add(conjugation, out_rounding)
    # out + conj(out) doubles the real part of the tree error; only the real
    # part of the conjugation switch error reaches a real slot.
    circle_v = combine_variances((tree_v + 1, drop_v + 1, conjugation - 1), model)
    linear_v, linear_parts = linear_variance(p, model)
    total_v = combine_variances((circle_v, linear_v), model)
    rms = 2 ** (total_v / 2)
    slots = n // 2
    expected_max = rms * math.sqrt(2 * math.log(2 * slots)) + sine_error
    limbs = (p.selector_bits + 63) // 64
    digit_bits = min(60, 118 - p.plain_bits - 2 - (n.bit_length() - 1)
                     - max(1, (p.candidates - 1).bit_length()))
    digits = (p.selector_bits + digit_bits - 1) // digit_bits
    return {
        "scope": __doc__.splitlines()[0], "model": model, "security_validated": False,
        "budget": budget(p),
        "gap": p.gap, "sine_scale": sine_scale, "sine_error_bound": sine_error,
        "unit_circle_amplitude": beta, "per_factor_amplitude": delta,
        "factor_log2_variances": {"indicator_noise": None if fresh == NEG_INF else fresh,
                                  "root_encoding": encoding, "accmul_rescale": rounding,
                                  "total": factor_v},
        "levels": levels,
        "tree_log2_variance": tree_v,
        "output_drop_log2_variance": None if drop_v == NEG_INF else drop_v,
        "conjugation_log2_variance": conjugation,
        "unit_circle_real_log2_variance": circle_v,
        "linear_log2_variances": {k: None if v == NEG_INF else v
                                  for k, v in linear_parts.items()},
        "linear_log2_variance": linear_v,
        "real_refresh_rms": rms,
        "rms_precision_bits": -total_v / 2,
        "expected_worst_slot_error": expected_max,
        "expected_worst_slot_precision_bits": -math.log2(expected_max),
        "eight_rms_plus_sine_screen": 8 * rms + sine_error,
        "screen_note": "Expected maximum assumes independent Gaussian slot errors; "
                       "neither figure is a tail or failure bound",
        "accmul": {"indicators": p.indicators, "candidates_per_block": p.candidates,
                   "selector_digit_bits": digit_bits, "selector_digits": digits,
                   "raw_indicator_bytes": p.indicators * 2 * n * limbs * 8,
                   "ntt_cache_bytes": p.indicators * 2 * digits * 4 * n * 4,
                   "forward_ntts": p.indicators * 4,
                   "inverse_ntts": (p.blocks - 1) * 2 * digits * 4},
        "excluded": ["FFT error", "security"],
    }


def _switch(p):
    return NarrowParams(n=p.n, sigma=p.sigma, scale_bits=p.scale_bits, bbar_bits=p.bbar_bits,
                        primary_bits=p.primary_bits)


def search(base, maximum_dense_bits, application_levels=1, gap_bits=range(3, 8),
           s2c_bits=range(14, 21), plain_bits=range(14, 25), indicator_gain_bits=range(0, 11),
           model="independent"):
    """Spend a dense-key budget; rank by the expected worst-slot precision.

    The ring, sparse key and base modulus are taken from `base`. The output
    keeps `application_levels` multiplications that return to the refresh
    input level, and the tree scale takes what the other choices leave.
    """
    candidates = []
    for gap in gap_bits:
        input_scale = base.base_bits - gap
        if input_scale <= 0 or base.message_bound >= 2.0 ** gap / 4:
            continue
        for output_scale in range(max(12, input_scale - 6), input_scale + 13):
            drop = 2 * output_scale - input_scale
            if drop <= 0:
                continue
            for s2c in s2c_bits:
                output = base.base_bits + s2c + application_levels * drop
                for plain in plain_bits:
                    for gain in indicator_gain_bits:
                        rest = maximum_dense_bits - plain - gain - output + output_scale
                        scale = min(rest // (base.depth + 1), 60)
                        if scale < output_scale:
                            continue
                        p = replace(base, input_scale_bits=input_scale, s2c_bits=s2c,
                                    scale_bits=scale, indicator_scale_bits=scale + gain,
                                    plain_bits=plain, output_bits=output,
                                    output_scale_bits=output_scale)
                        try:
                            result = estimate(p, model)
                        except ValueError:
                            continue
                        candidates.append(
                            (result["expected_worst_slot_precision_bits"], p, result))
    candidates.sort(key=lambda item: (-round(item[0], 2), item[1].selector_bits))
    return [{"parameters": asdict(p), "expected_worst_slot_precision_bits": bits,
             "rms_precision_bits": result["rms_precision_bits"],
             "sine_error_bound": result["sine_error_bound"],
             "selector_bits": p.selector_bits} for bits, p, result in candidates]


# TFHEpp schedules. Noise is the integer RMS of continuous sigma 3.6
# truncated toward zero, which is at least the sigma 3.2 of the paper's caps.
_RMS = truncated_normal_rms(3.6)
PROFILES = {
    # TFHEpp CKKSTuttiLL13Config: search(...) optimum at the weight-512 cap.
    "ll13": TuttiParams(n=8192, blocks=64, stride=4, base_bits=24, input_scale_bits=20,
                        s2c_bits=14, scale_bits=22, indicator_scale_bits=24, plain_bits=18,
                        output_bits=58, output_scale_bits=20, auxiliary_bits=20,
                        primary_bits=12, bbar_bits=8, sparse_digit_bits=1, sigma=_RMS,
                        dense_weight=512),
    # TFHEpp CKKSTuttiBooleanLL13Config. The noise screen of this profile
    # covers the unit-circle stage only; Boolean values sit on the sine peaks.
    "ll13_boolean": TuttiParams(n=8192, blocks=64, stride=4, base_bits=24,
                                input_scale_bits=20, s2c_bits=14, scale_bits=22,
                                indicator_scale_bits=24, plain_bits=18, output_bits=46,
                                output_scale_bits=20, auxiliary_bits=20, primary_bits=12,
                                bbar_bits=8, sparse_digit_bits=1, sigma=_RMS,
                                dense_weight=512, boolean_refresh=True),
}
