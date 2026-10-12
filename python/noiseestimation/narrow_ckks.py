"""Average-case factor/tree screen for ePrint 2026/2119's TFHEpp DD port.

This supports TFHEpp DD multiplication with rescaling before or after
relinearization, without a hybrid P modulus. Bbar limbs are exact arithmetic and introduce no
independent encryption error. It reuses the GL estimator's variance arithmetic.
Optional auxiliary switching lifts its modulus and scale before DD switching
and divides afterward; S2C uses a whole-pipeline lift screen.
S2C/encapsulation use an average-case BSGS variance screen. FFT error and
structured-secret security are excluded. This is not a correctness certificate.
"""

from dataclasses import dataclass, replace
import math

from .gl import ProductState, combine_variances, log2_add, NEG_INF


@dataclass(frozen=True)
class NarrowParams:
    n: int = 8192
    h: int = 31
    selector_bits: int = 184
    scale_bits: int = 20
    plain_bits: int = 20
    source_bits: int = 25
    input_scale_bits: int = 21
    s2c_bits: int = 20
    primary_bits: int = 4
    bbar_bits: int = 4
    sigma: float = 3.2
    message_bound: float = 0.125

    @property
    def depth(self):
        return self.h.bit_length()

    @property
    def output_bits(self):
        return self.selector_bits - self.plain_bits - self.depth * self.scale_bits

    @property
    def packed_bits(self):
        return self.selector_bits + self.plain_bits + self.n.bit_length() - 1

    def validate(self):
        if self.n < 16 or self.n & (self.n - 1) or not 1 < self.h < self.n:
            raise ValueError("require power-of-two N >= 16 and 1 < h < N")
        width = (self.n + self.h - 1) // self.h
        if (self.h - 1) * width >= self.n or width - 1 > self.n // 4:
            raise ValueError("invalid final block or tail packing geometry")
        if min(self.scale_bits, self.plain_bits, self.primary_bits, self.bbar_bits) <= 0:
            raise ValueError("scale and decomposition bit counts must be positive")
        if not 0 < self.input_scale_bits < self.source_bits or self.s2c_bits <= 0:
            raise ValueError("invalid input scale/S2C schedule")
        if not math.isfinite(self.sigma) or self.sigma < 0:
            raise ValueError("sigma must be finite and nonnegative")
        if not math.isfinite(self.message_bound) or self.message_bound <= 0:
            raise ValueError("message bound must be positive and finite")
        if self.message_bound >= 2**(self.source_bits - self.input_scale_bits - 1):
            raise ValueError("message domain must lie inside half the phase gap")
        if self.output_bits <= max(self.scale_bits, self.source_bits + self.s2c_bits):
            raise ValueError("output modulus does not permit a repeatable refresh")


def minimum_square_output_bits(scale_bits, component_bound, guard_bits=3):
    """Amplitude budget for arbitrary complex slots after a fixed-scale DD square.

    For real/imaginary components bounded by R, |z^2| <= 2R^2. Canonical
    inverse encoding bounds every polynomial coefficient by that magnitude.
    After a square Q loses scale_bits while scale remains fixed, so centered
    coefficient capacity is 2^(Qout-2*scale_bits-1). Require strict containment
    and add guard bits. Include approximation/evaluation error in R. This is
    an amplitude budget, not a noise or security estimate.
    """
    if not isinstance(scale_bits, int) or scale_bits <= 0:
        raise ValueError("scale bits must be a positive integer")
    if not isinstance(guard_bits, int) or guard_bits < 0:
        raise ValueError("guard bits must be a nonnegative integer")
    if not math.isfinite(component_bound) or component_bound <= 0:
        raise ValueError("component bound must be positive and finite")
    return max(scale_bits + 1,
               math.floor(2 * scale_bits + 2 + 2 * math.log2(component_bound)) + 1) + guard_bits


def retained_square_range(output_bits, scale_bits, component_bound, guard_bits=3):
    minimum = minimum_square_output_bits(scale_bits, component_bound, guard_bits)
    if not isinstance(output_bits, int) or output_bits <= scale_bits:
        raise ValueError("output modulus must permit a DD square")
    required_log2 = 1 + 2 * math.log2(component_bound)
    capacity_log2 = output_bits - 2 * scale_bits - 1
    return {"scope": "Worst-case complex slot amplitude after one fixed-scale DD square; excludes noise/security",
            "component_bound": component_bound, "guard_bits": guard_bits,
            "square_log_q": output_bits - scale_bits, "square_log_delta": scale_bits,
            "square_coefficient_bound_log2": required_log2,
            "centered_coefficient_capacity_log2": capacity_log2,
            "amplitude_fits": capacity_log2 > required_log2,
            "minimum_output_bits": minimum,
            "guarded_range_passed": output_bits >= minimum}


def switch_slot_variance(p, modulus_bits, digit_bits, auxiliary_bits=0, switch_count=1):
    """Raw gadget switch, optionally lifted then divided to its data level.

    Count all raised-level rows conservatively, including the zero lower rows
    of a power-of-two lift. A binary rotation chain has switch_count key
    switches but only one final division of its two ciphertext components.
    Key errors within the chain use the independent average-case screen.
    """
    if not isinstance(switch_count, int) or switch_count < 0:
        raise ValueError("switch count must be a nonnegative integer")
    if switch_count == 0:
        return NEG_INF
    if p.sigma == 0:
        noise = NEG_INF
    else:
        rows = (modulus_bits + auxiliary_bits + digit_bits - 1) // digit_bits
        noise = (math.log2(switch_count * rows) + 2 * math.log2(p.n) +
                 log2_add(2 * digit_bits, 1) - math.log2(12) +
                 2 * math.log2(p.sigma) - 2 * (p.scale_bits + auxiliary_bits))
    if not auxiliary_bits:
        return noise
    rounding = math.log2(p.n * (1 + p.n * 2 / 3) / 12) - 2 * p.scale_bits
    return log2_add(noise, rounding)


def retained_unit_product_chain(p, levels, product_scale_bits, component_bound,
                                auxiliary_bits=0, model="independent", guard_bits=3):
    """Screen HT's retained levels with fresh encrypted unit operands.

    Each multiplication keeps the signal amplitude, so tiny powers cannot
    make this capacity/noise screen pass trivially. DD switches before
    rescaling. A preliminary scale reduction preserves the input error and
    adds two-component rounding at the new scale. FFT/security are excluded.
    """
    p.validate()
    if not isinstance(levels, int) or levels <= 0:
        raise ValueError("retained levels must be a positive integer")
    if not isinstance(product_scale_bits, int) or not 0 < product_scale_bits <= p.scale_bits:
        raise ValueError("product scale must be positive and no greater than the refresh scale")
    if not math.isfinite(component_bound) or component_bound <= 0:
        raise ValueError("component bound must be positive and finite")
    if not isinstance(guard_bits, int) or guard_bits < 0:
        raise ValueError("guard bits must be a nonnegative integer")
    drop = p.scale_bits - product_scale_bits
    initial_q = p.output_bits - drop
    final_q = initial_q - levels * product_scale_bits
    if final_q <= product_scale_bits:
        raise ValueError("output modulus cannot retain the requested product levels")
    magnitude_log2 = math.log2(component_bound) + 0.5
    capacity_log2 = final_q - product_scale_bits - 1
    required_output = max(drop + (levels + 1) * product_scale_bits + 1,
        drop + math.floor((levels + 1) * product_scale_bits + 1 + magnitude_log2) + 1 + guard_bits)
    base = estimate(p, model, True, True, auxiliary_bits)
    variance = 2 * math.log2(base["complex_refresh_rms"])
    rounding = math.log2(p.n * (1 + p.n * 2 / 3) / 12) - 2 * product_scale_bits
    if drop:
        variance = combine_variances((variance, rounding), model)
    fresh_unit = NEG_INF if p.sigma == 0 else math.log2(p.n * p.sigma**2) - 2 * product_scale_bits
    steps = []
    for level in range(levels):
        q = initial_q - level * product_scale_bits
        switch = switch_slot_variance(replace(p, scale_bits=2 * product_scale_bits),
                                      q, p.primary_bits)
        variance = combine_variances((variance,
            fresh_unit + 2 * magnitude_log2, variance + fresh_unit,
            rounding, switch), model)
        steps.append({"input_bits": q, "output_bits": q - product_scale_bits,
                      "rms": 2**(variance / 2)})
    return {"scope": "Average-case late-DD retained chain with encrypted unit operands; FFT/security excluded",
            "retained_levels": levels, "scale_reduction_bits": drop,
            "product_scale_bits": product_scale_bits, "final_output_bits": final_q,
            "component_bound": component_bound, "guard_bits": guard_bits,
            "centered_coefficient_capacity_log2": capacity_log2,
            "coefficient_bound_log2": magnitude_log2,
            "guarded_range_passed": capacity_log2 > magnitude_log2 + guard_bits,
            "minimum_refresh_output_bits": required_output,
            "steps": steps, "eight_rms_plus_sine_screen":
                8 * 2**(variance / 2) + base["complex_sine_error_bound"],
            "security_validated": False}


def s2c_variance(p, model, auxiliary_bits=0):
    """Natural-order one-level BSGS coefficient variance.

    An auxiliary-scale variant lifts the entire S2C/encapsulation pipeline,
    then divides once under the structured output key. Input encryption error
    is preserved. The optional TFHEpp refresh schedule uses this ordering;
    this average-case screen still excludes FFT error and security analysis.
    """
    if auxiliary_bits:
        lifted = replace(p, source_bits=p.source_bits + auxiliary_bits,
                         input_scale_bits=p.input_scale_bits + auxiliary_bits)
        _, parts = s2c_variance(lifted, model)
        parts["fresh"] = NEG_INF if p.sigma == 0 else (
            math.log2(p.n * p.sigma**2) - 2 * p.input_scale_bits)
        # These are coefficient errors: no extra slot-embedding factor N.
        parts["auxiliary_final_rounding"] = math.log2((1 + p.h) / 12) - 2 * p.input_scale_bits
        return combine_variances(parts.values(), model), parts
    input_p = replace(p, scale_bits=p.input_scale_bits)
    n = p.n
    slots = n // 2
    baby = 1 << ((n.bit_length() - 2) // 2)
    giants = (slots + baby - 1) // baby
    before = switch_slot_variance(input_p, p.source_bits + p.s2c_bits, p.bbar_bits) - math.log2(n)
    after = switch_slot_variance(input_p, p.source_bits, p.bbar_bits) - math.log2(n)
    baby_switches = sum(i.bit_count() for i in range(baby))
    giant_switches = sum((i * baby).bit_count() for i in range(giants))
    fresh = NEG_INF if p.sigma == 0 else math.log2(n * p.sigma**2) - 2 * p.input_scale_bits
    parts = {"fresh": fresh,
        "input_conjugation": before + math.log2(n / 2),
        "baby_rotations": before + math.log2(max(1, giants * baby_switches)),
        "giant_rotations": after + math.log2(max(1, 2 * giant_switches)),
        "encapsulation": after,
        "rescale": math.log2(2 * giants * (1 + n * 2 / 3) / 12) - 2 * p.input_scale_bits,
        "diagonal_encoding": math.log2(n * p.message_bound**2 / 18) - 2 * p.s2c_bits}
    return combine_variances(parts.values(), model), parts


def estimate(p: NarrowParams, model="independent", streamed=False, relinearize_before_rescale=False, auxiliary_switch_bits=0,
             selector_chunk_size=0, batch_columns=1):
    p.validate()
    if not isinstance(auxiliary_switch_bits, int) or auxiliary_switch_bits < 0:
        raise ValueError("auxiliary switch bits must be a nonnegative integer")
    if not isinstance(selector_chunk_size, int) or selector_chunk_size < 0:
        raise ValueError("selector chunk size must be a nonnegative integer")
    if not isinstance(batch_columns, int) or batch_columns < 1:
        raise ValueError("batch columns must be a positive integer")
    n, width = p.n, (p.n + p.h - 1) // p.h
    factor_q = p.selector_bits - p.plain_bits
    dense_gain = n * 2 / 3
    fresh = NEG_INF if p.sigma == 0 else (
        math.log2(2 * width * n) + 2 * math.log2(p.sigma) - 2 * p.scale_bits)
    encoding = math.log2(4 * n / 12) - 2 * p.plain_bits
    rescale = math.log2(2 * n * (1 + dense_gain) / 12) - 2 * p.scale_bits
    conjugate = switch_slot_variance(p, factor_q, p.bbar_bits, auxiliary_switch_bits)
    tail = width - 1
    pack_blocks = min(p.h, n // (4 * tail))
    rotation_switches = max(((2 * (ell + 1) * tail) % (n // 2)).bit_count()
                           for ell in range(pack_blocks))
    rotation = switch_slot_variance(p, factor_q, p.bbar_bits, auxiliary_switch_bits,
                                   rotation_switches)
    # Screen the noisier tail-aligned factor at every leaf conservatively.
    # Head factors do not perform this additional rotation chain.
    factor_v = combine_variances((fresh, encoding, rescale, conjugate, rotation), model)
    gap = 2.0 ** (p.source_bits - p.input_scale_bits)
    beta = gap / (4 * math.pi)
    beta_v = math.log2(n / 12) - 2 * p.scale_bits
    # E||s^2||^2 <= 2*N^2*(2/3)^2 + N*(2/3) for dense ternary s.
    # Rounding the tensor's s^2 component before relin amplifies this term.
    tensor_v = math.log2(n * (1 + dense_gain + 2 * dense_gain**2 + dense_gain) / 12) - 2 * p.scale_bits
    nodes = [ProductState(factor_v, 1)] * p.h + [ProductState(beta_v, beta)]
    levels = []
    for level in range(p.depth):
        q = factor_q - (level + 1) * p.scale_bits
        if relinearize_before_rescale:
            # Switch the exact modular tensor at doubled fixed-point scale,
            # then round only the two-component ciphertext. Rounding is
            # amplified by s, rather than s^2. Q for the switch is the input
            # level; this requires different DD keys from the default path.
            rounding = math.log2(n * (1 + dense_gain) / 12) - 2 * p.scale_bits
            relin = switch_slot_variance(replace(p, scale_bits=2 * p.scale_bits),
                                        q + p.scale_bits, p.primary_bits)
        else:
            rounding = tensor_v
            relin = switch_slot_variance(p, q, p.primary_bits)
        arithmetic = combine_variances((rounding, relin), model)
        next_nodes = []
        for i in range(0, len(nodes), 2):
            if i + 1 == len(nodes):
                next_nodes.append(nodes[i])
                continue
            a, b = nodes[i:i + 2]
            variance = combine_variances((
                a.log2_variance + 2 * math.log2(b.value_bound),
                b.log2_variance + 2 * math.log2(a.value_bound),
                a.log2_variance + b.log2_variance, arithmetic), model)
            next_nodes.append(ProductState(variance, a.value_bound * b.value_bound))
        levels.append({"output_bits": q, "tensor_rounding_log2_variance": rounding,
                       "relin_log2_variance": None if relin == NEG_INF else relin})
        nodes = next_nodes
    final_v = log2_add(nodes[0].log2_variance + 2,
                       switch_slot_variance(p, p.output_bits, p.bbar_bits, auxiliary_switch_bits))
    # Both real/imaginary channels, and conjugation reuses the same error.
    final_v += 1
    packs = (p.h + n // (4 * tail) - 1) // (n // (4 * tail))
    selectors = width * (p.h + packs)
    # Same conservative convolution bound as TFHEpp's 16-bit NTT cache.
    # First degree65536_primes entry in TFHEpp/include/modular_ntt.hpp.
    primes = 1 if (2**p.plain_bits + 1) * 65535 * n * width < 2301972608560791553 // 2 else 2
    cache_bytes = selectors * 2 * ((p.selector_bits + 15) // 16) * n * primes * 8
    streamed_bytes = cache_bytes * width // selectors
    chunk_count = min(selector_chunk_size, width) if selector_chunk_size else width
    chunk_cache_bytes = streamed_bytes * chunk_count // width
    chunk_sum_bytes = (2 * 2 * ((p.selector_bits + 15) // 16) * n * primes * 8 * batch_columns
                       if selector_chunk_size else 0)
    s2c_v, s2c_parts = s2c_variance(p, model, auxiliary_switch_bits)
    full_v = combine_variances((final_v, s2c_v + 1), model)
    sine_bound = (2 * math.pi)**2 * p.message_bound**3 / (6 * gap**2)
    key_bits = max(p.packed_bits, p.selector_bits, factor_q + auxiliary_switch_bits,
                   p.output_bits + auxiliary_switch_bits,
                   p.source_bits + p.s2c_bits + auxiliary_switch_bits)
    storage_bits = 64 * ((key_bits + 63) // 64)
    selector_bytes = 2 * n * (storage_bits // 8)
    first_key_level = 0 if relinearize_before_rescale else 1
    dd_rows = sum(((factor_q - (level + first_key_level) * p.scale_bits + p.primary_bits - 1) // p.primary_bits) *
                  ((factor_q - (level + first_key_level) * p.scale_bits + p.bbar_bits - 1) // p.bbar_bits)
                  for level in range(p.depth))
    compact_digit_bytes = 1 if p.bbar_bits <= 8 else 2 if p.bbar_bits <= 16 else 4 if p.bbar_bits <= 32 else None
    factor_fft_bytes = 2 * n * 8 * ((factor_q + auxiliary_switch_bits + p.bbar_bits - 1) // p.bbar_bits)**2
    output_fft_bytes = 2 * n * 8 * ((p.output_bits + auxiliary_switch_bits + p.bbar_bits - 1) // p.bbar_bits)**2
    rotation_mask = 0
    for ell in range(min(p.h, (n // 2) // (2 * (width - 1)))):
        rotation_mask |= (2 * (ell + 1) * (width - 1)) % (n // 2)
    rotation_keys = rotation_mask.bit_count()
    return {"scope": __doc__, "model": model, "security_validated": False,
            "multiplication_order": "relinearize-then-rescale" if relinearize_before_rescale else "rescale-then-relinearize",
            "output_bits": p.output_bits, "packed_bits": p.packed_bits,
            "auxiliary_switch_bits": auxiliary_switch_bits,
            "structured_switch_key_bits": p.source_bits + auxiliary_switch_bits,
            "maximum_evaluation_key_bits": key_bits,
            "selector_memory": {
                "scope": "Selector torus/NTT storage only; excludes evaluation keys, their FFT caches, batches, and allocator overhead",
                "inferred_torus_storage_bits": storage_bits,
                "expanded_torus_bytes": selectors * selector_bytes,
                "one_group_torus_bytes": width * selector_bytes,
                "file_backed_torus_and_ntt_peak_bytes": width * selector_bytes + streamed_bytes,
                "in_memory_streamed_torus_and_ntt_peak_bytes": selectors * selector_bytes + streamed_bytes},
            "chunked_selector_memory": {
                "scope": "Loaded selector group, bounded NTT chunk and modular batch sums only; excludes evaluation keys, FFT caches and other workspaces",
                "selector_chunk_size": chunk_count,
                "batch_columns": batch_columns,
                "ntt_chunk_bytes": chunk_cache_bytes,
                "modular_sum_bytes": chunk_sum_bytes,
                "file_backed_group_chunk_and_sums_bytes": width * selector_bytes + chunk_cache_bytes + chunk_sum_bytes},
            "raw_dd_chain_torus_bytes": dd_rows * selector_bytes,
            "compact_dd_chain_digit_bytes": None if compact_digit_bytes is None else
                dd_rows * 2 * n * compact_digit_bytes,
            "dense_s2c_plaintext_polynomial_bytes": n * n * (storage_bits // 8),
            "compact_s2c_plaintext_coefficient_bytes": n * n * 8,
            "prepared_fft_memory": {
                "scope": "Derived key payloads only; excludes raw keys, selectors, transform plans and workspaces",
                "factor_automorphism_bytes": factor_fft_bytes,
                "output_automorphism_bytes": output_fft_bytes,
                "dd_chain_bytes": dd_rows * 2 * n * 8,
                "tail_rotation_binary_keys": rotation_keys,
                "tail_rotation_bytes": rotation_keys * factor_fft_bytes,
                "bounded_tail_rotation_peak_bytes": factor_fft_bytes if rotation_keys else 0,
                "staged_factor_phase_bytes": factor_fft_bytes * (2 if rotation_keys else 1),
                "staged_product_tree_phase_bytes": dd_rows * 2 * n * 8,
                "staged_output_phase_bytes": output_fft_bytes,
                "staged_peak_bytes": max(factor_fft_bytes * (2 if rotation_keys else 1),
                                         dd_rows * 2 * n * 8, output_fft_bytes)},
            "factor_log2_variance": factor_v, "levels": levels,
            "maximum_tail_rotation_switches": rotation_switches,
            "tail_rotation_log2_variance": None if rotation == NEG_INF else rotation,
            "complex_factor_tree_rms": 2**(final_v / 2),
            "factor_tree_rms_precision_bits": -final_v / 2,
            "ntt_cache_bytes_upper_bound": cache_bytes,
            "ntt_cache_peak_bytes": streamed_bytes if streamed else cache_bytes,
            "s2c_coefficient_rms": 2**(s2c_v / 2), "s2c_log2_variances": {
                k: None if v == NEG_INF else v for k, v in s2c_parts.items()},
            "complex_refresh_rms": 2**(full_v / 2), "real_sine_error_bound": sine_bound,
            "complex_sine_error_bound": math.sqrt(2) * sine_bound,
            "eight_rms_plus_sine_screen": 8 * 2**(full_v / 2) + math.sqrt(2) * sine_bound,
            "screen_note": "Eight RMS is an engineering screen, not a tail/failure proof",
            "excluded": ["FFT error", "security"]}
