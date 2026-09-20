"""Conditional error envelope for TFHEpp's scalar Boolean/BFV round trip.

All errors are normalized by their ciphertext modulus. Primitive variances
below are *assumed sub-Gaussian proxy scales*, not proved tail bounds. The
composition uses absolute envelopes, so shared keys and dependent inputs do
not require independent-error assumptions. DD integer FFT recovery is assumed
exact; non-DD bootstrap FFT error has a separate deterministic allowance.
"""

from dataclasses import asdict, dataclass, replace
import math


@dataclass(frozen=True)
class Parameters:
    products: int = 256
    target_failure_bits: int = 40
    ring_n: int = 4096
    qbits: int = 128
    plaintext_bits: int = 16
    ring_alpha_log2: int = -105
    ring_levels: int = 4
    ring_basebit: int = 21
    io_n: int = 1024
    io_alpha_log2: int = -25
    io_levels: int = 3
    io_basebit: int = 6
    half_n: int = 760
    half_alpha_log2: int = -17
    to_half_levels: int = 10
    to_half_basebit: int = 3
    to_io_levels: int = 7
    to_io_basebit: int = 2
    boolean_fft_error: float = 2.0**-20
    forward_digit_bits: int = 1

    def validate(self):
        if self.forward_digit_bits not in (1, 2) or (self.plaintext_bits // 2) % self.forward_digit_bits:
            raise ValueError("forward digits must contain one or two complete bits")
        if (self.qbits, self.io_n, self.half_n) != (128, 1024, 760) or self.ring_n not in (4096, 8192):
            raise ValueError("this envelope models the fixed scalar source dimensions")
        if not 2 <= self.plaintext_bits <= 64 or self.plaintext_bits % 2:
            raise ValueError("plaintext_bits must be even and in 2..64")
        if self.products < 1 or self.target_failure_bits < 1:
            raise ValueError("products and target_failure_bits must be positive")
        if self.boolean_fft_error < 0 or not math.isfinite(self.boolean_fft_error):
            raise ValueError("FFT error must be finite and nonnegative")
        if (not 1 <= self.ring_levels <= 32 or not 1 <= self.ring_basebit <= 21
                or self.ring_levels * self.ring_basebit > self.qbits
                or self.ring_basebit + 16 + (self.ring_n.bit_length() - 1) + 3 >= 53):
            raise ValueError("invalid ring gadget or unsafe DD FFT digit budget")
        if not -124 <= self.ring_alpha_log2 <= -1:
            raise ValueError("ring noise must have integer-unit sigma >=16")
        for levels, basebit in ((self.to_half_levels, self.to_half_basebit),
                               (self.to_io_levels, self.to_io_basebit)):
            if not 1 <= levels <= 10 or not 1 <= basebit <= 8 or levels * basebit > 32:
                raise ValueError("switch decomposition must fit the 32-bit target")


def parameters_for_width(bits=8, *, legacy=False, **kwargs):
    """Compiled width profiles; legacy reproduces the original failed sweep.

    Security must be screened separately: allowing a value in Parameters is
    not a claim that its LWE instance is secure.
    """
    if not 1 <= bits <= 32:
        raise ValueError("operand width must be in 1..32")
    p = Parameters(plaintext_bits=2 * bits, **kwargs)
    if bits == 16 and not legacy:
        p = replace(p, ring_n=8192, ring_alpha_log2=-124,
                    ring_levels=12, ring_basebit=10,
                    to_io_levels=6, to_io_basebit=3, forward_digit_bits=2)
    return p


def iks_std(n, secret_second, levels, basebit, alpha):
    # Rounding of n source coefficients; up to n*levels selected encrypted
    # digits. No erroneous multiplication of key noise by the digit value.
    variance = n * secret_second * 2.0**(-2 * levels * basebit) / 12
    variance += n * levels * alpha**2
    # Torus-width conversion and integer rounding allowance, in 32-bit units.
    return math.sqrt(variance + (n * secret_second + 1) * 2.0**-64 / 12)


def pbs_std(n, levels, basebit, alpha, half_n, qbits):
    key_noise = 2 * levels * n * (2.0**(2 * basebit) + 2) / 12 * alpha**2
    decomposition = (n * (2 / 3) + 1) * 2.0**(-2 * levels * basebit) / 12
    # Binary selector: E[s^2]=1/2. A blind rotation has half_n CMUXes.
    integer_terms = (n / 6 + 1 / 16) * 2.0**(-2 * qbits)
    return math.sqrt(half_n * (key_noise + decomposition / 2 + integer_terms))


def automorphism_std(p):
    # Same proxy for a half-TRGSW switch to s(X^d) or s^2. The n^2 term
    # includes multiplication of decomposition residues by the key polynomial.
    n, levels, b = p.ring_n, p.ring_levels, p.ring_basebit
    key_noise = levels * n * (2.0**(2 * b) + 2) / 12 * 2.0**(2 * p.ring_alpha_log2)
    truncation = (n * 2 / 3)**2 * 2.0**(-2 * levels * b) / 12
    integer_terms = n * n / 9 * 2.0**(-2 * p.qbits)
    return math.sqrt(key_noise + truncation + integer_terms)


def primitive_scales(p):
    return {
        "input": 2.0**p.io_alpha_log2,
        "ks1h": iks_std(p.io_n, 2 / 3, p.to_half_levels, p.to_half_basebit,
                        2.0**p.half_alpha_log2),
        "ks31": iks_std(p.ring_n, 2 / 3, p.to_io_levels, p.to_io_basebit,
                        2.0**p.io_alpha_log2),
        "br3": math.sqrt((p.half_n / 4 + 1) / 12) / (2 * p.ring_n),
        "br1": math.sqrt((p.half_n / 4 + 1) / 12) / (2 * p.io_n),
        "pbs3": pbs_std(p.ring_n, p.ring_levels, p.ring_basebit,
                        2.0**p.ring_alpha_log2, p.half_n, p.qbits),
        "pbs1": pbs_std(p.io_n, p.io_levels, p.io_basebit,
                        2.0**p.io_alpha_log2, p.half_n, 32),
        "auto3": automorphism_std(p),
        "relin3": automorphism_std(p),
    }


def event_counts(p):
    # All output coefficients are covered during each packing automorphism
    # and relinearization. Only extracted PBS coefficients enter the circuit.
    bits = p.plaintext_bits // 2
    return {name: count * p.products for name, count in {
        "input": 2 * bits, "ks1h": (5 if p.forward_digit_bits == 2 else 4) * bits + 1, "ks31": bits + 1,
        "br3": 2 * bits // p.forward_digit_bits,
        "pbs3": 2 * bits // p.forward_digit_bits,
        "br1": (5 if p.forward_digit_bits == 2 else 3) * bits,
        "pbs1": (5 if p.forward_digit_bits == 2 else 3) * bits,
        "auto3": 2 * (p.ring_n.bit_length() - 1) * p.ring_n, "relin3": p.ring_n,
    }.items()}


def product_envelope(input_error, relin_error, p):
    """Worst-case coefficient envelope for the exact power-of-two BFV tensor.

    For integral lifts v_i = Delta*m_i + e_i + Q*k_i, retain t*k_i*e_j;
    dropping these terms would severely understate BFV scale-and-round noise.
    The digit representation has absolute coefficient < Q, hence |k_i| <
    n+2 for constant messages <sqrt(t) and errors < Delta/2. Negacyclic
    convolution obeys ||a*b||_inf <= n||a||_inf||b||_inf.
    """
    n, t, q = p.ring_n, 2**p.plaintext_bits, 2**p.qbits
    if not 0 <= input_error < 1 / (2 * t):
        return math.inf
    lifts = n + 2
    linear = 2 * (2**(p.plaintext_bits // 2) - 1) * input_error
    cross = 2 * t * n * lifts * input_error
    quadratic = n * t * input_error**2
    # Each tensor component is rounded by <1; the cross and square terms
    # decrypt against s and s^2, with l1 norms <=n and <=n^2.
    rounding = (1 + n + n * n) / q
    return linear + cross + quadratic + rounding + relin_error


def envelope(p, cutoff):
    p.validate()
    if cutoff < 0 or not math.isfinite(cutoff):
        raise ValueError("cutoff must be finite and nonnegative")
    errors = {name: cutoff * sigma for name, sigma in primitive_scales(p).items()}
    errors["pbs1"] += p.boolean_fft_error
    # Trace halves, applies an automorphism and adds. Its norm is <=1;
    # truncating the two ciphertext components costs at most (n+1)/Q.
    bits = p.plaintext_bits // 2
    forward = (bits // p.forward_digit_bits) * errors["pbs3"] + (p.ring_n.bit_length() - 1) * (errors["auto3"] + (p.ring_n + 1) * 2.0**-p.qbits)
    product = product_envelope(forward, errors["relin3"], p)
    checks = {
        "forward_boolean_address": (errors["input"] + errors["ks1h"] + errors["br3"], 1 / 8),
        "bfv_plaintext_rounding": (product, 2.0**(-p.plaintext_bits - 1)),
        "guard_digit_address": (product * 2**p.plaintext_bits + errors["ks31"] + errors["ks1h"] + errors["br1"], 1 / 4),
        "output_boolean_decryption": (errors["pbs1"], 1 / 8),
    }
    if p.forward_digit_bits == 2:
        checks["forward_boolean_address"] = (errors["input"] + errors["ks1h"] + errors["br1"], 1 / 8)
        checks["forward_packed_address"] = (3 * errors["pbs1"] + errors["ks1h"] + errors["br3"], 1 / 16)
    # Exact scalar reference gives margin >=1/16 for every data digit.
    # Its true minimum is slightly larger; retain the conservative uniform
    # margin. Guard PBS output replaces, rather than propagates, old noise.
    for digit in range(1, bits + 1):
        digit_error = product * 2**(p.plaintext_bits - 2 * digit) + errors["ks31"] + errors["pbs1"]
        small_error = digit_error + errors["ks1h"]
        if digit < bits:
            checks[f"digit_{digit}_correction"] = (small_error + errors["br1"], 1 / 16)
        checks[f"digit_{digit}_low_bit"] = (2 * small_error + errors["br1"], 1 / 8)
        checks[f"digit_{digit}_high_bit"] = (small_error + errors["br1"], 1 / 16)
    return {
        "errors": errors, "forward_error": forward, "product_error": product,
        "checks": {name: {"error": error, "radius": radius, "passes": error < radius}
                   for name, (error, radius) in checks.items()},
        "passes": all(error < radius for error, radius in checks.values()),
    }


def key_bytes(p):
    # FFT double rows, k=1; identity switches store balanced positive digits.
    bfv_bk = p.half_n * 2 * p.ring_levels * 8 * 2 * p.ring_n * 8
    packing = (p.ring_n.bit_length() - 1) * p.ring_levels * 8 * 2 * p.ring_n * 8
    relin = p.ring_levels * 8 * 2 * p.ring_n * 8
    boolean_bk = p.half_n * 2 * p.io_levels * 2 * p.io_n * 8
    ks1h = p.io_n * p.to_half_levels * 2**(p.to_half_basebit - 1) * (p.half_n + 1) * 4
    ks31 = p.ring_n * p.to_io_levels * 2**(p.to_io_basebit - 1) * (p.io_n + 1) * 4
    return bfv_bk + packing + relin + boolean_bk + ks1h + ks31


def analyse(p=Parameters()):
    p.validate()
    lo, hi = 0.0, 64.0
    for _ in range(64):
        mid = (lo + hi) / 2
        if envelope(p, mid)["passes"]:
            lo = mid
        else:
            hi = mid
    cutoff = lo * 0.99
    result = envelope(p, cutoff)
    count = sum(event_counts(p).values())
    # Assumption: each primitive centered error has tail <=2 exp(-z^2/2)
    # at the given proxy scale, uniformly over the inputs used in this run.
    # Union bound remains valid with correlated/reused evaluation keys.
    log_failure = min(0.0, math.log2(2 * count) - cutoff**2 / (2 * math.log(2)))
    return {
        "parameters": asdict(p), "primitive_assumed_subgaussian_scales": primitive_scales(p),
        "events": event_counts(p), "total_events": count,
        "maximum_uniform_cutoff": lo, "cutoff_with_slack": cutoff,
        "conditional_log2_failure": log_failure,
        "meets_target_conditionally": result["passes"] and log_failure <= -p.target_failure_bits,
        "evaluation_key_bytes": key_bytes(p),
        "assumptions": [
            "Each modeled primitive error obeys the stated sub-Gaussian tail uniformly over this circuit's inputs; variance estimates alone do not prove this.",
            "All DD FFT integer digit convolutions, including key generation and BFV tensor multiplication, round exactly; this is not a proved floating-point bound.",
            "Non-DD Boolean PBS numerical error is bounded by boolean_fft_error in normalized torus units.",
            "Primitive marginal bounds apply with reused keys and dependent ciphertexts; no independence between circuit errors is assumed in composition.",
            "Uniform ternary ring/IO secrets and uniform binary half-level secrets, using the source parameter family.",
        ],
        "envelope": result,
    }


def search(p=Parameters(), cap_bytes=24 * 1024**3):
    """Keep the default if it passes; otherwise minimize stored switching keys.

    Four GiB is reserved for buffers, OpenMP stacks, and allocator overhead.
    Both switches may change their gadget, never their security instance.
    """
    baseline = analyse(p)
    if baseline["meets_target_conditionally"] and key_bytes(p) + 4 * 1024**3 <= cap_bytes:
        return baseline
    choices = [(l, b) for l in range(1, 11) for b in range(1, 9) if l * b <= 32]
    candidates = []
    for lh, bh in choices:
        for li, bi in choices:
            candidate = replace(p, to_half_levels=lh, to_half_basebit=bh,
                                to_io_levels=li, to_io_basebit=bi)
            size = key_bytes(candidate)
            if size + 4 * 1024**3 <= cap_bytes:
                candidates.append((size, lh + li, lh, bh, li, bi, candidate))
    for *_, candidate in sorted(candidates, key=lambda item: item[:-1]):
        result = analyse(candidate)
        if result["meets_target_conditionally"]:
            return result
    return dict(baseline, search_failure="No candidate meets both the conditional bound and memory cap")


def tune_16bit(products=256, target_failure_bits=40, cap_bytes=24 * 1024**3):
    """Bounded search in the separately screened N=8192, sigma=16 instance.

    Retain eight extra failure-exponent bits beyond the requested target.
    Minimize primary gadget rows first, then stored key bytes; this is a
    transparent cost proxy, not a claim of globally optimal measured runtime.
    Reserve four GiB for buffers and stacks. Full 8x16-bit DD is fixed.
    """
    # Reproduce the original parameter selection independently of the later
    # paired-forward optimization. New candidates require new measurements.
    base = replace(parameters_for_width(16, products=products,
                   target_failure_bits=target_failure_bits), forward_digit_bits=1)
    candidates = []
    tested = 0
    for levels in range(4, 17):
        for basebit in range(2, 20):
            if levels * basebit > 128 or basebit + 16 + 13 + 3 >= 53:
                continue
            for switch_levels, switch_basebit in ((8, 2), (6, 3), (7, 3), (5, 4), (4, 5)):
                p = replace(base, ring_levels=levels, ring_basebit=basebit,
                            to_io_levels=switch_levels, to_io_basebit=switch_basebit)
                size = key_bytes(p)
                if size + 4 * 1024**3 > cap_bytes:
                    continue
                tested += 1
                result = analyse(p)
                if result["meets_target_conditionally"] and result["conditional_log2_failure"] <= -target_failure_bits - 8:
                    candidates.append((levels, size, result))
    if not candidates:
        raise ValueError("no 16-bit candidate meets the bound, safety margin, and memory cap")
    result = min(candidates, key=lambda item: item[:2])[2]
    return dict(result, tuning={"candidates_tested": tested,
                                "passing_candidates": len(candidates),
                                "extra_failure_bits": 8,
                                "cap_bytes": cap_bytes,
                                "buffer_reserve_bytes": 4 * 1024**3,
                                "ranking": "primary gadget rows, then stored key bytes"})


def reverse_reference(message, error=0, qbits=128, injections=None, plaintext_bits=16):
    """Exact integer trace of the implemented guard/HomDecomp/bit LUTs.

    error and injected phase offsets are in 2^qbits torus units. Keys and
    ciphertexts are absent: this tests the semantic circuit, not cryptography.
    injections maps (stage, digit) to signed integer offsets.
    """
    if (not 2 <= plaintext_bits <= 64 or plaintext_bits % 2
            or not 0 <= message < 2**plaintext_bits or qbits < plaintext_bits + 4):
        raise ValueError("invalid message, plaintext width, or torus width")
    inj = injections or {}
    q = 1 << qbits
    phase = (message * (q >> plaintext_bits) + (q >> (plaintext_bits + 3)) + error) % q
    previous = 0
    output = 0
    for digit in range(plaintext_bits // 2 + 1):
        c = ((phase << (plaintext_bits - 2 * digit)) + inj.get(("ks31", digit), 0)) % q
        if digit:
            c = (c + previous - q // 8) % q
        p = (c + q // 8 + inj.get(("ks1h", digit), 0)) % q
        address = (p + inj.get(("br1", digit), 0)) % q
        if digit:
            # Extraction uses its own identity switch and each PBS its own
            # address rounding. It does not reuse the correction's rounding.
            bits_phase = (c + q // 8 + inj.get(("extract_ks", digit), 0)) % q
            low = (2 * bits_phase + inj.get(("low_br", digit), 0)) % q
            high = (bits_phase + inj.get(("high_br", digit), 0)) % q
            output |= int(low >= q // 2) << (2 * (digit - 1))
            output |= int(high >= q // 2) << (2 * (digit - 1) + 1)
        previous = (q // 16 if address < q // 2 else -q // 16)
        previous += inj.get(("pbs1", digit), 0)
    return output
