"""Public frequency plans and level/amplitude budgets for ePrint 2026/2119.

These offline checks do not certify ciphertext noise, FFT accuracy or security.
"""

from dataclasses import dataclass
import math

from .narrow_ckks import retained_square_range


@dataclass(frozen=True)
class FrequencyBasis:
    degree: int
    basis: tuple[int, ...]
    # Entry k-1 is (k,) for a directly generated root or (a,b) for a product.
    reconstruction: tuple[tuple[int, ...], ...]
    search_nodes: int


def minimum_positive_frequency_basis(degree: int) -> FrequencyBasis:
    """Find a smallest A with {1,...,degree} contained in A union (A+A).

    Equal summands are allowed. This deliberately excludes signed summands:
    negative frequencies are obtained after weighting by conjugation, as in
    the paper. The degree limit bounds this offline exhaustive search.
    """
    if isinstance(degree, bool) or not isinstance(degree, int) or not 1 <= degree <= 40:
        raise ValueError("basis planning requires an integer degree from 1 to 40")
    target = (1 << (degree + 1)) - 2
    nodes = 0

    def search(basis, covered, size):
        nonlocal nodes
        nodes += 1
        missing = target & ~covered
        if not missing:
            return basis
        remaining = size - len(basis)
        # Each new member contributes one singleton, pairs with existing
        # members, and pairs (including squares) among the new members.
        if not remaining or missing.bit_count() > (
            remaining * (len(basis) + 1) + remaining * (remaining + 1) // 2
        ):
            return None
        gap = (missing & -missing).bit_length() - 1
        # A next member larger than the first gap cannot fill that gap:
        # all later singletons and all new positive sums would be larger.
        for value in range(basis[-1] + 1 if basis else 1, gap + 1):
            added = 1 << value
            if 2 * value <= degree:
                added |= 1 << (2 * value)
            for previous in basis:
                if previous + value <= degree:
                    added |= 1 << (previous + value)
            found = search(basis + (value,), covered | added, size)
            if found is not None:
                return found
        return None

    # At most size singletons and size*(size+1)/2 sums can be distinct.
    for size in range(1, degree + 1):
        if size * (size + 3) // 2 < degree:
            continue
        basis = search((), 0, size)
        if basis is not None:
            break
    reconstruction = []
    for frequency in range(1, degree + 1):
        if frequency in basis:
            reconstruction.append((frequency,))
        else:
            reconstruction.append(next(
                (a, b) for i, a in enumerate(basis) for b in basis[i:]
                if a + b == frequency
            ))
    return FrequencyBasis(degree, basis, tuple(reconstruction), nodes)


def frequency_basis_schedule_budget(params, degree, coefficient_bits=32,
                                    component_bound=4.01, auxiliary_bits=0):
    """Exact level/scale and conservative amplitude budgets, excluding noise.

    Coefficient quantization bounds the changed public trigonometric series
    at exact unit roots. It does not bound ciphertext reconstruction errors.
    Normalization checks use the existing rescale-only refresh input adapter;
    they do not assume an ability to raise the scale or restore lost modulus.
    """
    params.validate()
    plan = minimum_positive_frequency_basis(degree)
    if (isinstance(coefficient_bits, bool) or not isinstance(coefficient_bits, int)
            or not 1 <= coefficient_bits <= 40):
        raise ValueError("coefficient bits must be an integer from 1 to 40")
    if isinstance(auxiliary_bits, bool) or not isinstance(auxiliary_bits, int) or auxiliary_bits < 0:
        raise ValueError("auxiliary bits must be a nonnegative integer")
    root_bits = params.output_bits
    output_bits = root_bits - params.scale_bits
    output_scale = params.scale_bits + coefficient_bits
    if root_bits <= 2 * params.scale_bits or output_bits <= output_scale:
        raise ValueError("basis schedule lacks unscaled DD product or coefficient-scale headroom")
    square = retained_square_range(output_bits, output_scale, component_bound)
    input_bits = params.source_bits + params.s2c_bits
    normalization_drop = output_scale - params.input_scale_bits
    max_evaluation_key = max(params.selector_bits,
        params.selector_bits - params.plain_bits + auxiliary_bits,
        input_bits + auxiliary_bits, root_bits, output_bits)
    max_key = max(max_evaluation_key, params.packed_bits)
    return {
        "scope": "Basis levels, fixed-point coefficient quantization and slot-amplitude capacity; excludes noise/FFT/security",
        "degree": degree, "basis": list(plan.basis), "positive_blind_rotations": len(plan.basis),
        "reconstructed_products": sum(len(term) == 2 for term in plan.reconstruction),
        "extra_product_input_bits": root_bits, "output_bits": output_bits,
        "coefficient_bits": coefficient_bits, "output_scale_bits": output_scale,
        "coefficient_quantization_uniform_bound_at_exact_roots":
            math.ldexp(degree * math.sqrt(2) + 0.5, -coefficient_bits),
        "retained_square": square,
        "minimum_selector_bits_for_guarded_square":
            params.selector_bits + max(0, square["minimum_output_bits"] - output_bits),
        "refresh_input_bits": input_bits, "normalization_scale_drop": normalization_drop,
        "output_can_normalize_for_refresh":
            normalization_drop >= 0 and output_bits >= input_bits + normalization_drop,
        "retained_square_can_normalize_for_refresh":
            normalization_drop >= 0 and output_bits - output_scale >= input_bits + normalization_drop,
        "retained_square_worst_case_component_bound": 2 * component_bound**2,
        "refresh_phase_component_limit": 2**(params.source_bits - params.input_scale_bits - 1),
        "retained_square_fits_refresh_phase_domain":
            2 * component_bound**2 < 2**(params.source_bits - params.input_scale_bits - 1),
        "retained_square_fits_declared_message_domain": 2 * component_bound**2 <= params.message_bound,
        "maximum_evaluation_key_bits": max_evaluation_key,
        "packed_key_bits": params.packed_bits, "maximum_key_bits": max_key,
        "minimum_native_storage_bits_including_optional_packed_key": 64 * ((max_key + 63) // 64),
        "noise_validated": False, "security_validated": False,
    }
