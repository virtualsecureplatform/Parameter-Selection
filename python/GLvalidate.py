#!/usr/bin/env python3
"""Deterministic validation checks for the GL-SHIP noise estimator."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path
import sys


_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent
sys.path[:0] = [str(_REPO_ROOT), str(_SCRIPT_DIR)]

try:
    from python.noiseestimation.gl import (  # noqa: E402
        GLNoiseParams,
        base_ring_convolution_gain,
        dd_key_switch_noise,
        encoded_real_coefficient_variance,
        estimate_gl_ship,
        log2_add,
    )
    from python.noiseestimation.params.gl import PRESETS  # noqa: E402
except ModuleNotFoundError as exc:
    if exc.name != "python":
        raise
    from noiseestimation.gl import (  # type: ignore[no-redef] # noqa: E402
        GLNoiseParams,
        base_ring_convolution_gain,
        dd_key_switch_noise,
        encoded_real_coefficient_variance,
        estimate_gl_ship,
        log2_add,
    )
    from noiseestimation.params.gl import PRESETS  # type: ignore[no-redef] # noqa: E402


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def brute_w_convolution_gain(p: int) -> int:
    phi = p - 1
    counts = [0] * phi
    for left in range(phi):
        for right in range(phi):
            degree = left + right
            if degree < phi:
                counts[degree] += 1
            elif degree == phi:
                for output in range(phi):
                    counts[output] += 1
            else:
                counts[degree - p] += 1
    return max(counts)


def brute_w_inverse_row_norm(p: int, w: int) -> float:
    root = complex(math.cos(2.0 * math.pi / p), math.sin(2.0 * math.pi / p))
    return sum(
        abs(root ** ((-a * w) % p) - root**a) ** 2
        for a in range(1, p)
    ) / (p * p)


def tfhepp_toy_profile() -> GLNoiseParams:
    """Schedule from TFHEpp test/gl/gl_bootstrap.cpp."""
    return GLNoiseParams(
        tag="tfhepp-toy",
        n=2,
        p=5,
        log_q=90,
        log_p=20,
        q0_bits=24,
        gap_bits=4,
        stc_bits=16,
        x_scale_bits=8,
        w_scale_bits=8,
        tree_scale_bits=18,
        outside_multiplicative_depth=0,
        outside_scale_bits=1,
        dnum=1,
        theta=1,
        window_width=3,
        masked_column_count=14,
        primary_bit=8,
        bbar_bit=8,
        storage_bits=128,
        security_limit_log_pq=110,
        paper_precision_bits=0.0,
        sparse_hamming_weight=3,
        hmux_radix=2,
        w_baby_step=2,
        message_bound=0.11,
    )


def main() -> int:
    print("GL-SHIP noise estimator validation")
    check(abs(log2_add(0.0, 0.0) - 1.0) < 1e-12, "log2_add failed")

    expected_margins = {
        "n256p17": (42, 0),
        "n512p17": (18, 0),
        "n1024p17": (35, 7),
    }
    estimates = {}
    for name, params in PRESETS.items():
        params.validate()
        storage_margin, security_margin = expected_margins[name]
        check(params.storage_margin_bits == storage_margin, f"{name} storage margin")
        check(params.security_margin_bits == security_margin, f"{name} security margin")
        check(params.factor_count == 32 and params.tree_depth == 5, f"{name} tree")
        check(params.output_depth_margin_bits >= 0, f"{name} outside depth")
        max_tree_scale = replace(
            params,
            tree_scale_bits=(
                params.tree_scale_bits + params.tree_scale_headroom_bits
            ),
        )
        max_tree_scale.validate()
        try:
            replace(
                max_tree_scale,
                tree_scale_bits=max_tree_scale.tree_scale_bits + 1,
            ).validate()
        except ValueError:
            pass
        else:
            raise AssertionError(f"{name} tree-scale headroom")
        check(params.hmux_stages == 3, f"{name} HMux stages")
        logical_hmux = (
            params.sparse_hamming_weight
            * params.hmux_stages
            * params.hmux_radix
        )
        check(logical_hmux == 372, f"{name} HMux branch count")
        stc_small_keys = (params.w_baby_step - 1) + (params.w_giant_steps - 1)
        check(stc_small_keys == 6, f"{name} StC small-key count")
        check(params.bbar_cover_bits <= params.storage_bits, f"{name} Bbar cover")
        check(
            params.primary_cover_bits <= params.storage_bits,
            f"{name} primary cover",
        )
        check(
            base_ring_convolution_gain(params)
            == 2 * params.n * brute_w_convolution_gain(params.p),
            f"{name} GL convolution gain",
        )
        for w in (0, params.phi // 2, params.phi - 1):
            check(
                abs(brute_w_inverse_row_norm(params.p, w) - 2.0 / params.p)
                < 1e-12,
                f"{name} inverse-W Parseval factor",
            )
        check(
            abs(
                encoded_real_coefficient_variance(params)
                - params.message_bound**2
                * 2.0
                / (3.0 * params.p * params.n * params.n)
            )
            < 1e-30,
            f"{name} encoded coefficient variance",
        )

        fused = estimate_gl_ship(params)
        legacy = estimate_gl_ship(params, arithmetic_mode="legacy")
        correlated = estimate_gl_ship(params, model="correlated")
        check(math.isfinite(fused.full_precision_bits), f"{name} finite result")
        check(fused.phase_wrap_margin_bits > 0.0, f"{name} phase margin")
        check(
            correlated.full_he_log2_variance >= fused.full_he_log2_variance,
            f"{name} correlated model monotonicity",
        )
        check(
            fused.full_he_log2_variance <= legacy.full_he_log2_variance,
            f"{name} fused-DD monotonicity",
        )
        check(
            legacy.stages["tree/level-1-tensor-rescale"]
            >= legacy.stages["tree/level-1-relinearization"],
            f"{name} rescale-before-relinearization floor",
        )
        check(
            fused.stages["tree/level-1-rescale"]
            <= legacy.stages["tree/level-1-tensor-rescale"],
            f"{name} relinearize-before-rescale floor",
        )
        check(
            fused.stages["half/HMux-moddown"]
            <= legacy.stages["half/HMux-moddown"],
            f"{name} fused HMux ModDown",
        )
        check(
            abs(
                legacy.stages["half/HMux-moddown"]
                - fused.stages["half/HMux-moddown"]
                - math.log2(2 * params.hmux_radix)
            )
            < 1e-12,
            f"{name} one ModDown per HMux stage",
        )

        noisier = estimate_gl_ship(replace(params, error_stddev=6.4))
        check(
            noisier.stc_phase_log2_variance >= fused.stc_phase_log2_variance,
            f"{name} encryption-noise monotonicity",
        )

        switch = dd_key_switch_noise(
            params, log_q=params.log_q, destination="dense", full_ring=False
        )
        more_aux = replace(
            params,
            log_p=params.log_p + 8,
            storage_bits=params.storage_bits + 8,
        )
        quieter_switch = dd_key_switch_noise(
            more_aux,
            log_q=more_aux.log_q,
            destination="dense",
            full_ring=False,
        )
        check(
            quieter_switch.log2_total_variance <= switch.log2_total_variance + 1e-12,
            f"{name} auxiliary-modulus monotonicity",
        )
        alternate_bbar = replace(params, bbar_bit=4)
        alternate_bbar.validate()
        exact_limb_switch = dd_key_switch_noise(
            alternate_bbar,
            log_q=alternate_bbar.log_q,
            destination="dense",
            full_ring=False,
        )
        check(
            abs(exact_limb_switch.log2_total_variance - switch.log2_total_variance)
            < 1e-12,
            f"{name} exact Bbar-limb invariance",
        )
        estimates[name] = fused
        print(
            f"  {name}: fused/legacy precision="
            f"{fused.full_precision_bits:.2f}/{legacy.full_precision_bits:.2f} "
            f"bits, phase margin={fused.phase_wrap_margin_bits:.2f} bits, "
            f"storage/security={params.storage_margin_bits:+d}/"
            f"{params.security_margin_bits:+d}"
        )

    # The paper omits its individual RNS primes.  Starting from the
    # reconstructed tree limbs, use only the already available output-depth
    # headroom and find the first uniform tree scale that reaches Table 4.
    expected_recovery = {"n512p17": 47, "n1024p17": 50}
    for name, expected_scale in expected_recovery.items():
        params = PRESETS[name]
        recovered = None
        for scale in range(
            params.tree_scale_bits,
            params.tree_scale_bits + params.tree_scale_headroom_bits + 1,
        ):
            candidate = replace(params, tree_scale_bits=scale)
            candidate.validate()
            estimate = estimate_gl_ship(candidate)
            if estimate.full_precision_bits >= params.paper_precision_bits:
                recovered = estimate
                break
        check(recovered is not None, f"{name} paper-precision recovery")
        check(
            recovered.params.tree_scale_bits == expected_scale,
            f"{name} minimum recovered tree scale",
        )
        print(
            f"  {name}: target recovered at tree scale "
            f"{recovered.params.tree_scale_bits} with "
            f"{recovered.full_precision_bits:.2f} bits"
        )

    n256 = PRESETS["n256p17"]
    n256_max = replace(
        n256,
        tree_scale_bits=n256.tree_scale_bits + n256.tree_scale_headroom_bits,
    )
    n256_estimate = estimate_gl_ship(n256_max)
    check(
        n256_estimate.full_precision_bits < n256.paper_precision_bits,
        "n256 reconstructed schedule must remain explicitly unresolved",
    )

    # The C++ regression uses deliberately small messages and accepts 0.08 for
    # the direct half path and 0.12 for the full path.  This is a schedule
    # calibration, not a fit: the estimator must independently stay within the
    # same guards while retaining a positive phase margin.
    toy = tfhepp_toy_profile()
    toy.validate()
    toy_estimate = estimate_gl_ship(toy)
    toy_half_sigma = 2.0 ** (0.5 * toy_estimate.half_he_log2_variance)
    toy_full_sigma = 2.0 ** (0.5 * toy_estimate.full_he_log2_variance)
    check(toy_estimate.half_total_error_bound < 0.08, "toy half-bootstrap guard")
    check(toy_estimate.full_total_error_bound < 0.12, "toy full-bootstrap guard")
    check(toy_estimate.phase_wrap_margin_bits > 0.0, "toy phase margin")
    check(
        toy_estimate.full_he_log2_variance >= toy_estimate.half_he_log2_variance,
        "toy full-path accumulation",
    )
    smaller_message = estimate_gl_ship(replace(toy, message_bound=0.055))
    check(
        abs(
            smaller_message.half_sine_error_bound
            * 8.0
            / toy_estimate.half_sine_error_bound
            - 1.0
        )
        < 1e-12,
        "SHIP cubic sine bound",
    )
    print(
        "  tfhepp-toy: "
        f"half/full HE sigma={toy_half_sigma:.3e}/{toy_full_sigma:.3e}, "
        f"D-sigma totals={toy_estimate.half_total_error_bound:.3e}/"
        f"{toy_estimate.full_total_error_bound:.3e}"
    )

    check(
        estimates["n512p17"].full_precision_bits
        > estimates["n256p17"].full_precision_bits,
        "higher-precision profile ordering",
    )
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
