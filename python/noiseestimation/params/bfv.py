#!/usr/bin/env python3
"""BFV parameter presets for the average-case noise estimator."""

from __future__ import annotations

import math

try:
    from python.noiseestimation.bfv import BFVParams
except ModuleNotFoundError as exc:
    if exc.name != "python":
        raise
    from noiseestimation.bfv import BFVParams


TFHEPP_LVL3SIMD_P = 114_689
TFHEPP_LVL3SIMD_Q_BITS = 128
TFHEPP_LVL3SIMD_ALPHA_LOG2 = -105.0
TFHEPP_LVL5_P = 786_433
TFHEPP_LVL5_Q_BITS = 448
TFHEPP_LVL5_ALPHA_LOG2 = -425.0


def paper_u3(
    *,
    nbit: int = 13,
    q_bits: int = 149,
    t: int = 65_537,
    error_std: float = 3.19,
) -> BFVParams:
    """Generic public-key BFV preset matching the U3 examples in 600.pdf."""
    return BFVParams(
        nbit=nbit,
        t=t,
        q_bits=q_bits,
        secret_variance=2.0 / 3.0,
        encryption_u_variance=2.0 / 3.0,
        error_log2_std=math.log2(error_std),
        key_switch="hybrid-rns",
        rns_digits=4,
        hybrid_omega=4,
        fresh="public",
        tag="paper-u3",
    )


def openfhe_paper(
    *,
    nbit: int = 13,
    q_bits: int = 60,
    t: int = 65_537,
    error_std: float = 3.19,
) -> BFVParams:
    """OpenFHE-style BFV parameters used in 600.pdf validation tables."""
    params = paper_u3(nbit=nbit, q_bits=q_bits, t=t, error_std=error_std)
    return BFVParams(
        nbit=params.nbit,
        t=params.t,
        q_bits=params.q_bits,
        secret_variance=params.secret_variance,
        encryption_u_variance=params.encryption_u_variance,
        error_log2_std=params.error_log2_std,
        key_switch="hybrid-rns",
        rns_digits=4,
        hybrid_omega=4,
        fresh="public",
        correction_fit=params.correction_fit,
        tag="openfhe-paper",
    )


def tfhepp_lvl3simd_base(
    *,
    q_bits: int = TFHEPP_LVL3SIMD_Q_BITS,
    alpha_log2: float = TFHEPP_LVL3SIMD_ALPHA_LOG2,
) -> BFVParams:
    """TFHEpp lvl3simdparam, plaintext modulus p=114689."""
    return BFVParams(
        nbit=12,
        t=TFHEPP_LVL3SIMD_P,
        q_bits=q_bits,
        secret_variance=2.0 / 3.0,
        encryption_u_variance=2.0 / 3.0,
        error_log2_std=q_bits + alpha_log2,
        key_switch="hybrid-rns",
        rns_digits=4,
        hybrid_omega=8,
        fresh="symmetric",
        tag="tfhepp-lvl3simd-base",
    )


def tfhepp_lvl3simd_boot(
    *,
    q_bits: int = TFHEPP_LVL3SIMD_Q_BITS,
    alpha_log2: float = TFHEPP_LVL3SIMD_ALPHA_LOG2,
) -> BFVParams:
    """TFHEpp BFV bootstrap PrimePower2Param, plaintext modulus p^2."""
    return BFVParams(
        nbit=12,
        t=TFHEPP_LVL3SIMD_P * TFHEPP_LVL3SIMD_P,
        q_bits=q_bits,
        secret_variance=2.0 / 3.0,
        encryption_u_variance=2.0 / 3.0,
        error_log2_std=q_bits + alpha_log2,
        key_switch="hybrid-rns",
        rns_digits=4,
        hybrid_omega=8,
        fresh="symmetric",
        tag="tfhepp-lvl3simd-boot",
    )


def tfhepp_lvl5_base(
    *,
    q_bits: int = TFHEPP_LVL5_Q_BITS,
    alpha_log2: float = TFHEPP_LVL5_ALPHA_LOG2,
) -> BFVParams:
    """TFHEpp lvl5param multi-limb BFV scaffold, plaintext modulus p=786433."""
    return BFVParams(
        nbit=14,
        t=TFHEPP_LVL5_P,
        q_bits=q_bits,
        secret_variance=2.0 / 3.0,
        encryption_u_variance=2.0 / 3.0,
        error_log2_std=q_bits + alpha_log2,
        key_switch="hybrid-rns",
        rns_digits=4,
        hybrid_omega=28,
        fresh="symmetric",
        tag="tfhepp-lvl5-base",
    )


def tfhepp_lvl5_boot(
    *,
    q_bits: int = TFHEPP_LVL5_Q_BITS,
    alpha_log2: float = TFHEPP_LVL5_ALPHA_LOG2,
) -> BFVParams:
    """TFHEpp lvl5param BFV bootstrap PrimePower2Param, plaintext modulus p^2."""
    return BFVParams(
        nbit=14,
        t=TFHEPP_LVL5_P * TFHEPP_LVL5_P,
        q_bits=q_bits,
        secret_variance=2.0 / 3.0,
        encryption_u_variance=2.0 / 3.0,
        error_log2_std=q_bits + alpha_log2,
        key_switch="hybrid-rns",
        rns_digits=4,
        hybrid_omega=28,
        fresh="symmetric",
        tag="tfhepp-lvl5-boot",
    )


PRESETS = {
    "openfhe-paper": openfhe_paper,
    "paper-u3": paper_u3,
    "tfhepp-lvl3simd-base": tfhepp_lvl3simd_base,
    "tfhepp-lvl3simd-boot": tfhepp_lvl3simd_boot,
    "tfhepp-lvl5-base": tfhepp_lvl5_base,
    "tfhepp-lvl5-boot": tfhepp_lvl5_boot,
}
