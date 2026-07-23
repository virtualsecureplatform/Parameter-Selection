#!/usr/bin/env python3
"""GL-SHIP profiles from ePrint 2026/811 and the TFHEpp DD backend."""

from __future__ import annotations

try:
    from python.noiseestimation.gl import GLNoiseParams
except ModuleNotFoundError as exc:
    if exc.name != "python":
        raise
    from noiseestimation.gl import GLNoiseParams  # type: ignore[no-redef]


# q0 and the product-tree scales are reconstructed from the matched SHIP
# profiles in ePrint 2025/784, Table 2.  ePrint 2026/811 publishes only the GL
# total Q/P sizes, grouped StC limb, and gap.  They are intentionally exposed as
# CLI overrides by GLnoise.py.
N256P17 = GLNoiseParams(
    tag="n256p17",
    n=256,
    p=17,
    log_q=180,
    log_p=34,
    q0_bits=25,
    gap_bits=5,
    stc_bits=26,
    x_scale_bits=13,
    w_scale_bits=13,
    tree_scale_bits=25,
    outside_multiplicative_depth=1,
    outside_scale_bits=20,
    dnum=8,
    theta=4,
    window_width=384,
    masked_column_count=848,
    primary_bit=16,
    bbar_bit=8,
    storage_bits=256,
    security_limit_log_pq=214,
    paper_precision_bits=4.60,
)

N512P17 = GLNoiseParams(
    tag="n512p17",
    n=512,
    p=17,
    log_q=338,
    log_p=92,
    q0_bits=48,
    gap_bits=11,
    stc_bits=37,
    x_scale_bits=18,
    w_scale_bits=19,
    # Four wide primary rows mirror the paper's dnum=4 partition count while
    # the independently decomposed evaluation-key rows remain signed 16-bit
    # data.  This is a TFHEpp DD mapping; dnum itself describes the paper's RNS
    # basis, not an 85-bit gadget base.  A 49-bit tree scale gives a useful
    # margin over the paper's reported 14.94-bit precision.
    tree_scale_bits=49,
    outside_multiplicative_depth=1,
    outside_scale_bits=37,
    dnum=4,
    theta=8,
    window_width=576,
    masked_column_count=1504,
    primary_bit=85,
    bbar_bit=16,
    storage_bits=448,
    security_limit_log_pq=430,
    paper_precision_bits=14.94,
)

N1024P17 = GLNoiseParams(
    tag="n1024p17",
    n=1024,
    p=17,
    log_q=641,
    log_p=220,
    q0_bits=50,
    gap_bits=11,
    stc_bits=39,
    x_scale_bits=19,
    w_scale_bits=20,
    tree_scale_bits=46,
    outside_multiplicative_depth=8,
    outside_scale_bits=39,
    dnum=4,
    theta=16,
    window_width=1056,
    masked_column_count=2880,
    primary_bit=16,
    bbar_bit=7,
    storage_bits=896,
    security_limit_log_pq=868,
    paper_precision_bits=15.88,
)


PRESETS = {
    profile.tag: profile for profile in (N256P17, N512P17, N1024P17)
}

for _profile in PRESETS.values():
    _profile.validate()
