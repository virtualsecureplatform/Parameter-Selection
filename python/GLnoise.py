#!/usr/bin/env python3
"""Command-line noise estimator for TFHEpp's GL-SHIP bootstrap."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import sys


_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent
sys.path[:0] = [str(_REPO_ROOT), str(_SCRIPT_DIR)]

try:
    from python.noiseestimation.gl import (  # noqa: E402
        GLNoiseEstimate,
        GLNoiseParams,
        dd_key_switch_noise,
        estimate_as_dict,
        estimate_gl_ship,
        format_log2_variance,
    )
    from python.noiseestimation.params.gl import PRESETS  # noqa: E402
except ModuleNotFoundError as exc:
    if exc.name != "python":
        raise
    from noiseestimation.gl import (  # type: ignore[no-redef] # noqa: E402
        GLNoiseEstimate,
        GLNoiseParams,
        dd_key_switch_noise,
        estimate_as_dict,
        estimate_gl_ship,
        format_log2_variance,
    )
    from noiseestimation.params.gl import PRESETS  # type: ignore[no-redef] # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate coefficient and decoded-slot noise for TFHEpp's "
            "Double-Decomposition GL-SHIP bootstrap."
        )
    )
    parser.add_argument(
        "--preset",
        choices=("all", *sorted(PRESETS)),
        default="all",
        help="paper/TFHEpp profile to estimate",
    )
    parser.add_argument(
        "--model",
        choices=("independent", "correlated"),
        default="independent",
        help="independent variance sums or worst-aligned correlated sums",
    )
    parser.add_argument(
        "--arithmetic",
        choices=("fused-dd", "legacy"),
        default="fused-dd",
        help=(
            "fused-dd matches hybrid-RNS operation boundaries and updated "
            "TFHEpp; legacy models the former per-term ModDown and "
            "rescale-before-relinearization path"
        ),
    )
    parser.add_argument(
        "--masked-moddown",
        choices=("per-candidate", "fused"),
        default=None,
        help=(
            "override only the masked-column boundary; by default it follows "
            "--arithmetic"
        ),
    )
    parser.add_argument(
        "--D",
        type=float,
        default=6.0,
        help="Gaussian tail multiplier applied to the modeled HE standard deviation",
    )
    parser.add_argument("--q-bits", type=int, help="override total ciphertext log Q")
    parser.add_argument("--p-bits", type=int, help="override auxiliary log P")
    parser.add_argument("--q0-bits", type=int, help="override level-zero modulus bits")
    parser.add_argument("--gap-bits", type=int, help="override log2(gamma)")
    parser.add_argument("--stc-bits", type=int, help="override grouped StC limb bits")
    parser.add_argument("--x-scale-bits", type=int, help="override X transform scale")
    parser.add_argument("--w-scale-bits", type=int, help="override W transform scale")
    parser.add_argument(
        "--tree-scale-bits", type=int, help="override each of the five product levels"
    )
    parser.add_argument(
        "--optimize-tree-scale",
        action="store_true",
        help=(
            "select the smallest uniform tree scale that reaches the target "
            "using the schedule's existing output-depth headroom"
        ),
    )
    parser.add_argument("--primary-bit", type=int, help="override DD primary base bits")
    parser.add_argument("--bbar-bit", type=int, help="override DD auxiliary limb bits")
    parser.add_argument(
        "--storage-bits", type=int, help="override TFHEpp coefficient storage width"
    )
    parser.add_argument(
        "--error-stddev",
        type=float,
        help="absolute coefficient Gaussian standard deviation",
    )
    parser.add_argument(
        "--message-bound",
        type=float,
        help=(
            "absolute bound for each modeled real/imaginary component "
            "entering the SHIP phase map"
        ),
    )
    parser.add_argument(
        "--target-precision",
        type=float,
        help=(
            "required full-output precision in bits; defaults to the paper's "
            "measured precision for each preset"
        ),
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit machine-readable JSON instead of tables",
    )
    parser.add_argument(
        "--stages", action="store_true", help="print every modeled stage contribution"
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "exit nonzero on an exceeded storage/security/phase bound or "
            "precision below the selected target"
        ),
    )
    return parser.parse_args()


def apply_overrides(params: GLNoiseParams, args: argparse.Namespace) -> GLNoiseParams:
    overrides: dict[str, object] = {}
    mapping = {
        "q_bits": "log_q",
        "p_bits": "log_p",
        "q0_bits": "q0_bits",
        "gap_bits": "gap_bits",
        "stc_bits": "stc_bits",
        "x_scale_bits": "x_scale_bits",
        "w_scale_bits": "w_scale_bits",
        "tree_scale_bits": "tree_scale_bits",
        "primary_bit": "primary_bit",
        "bbar_bit": "bbar_bit",
        "storage_bits": "storage_bits",
        "error_stddev": "error_stddev",
        "message_bound": "message_bound",
    }
    for argument, field in mapping.items():
        value = getattr(args, argument)
        if value is not None:
            overrides[field] = value

    if "stc_bits" in overrides:
        stc_bits = int(overrides["stc_bits"])
        x_explicit = "x_scale_bits" in overrides
        w_explicit = "w_scale_bits" in overrides
        if not x_explicit and not w_explicit:
            overrides["x_scale_bits"] = stc_bits // 2
            overrides["w_scale_bits"] = stc_bits - stc_bits // 2
        elif x_explicit and not w_explicit:
            overrides["w_scale_bits"] = stc_bits - int(overrides["x_scale_bits"])
        elif w_explicit and not x_explicit:
            overrides["x_scale_bits"] = stc_bits - int(overrides["w_scale_bits"])

    result = replace(params, **overrides)
    result.validate()
    return result


def _format_error(value: float) -> str:
    if value == 0.0:
        return "0"
    return f"{value:.4e}"


def estimate_with_args(
    params: GLNoiseParams, args: argparse.Namespace
) -> GLNoiseEstimate:
    return estimate_gl_ship(
        params,
        model=args.model,
        arithmetic_mode=args.arithmetic,
        masked_moddown=args.masked_moddown,
        tail_bound=args.D,
    )


def optimize_tree_scale(
    params: GLNoiseParams, args: argparse.Namespace
) -> GLNoiseEstimate:
    target = (
        args.target_precision
        if args.target_precision is not None
        else params.paper_precision_bits
    )
    best = estimate_with_args(params, args)
    if best.full_precision_bits >= target:
        return best
    maximum = params.tree_scale_bits + params.tree_scale_headroom_bits
    for scale in range(params.tree_scale_bits + 1, maximum + 1):
        candidate = replace(params, tree_scale_bits=scale)
        candidate.validate()
        best = estimate_with_args(candidate, args)
        if best.full_precision_bits >= target:
            return best
    return best


def _margin_status(margin: float) -> str:
    if margin > 0.0:
        return "OK"
    if margin == 0.0:
        return "AT LIMIT"
    return "EXCEEDED"


def _target_precision(
    estimate: GLNoiseEstimate, requested_target: float | None
) -> float:
    if requested_target is not None:
        return requested_target
    return estimate.params.paper_precision_bits


def _screen_status(estimate: GLNoiseEstimate, target: float) -> str:
    if (
        estimate.params.storage_margin_bits < 0
        or estimate.params.security_margin_bits < 0
        or estimate.params.output_depth_margin_bits < 0
    ):
        return "INFEASIBLE"
    if (
        not math.isfinite(estimate.full_precision_bits)
        or estimate.phase_wrap_margin_bits < 0.0
        or estimate.full_total_error_bound >= 1.0
    ):
        return "INSUFFICIENT"
    if estimate.full_precision_bits < target:
        return "BELOW TARGET"
    return "OK"


def print_summary(
    estimates: list[GLNoiseEstimate], requested_target: float | None
) -> None:
    print("GL-SHIP noise estimate (TFHEpp coefficient-domain Double Decomposition)")
    print(
        "  epsilon_HE is modeled statistically; sine and q0-encoding bounds are "
        "reported separately."
    )
    first = estimates[0]
    print(
        f"  arithmetic={first.arithmetic_mode}, aggregation={first.model}, "
        f"masked ModDown={first.masked_moddown}, tail D={first.tail_bound:g}."
    )
    print(
        "  q0/tree scale defaults are reconstructed from the matched SHIP profiles; "
        "override them for a concrete modulus schedule."
    )
    if requested_target is None:
        print("  target defaults to the paper's hybrid-RNS measured precision.")
    else:
        print(f"  requested full-output precision target={requested_target:g} bits.")
    print()
    header = (
        "profile       N  logQ logP storage sec-margin phase-margin "
        "HE-bound total-bound precision target delta status"
    )
    print(header)
    print("-" * len(header))
    for estimate in estimates:
        params = estimate.params
        target = _target_precision(estimate, requested_target)
        delta = estimate.full_precision_bits - target
        print(
            f"{params.tag:<11} {params.ring_degree:5d} {params.log_q:5d} "
            f"{params.log_p:4d} {params.storage_bits:7d} "
            f"{params.security_margin_bits:10d} "
            f"{estimate.phase_wrap_margin_bits:12.2f} "
            f"{_format_error(estimate.full_he_error_bound):>8} "
            f"{_format_error(estimate.full_total_error_bound):>10} "
            f"{estimate.full_precision_bits:9.2f} "
            f"{target:6.2f} {delta:+6.2f} "
            f"{_screen_status(estimate, target)}"
        )


def print_detail(
    estimate: GLNoiseEstimate,
    *,
    stages: bool,
    requested_target: float | None,
) -> None:
    params = estimate.params
    target = _target_precision(estimate, requested_target)
    print(f"\n{params.tag}")
    print(
        f"  algebra: N=2*n*phi(p)={params.ring_degree}, n={params.n}, "
        f"p={params.p}, phi(p)={params.phi}"
    )
    print(
        f"  moduli: logQ={params.log_q}, logP={params.log_p}, "
        f"log(PQ)={params.log_pq}, storage={params.storage_bits} "
        f"(margin {params.storage_margin_bits:+d})"
    )
    print(
        f"  security ceiling: {params.security_limit_log_pq} bits "
        f"(margin {params.security_margin_bits:+d}, "
        f"{_margin_status(params.security_margin_bits)})"
    )
    print(
        f"  schedule: q0={params.q0_bits}, gap=2^{params.gap_bits}, "
        f"StC={params.x_scale_bits}+{params.w_scale_bits}, "
        f"tree={params.tree_depth}x{params.tree_scale_bits}, "
        f"output logQ={params.output_log_q} "
        f"(outside-depth margin {params.output_depth_margin_bits:+d}, "
        f"tree-scale headroom {params.tree_scale_headroom_bits})"
    )
    print(
        f"  DD: primary={params.primary_bit} bits, Bbar={params.bbar_bit} bits, "
        f"covers={params.primary_cover_bits}/{params.bbar_cover_bits} of "
        f"{params.storage_bits} storage bits; paper dnum={params.dnum} "
        "(RNS metadata only)"
    )
    print(
        f"  noise inputs: coefficient sigma={params.error_stddev:g}, "
        f"uniform slot range=[-{params.message_bound:g},"
        f"{params.message_bound:g}], "
        f"dense-secret variance={params.dense_secret_variance:g}"
    )
    print(
        f"  selection: h={params.sparse_hamming_weight}, theta={params.theta}, "
        f"HMux={params.hmux_stages} stages radix {params.hmux_radix}, "
        f"masks={params.masked_column_count} "
        f"({params.average_candidates_per_sparse_term:.2f}/term)"
    )
    print(
        f"  phase-wrap {estimate.tail_bound:g}-sigma margin: "
        f"half={estimate.half_phase_wrap_margin_bits:.2f}, "
        f"full={estimate.full_phase_wrap_margin_bits:.2f} bits"
    )
    print(
        "  direct half bootstrap: "
        f"HE={_format_error(estimate.half_he_error_bound)}, "
        f"sine={_format_error(estimate.half_sine_error_bound)}, "
        f"q0-round={_format_error(estimate.half_quantization_error_bound)}, "
        f"total={_format_error(estimate.half_total_error_bound)} "
        f"({estimate.half_precision_bits:.2f} bits)"
    )
    print(
        f"  full bootstrap: HE={_format_error(estimate.full_he_error_bound)}, "
        f"sine={_format_error(estimate.full_sine_error_bound)}, "
        f"q0-round={_format_error(estimate.full_quantization_error_bound)}, "
        f"total={_format_error(estimate.full_total_error_bound)} "
        f"({estimate.full_precision_bits:.2f} bits)"
    )
    print(
        f"  paper measured precision: {params.paper_precision_bits:.2f} bits; "
        f"screening delta={estimate.precision_delta_from_paper:+.2f} bits"
    )
    print(
        f"  DD precision target: {target:.2f} bits; "
        f"margin={estimate.full_precision_bits - target:+.2f} bits "
        f"({_screen_status(estimate, target)})"
    )

    dense_switch = dd_key_switch_noise(
        params, log_q=params.log_q, destination="dense", full_ring=False
    )
    print(
        "  top-level DD switch coefficient log2(V): "
        f"eval-key={dense_switch.log2_eval_key_variance:.2f}, "
        f"moddown={dense_switch.log2_moddown_rounding_variance:.2f}, "
        f"total={dense_switch.log2_total_variance:.2f}, "
        f"rows={dense_switch.primary_rows}"
    )
    rescale_floor = estimate.stages["tree/level-1-rescale"]
    relin_floor = estimate.stages["tree/level-1-relinearization"]
    print(
        "  level-1 product decoded log2(V): "
        f"rescale={rescale_floor:.2f}, relinearization={relin_floor:.2f}"
    )

    if stages:
        print("  stage contributions (log2 variance in the stage's stated domain):")
        for name, value in estimate.stages.items():
            print(f"    {name:<38} {format_log2_variance(value):>9}")


def strict_failure(estimate: GLNoiseEstimate, target: float) -> bool:
    params = estimate.params
    return (
        params.storage_margin_bits < 0
        or params.security_margin_bits < 0
        or estimate.phase_wrap_margin_bits < 0.0
        or not math.isfinite(estimate.full_precision_bits)
        or estimate.full_precision_bits < target
    )


def main() -> int:
    args = parse_args()
    if args.optimize_tree_scale and args.tree_scale_bits is not None:
        print(
            "GLnoise.py: error: --optimize-tree-scale and --tree-scale-bits "
            "are mutually exclusive",
            file=sys.stderr,
        )
        return 2
    if args.target_precision is not None and (
        not math.isfinite(args.target_precision) or args.target_precision < 0.0
    ):
        print(
            "GLnoise.py: error: target precision must be a finite, nonnegative value",
            file=sys.stderr,
        )
        return 2
    names = list(PRESETS) if args.preset == "all" else [args.preset]
    try:
        configured = [apply_overrides(PRESETS[name], args) for name in names]
        estimates = [
            optimize_tree_scale(params, args)
            if args.optimize_tree_scale
            else estimate_with_args(params, args)
            for params in configured
        ]
    except (OverflowError, ValueError) as exc:
        print(f"GLnoise.py: error: {exc}", file=sys.stderr)
        return 2

    if args.json:
        payload = []
        for estimate in estimates:
            target = _target_precision(estimate, args.target_precision)
            item = estimate_as_dict(estimate)
            item["screening"] = {
                "target_precision_bits": target,
                "precision_margin_bits": estimate.full_precision_bits - target,
                "status": _screen_status(estimate, target),
            }
            payload.append(item)
        print(json.dumps(payload[0] if len(payload) == 1 else payload, indent=2))
    else:
        print_summary(estimates, args.target_precision)
        for estimate in estimates:
            print_detail(
                estimate,
                stages=args.stages,
                requested_target=args.target_precision,
            )
        below = [
            estimate.params.tag
            for estimate in estimates
            if _screen_status(
                estimate, _target_precision(estimate, args.target_precision)
            )
            != "OK"
        ]
        if below:
            print(
                "\nDD screen warning: the selected schedule is below its precision "
                f"target for {', '.join(below)}."
            )
        print(
            "\nCaution: this is a parameter-screening model, not an RLWE security "
            "estimate or a proof of bootstrap correctness."
        )

    if args.strict and any(
        strict_failure(
            estimate, _target_precision(estimate, args.target_precision)
        )
        for estimate in estimates
    ):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
