#!/usr/bin/env python3
"""Command-line CLPX scheme-switch noise estimator."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent
sys.path[:0] = [str(_REPO_ROOT), str(_SCRIPT_DIR)]

try:
    from python.noiseestimation.clpx import (  # noqa: E402
        NEG_INF,
        estimate_clpx_multiplication_depth,
        estimate_clpx_to_tlwes,
        estimate_tlwes_to_clpx,
        format_log2,
        pbs_input_margin_log2,
    )
    from python.noiseestimation.params import clpx as params  # noqa: E402
except ModuleNotFoundError as exc:
    if exc.name != "python":
        raise
    from noiseestimation.clpx import (  # type: ignore[no-redef] # noqa: E402
        NEG_INF,
        estimate_clpx_multiplication_depth,
        estimate_clpx_to_tlwes,
        estimate_tlwes_to_clpx,
        format_log2,
        pbs_input_margin_log2,
    )
    from noiseestimation.params import clpx as params  # type: ignore[no-redef] # noqa: E402


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Estimate TFHEpp CLPX scheme-switching noise by composing TFHE operations."
    )
    ap.add_argument(
        "--direction",
        choices=(
            "all",
            "both",
            "tlwes-to-clpx",
            "clpx-to-tlwes",
            "multiplication",
            "switched-multiplication",
        ),
        default="both",
    )
    ap.add_argument("--validbit", type=int, default=8)
    ap.add_argument("--num-multi", type=int, default=4)
    ap.add_argument(
        "--shift",
        type=int,
        help="TFHEpp template shift; paper shiftnum is shift+1",
    )
    ap.add_argument(
        "--shiftnum",
        type=int,
        help="paper shiftnum value; mutually exclusive with --shift",
    )
    ap.add_argument(
        "--w",
        type=int,
        help="Delta_b truncation width; default is the full torus width unless --paper-ss2clpx is used",
    )
    ap.add_argument(
        "--paper-ss2clpx",
        action="store_true",
        help="use Nagai et al. TFHE->CLPX settings: base=2, w=20, shiftnum=5",
    )
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--numdigit", type=int, default=4)
    ap.add_argument("--basebit", type=int, default=2)
    ap.add_argument(
        "--max-mults",
        type=int,
        default=16,
        help="maximum CLPXMult count to print in multiplication-depth mode",
    )
    ap.add_argument(
        "--mult-chain",
        choices=("fresh", "square"),
        default="fresh",
        help="fresh multiplies the accumulated ciphertext by a fresh switched ciphertext; square repeatedly squares",
    )
    ap.add_argument(
        "--message-bound",
        type=float,
        help="absolute bound on CLPX plaintext digit values during multiplication; default is p-1",
    )
    ap.add_argument(
        "--D",
        type=float,
        default=6.0,
        help="Gaussian tail parameter used for multiplication-depth OK/FAIL margins",
    )
    ap.add_argument(
        "--input-log2-var",
        type=float,
        help="override input variance for the selected direction; use separate runs for different overrides",
    )
    return ap.parse_args()


def resolve_switch_args(args: argparse.Namespace) -> None:
    if args.shift is not None and args.shiftnum is not None:
        raise SystemExit("--shift and --shiftnum are mutually exclusive")
    if args.shiftnum is not None:
        if args.shiftnum <= 0:
            raise SystemExit("--shiftnum must be positive")
        shift = args.shiftnum - 1
    elif args.shift is not None:
        if args.shift < 0:
            raise SystemExit("--shift must be non-negative")
        shift = args.shift
    elif args.paper_ss2clpx:
        shift = 4
    else:
        shift = 0

    if args.w is not None and args.w <= 0:
        raise SystemExit("--w must be positive")
    w = args.w
    if args.paper_ss2clpx and w is None:
        w = 20

    args.effective_shift = shift
    args.effective_w = w
    args.effective_shiftnum = shift + 1
    args.switch_iksP = params.lvl10param
    if args.paper_ss2clpx:
        args.switch_iksP = params.lvl1hparam
        args.switch_bkP = params.SS2CLPXlvlh2param
        args.switch_sskP = params.SS2CLPXlvl22param
        args.switch_clpxP = params.SS2CLPXlvl2param
        args.switch_label = "SS2CLPX paper path"
    else:
        args.switch_bkP = params.lvl02param
        args.switch_sskP = params.lvl22param
        args.switch_clpxP = params.lvl2param
        args.switch_label = "TFHEpp default CLPX path"


def _fmt_prob(log2_prob: float) -> str:
    if log2_prob == NEG_INF:
        return "2^-inf"
    return f"2^{log2_prob:.2f}"


def print_tlwes_to_clpx(args: argparse.Namespace) -> None:
    input_var = None
    if args.input_log2_var is not None and args.direction != "both":
        input_var = 2.0 ** args.input_log2_var
    est = estimate_tlwes_to_clpx(
        validbit=args.validbit,
        num_multi=args.num_multi,
        shift=args.effective_shift,
        w=args.effective_w,
        iksP=args.switch_iksP,
        bkP=args.switch_bkP,
        sskP=args.switch_sskP,
        input_variance=input_var,
    )
    print("TLWES2CLPXIKS")
    print(
        "  TFHEpp path: "
        f"{args.switch_iksP.domainP.__name__} TLWEs -> "
        f"{args.switch_iksP.targetP.__name__} IKS -> "
        f"{args.switch_bkP.__name__} many-LUT PBS -> "
        f"{args.switch_bkP.targetP.__name__} packing"
    )
    print(
        f"  parameter path: {args.switch_label}; CLPX base={args.switch_clpxP.plain_modulus}; "
        f"shift={args.effective_shift}; shiftnum={args.effective_shiftnum}; "
        f"w={args.effective_w or 'full'}"
    )
    print(f"  input log2(V):              {format_log2(est.input_variance)}")
    print(f"  after IKS log2(V):          {format_log2(est.iks_variance)}")
    iks_margin = pbs_input_margin_log2(
        args.switch_bkP, est.iks_variance, num_out=args.num_multi
    )
    print(f"  PBS bin margin:             {iks_margin:.2f} bits")
    print(f"  PBS output log2(V):         {format_log2(est.pbs_variance)}")
    print(
        f"  packing inputs:             {est.temp_tlwe_count} TLWEs, "
        f"max {est.max_terms_per_temp} PBS terms/TLWE"
    )
    print(f"  max packing input log2(V):  {format_log2(est.max_packed_input_variance)}")
    print(f"  CLPX coeff log2(V):         {format_log2(est.packed_variance)}")
    print(f"  CLPX digit-value log2(V):   {format_log2(est.clpx_value_variance)}")
    print(f"  estimated digit fail:       {_fmt_prob(est.log2_failure)}")


def print_clpx_to_tlwes(args: argparse.Namespace) -> None:
    input_var = None
    if args.input_log2_var is not None and args.direction != "both":
        input_var = 2.0 ** args.input_log2_var
    est = estimate_clpx_to_tlwes(
        validbit=args.validbit,
        batch_size=args.batch_size,
        numdigit=args.numdigit,
        basebit=args.basebit,
        input_variance=input_var,
    )
    internal_margin = pbs_input_margin_log2(
        params.lvlh2param, est.max_internal_pbs_input_variance, num_out=2
    )
    final_margin = pbs_input_margin_log2(
        params.lvlh1param, est.max_final_pbs_input_variance
    )
    print("CLPX2TLWESIKSanybit")
    print("  TFHEpp path: lvl2 CLPX -> lvlhalf/lvl2 PBS chain -> HomDecomp -> lvl1 PBS")
    print(f"  input CLPX coeff log2(V):   {format_log2(est.input_variance)}")
    print(f"  lvlh2 PBS output log2(V):   {format_log2(est.pbs02_variance)}")
    print(f"  sumpra log2(V):             {format_log2(est.sumpra_variance)}")
    print(f"  max HomDecomp sum log2(V):  {format_log2(est.max_homdecomp_sum_variance)}")
    print(f"  max internal PBS in log2(V): {format_log2(est.max_internal_pbs_input_variance)}")
    print(f"  internal PBS bin margin:    {internal_margin:.2f} bits")
    print(f"  max final PBS in log2(V):   {format_log2(est.max_final_pbs_input_variance)}")
    print(f"  final PBS bin margin:       {final_margin:.2f} bits")
    print(f"  produced TLWEs:             {est.produced_tlwes}")
    print(f"  output TLWE log2(V):        {format_log2(est.output_variance)}")
    print(f"  output decrypt fail:        {_fmt_prob(est.log2_failure)}")


def print_multiplication(args: argparse.Namespace) -> None:
    switch_est = None
    input_var = None
    if args.input_log2_var is not None and args.direction in {
        "multiplication",
        "switched-multiplication",
    }:
        input_var = 2.0 ** args.input_log2_var
        source = f"manual --input-log2-var={args.input_log2_var:g}"
    else:
        switch_est = estimate_tlwes_to_clpx(
            validbit=args.validbit,
            num_multi=args.num_multi,
            shift=args.effective_shift,
            w=args.effective_w,
            iksP=args.switch_iksP,
            bkP=args.switch_bkP,
            sskP=args.switch_sskP,
        )
        input_var = switch_est.packed_variance
        source = (
            "TLWES2CLPXIKS output "
            f"(validbit={args.validbit}, num_multi={args.num_multi}, "
            f"shift={args.effective_shift}, w={args.effective_w or 'full'})"
        )
    est = estimate_clpx_multiplication_depth(
        initial_coefficient_variance=input_var,
        validbit=args.validbit,
        num_multi=args.num_multi,
        shift=args.effective_shift,
        w=args.effective_w,
        max_multiplications=args.max_mults,
        chain=args.mult_chain,
        P=args.switch_clpxP,
        iksP=args.switch_iksP,
        bkP=args.switch_bkP,
        sskP=args.switch_sskP,
        message_bound=args.message_bound,
        d=args.D,
        model="paper" if args.paper_ss2clpx else "tfhepp",
    )
    print("CLPXMult depth estimate")
    print("  TFHEpp path: TRLWEMultWithoutRelinerizationCLPX + Relinearization")
    print(
        f"  parameter path: {args.switch_label}; CLPX base={args.switch_clpxP.plain_modulus}; "
        f"shift={args.effective_shift}; shiftnum={args.effective_shiftnum}; "
        f"w={args.effective_w or 'full'}"
    )
    if args.paper_ss2clpx:
        print("  multiplication noise model: Nagai et al. Eq. (44)-(48)")
    else:
        print("  caveat: approximate digit-value model with fixed plaintext bound")
    print(f"  initial source: {source}")
    if switch_est is not None:
        if args.paper_ss2clpx:
            print("  direct post-switch CLPXMult path: TFHEpp SS2CLPX base-2 test passes")
        else:
            print(
                "  caveat: default TFHEpp TLWES2CLPXIKS output is not a supported "
                "direct CLPXMult input; use --paper-ss2clpx for the paper path"
            )
        print(
            f"  switch output coeff log2(V): {format_log2(switch_est.packed_variance)}; "
            f"digit fail {_fmt_prob(switch_est.log2_failure)}"
        )
    print(f"  chain={est.chain} message_bound={est.message_bound:g} D={est.d:g}")
    print(f"  initial CLPX coeff log2(V): {format_log2(est.initial_coefficient_variance)}")
    if not args.paper_ss2clpx:
        print(f"  relin coeff log2(V):        {format_log2(est.relin_coefficient_variance)}")
    header = (
        "mults  log2(CLPX V)  log2(C2T V)  margin  log2(p_fail)  status"
        if args.paper_ss2clpx
        else "mults  log2(coeff V)  log2(value V)  margin  log2(p_fail)  status"
    )
    print(header)
    print("-" * len(header))
    for step in est.steps:
        print(
            f"{step.multiplication_count:5d} "
            f"{format_log2(step.coefficient_variance):>14} "
            f"{format_log2(step.digit_value_variance):>14} "
            f"{step.margin_bits:7.2f} "
            f"{_fmt_prob(step.log2_failure):>13}  "
            f"{step.status}"
        )
    if switch_est is not None and args.paper_ss2clpx:
        support_label = "supported post-switch multiplications"
    elif switch_est is not None:
        support_label = "supported bounded-screen multiplications"
    else:
        support_label = "supported multiplications"
    print(f"  {support_label}:  {est.supported_multiplications}")


def main() -> int:
    args = parse_args()
    resolve_switch_args(args)
    print("CLPX scheme-switch noise estimate")
    print("  source: TFHEpp/include/bfv-clpx.hpp operation sequence")
    print("  note: PBS output noise is estimated separately from PBS input-bin margins")
    print(
        f"  validbit={args.validbit} num_multi={args.num_multi} "
        f"shift={args.effective_shift} shiftnum={args.effective_shiftnum} "
        f"w={args.effective_w or 'full'} "
        f"batch_size={args.batch_size} numdigit={args.numdigit} basebit={args.basebit}"
    )
    if args.direction == "all":
        args.direction = "both"
        run_multiplication = True
    else:
        run_multiplication = args.direction in {"multiplication", "switched-multiplication"}
    if args.input_log2_var is not None and args.direction == "both":
        raise SystemExit("--input-log2-var requires a single --direction")

    if args.direction in {"both", "tlwes-to-clpx"}:
        print_tlwes_to_clpx(args)
    if args.direction == "both":
        print()
    if args.direction in {"both", "clpx-to-tlwes"}:
        print_clpx_to_tlwes(args)
    if run_multiplication:
        if args.direction == "both":
            print()
        print_multiplication(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
