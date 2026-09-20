#!/usr/bin/env python3
"""Conditional whole-run failure analysis for scalar TFHE/BFV switching."""

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

from noiseestimation.bfv_switch import parameters_for_width, analyse, search, tune_16bit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--products", type=int, default=256)
    parser.add_argument("--operand-bits", type=int, default=8)
    parser.add_argument("--legacy-profile", action="store_true",
                        help="reproduce the original untuned wider-operand experiment")
    parser.add_argument("--bitwise-forward", action="store_true",
                        help="model the original one-bit forward conversion")
    parser.add_argument("--failure-bits", type=int, default=40)
    parser.add_argument("--search", action="store_true")
    parser.add_argument("--tune-16bit", action="store_true",
                        help="reproduce the bounded 16-bit ring/gadget search")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--benchmark-report", type=Path,
                        help="check the modeled parameters against a successful measured BFV run")
    args = parser.parse_args()
    if not 1 <= args.operand_bits <= 32:
        parser.error("--operand-bits must be in 1..32")
    params = parameters_for_width(args.operand_bits, legacy=args.legacy_profile,
                                  products=args.products, target_failure_bits=args.failure_bits)
    if args.bitwise_forward:
        params = replace(params, forward_digit_bits=1)
    result = search(params) if args.search else analyse(params)
    if args.tune_16bit:
        if args.operand_bits != 16 or args.legacy_profile or args.search:
            parser.error("--tune-16bit requires --operand-bits 16 and no other search/profile override")
        result = tune_16bit(args.products, args.failure_bits)
    if args.benchmark_report:
        measured = json.loads(args.benchmark_report.read_text())
        if measured.get("products", 256) != params.products or measured.get("operand_bits", 8) != args.operand_bits:
            parser.error("model/benchmark product count or operand width mismatch")
        run = measured["runs"]["bfv"]
        if run.get("bfv_forward_digit_bits", 1) != result["parameters"]["forward_digit_bits"]:
            parser.error("model/benchmark forward digit width mismatch")
        if not run["passed"]:
            parser.error("BFV benchmark did not pass")
        mapping = {"ring_n": "bfv_ring_dimension", "qbits": "bfv_ciphertext_bits",
                   "ring_alpha_log2": "bfv_alpha_log2", "ring_levels": "bfv_gadget_levels",
                   "ring_basebit": "bfv_gadget_basebit", "to_io_levels": "bfv_to_tfhe_levels",
                   "to_io_basebit": "bfv_to_tfhe_basebit", "io_n": "tfhe_io_n",
                   "io_alpha_log2": "tfhe_io_alpha_log2", "io_levels": "tfhe_io_levels",
                   "io_basebit": "tfhe_io_basebit", "half_n": "tfhe_half_n",
                   "half_alpha_log2": "tfhe_half_alpha_log2", "to_half_levels": "tfhe_to_half_levels",
                   "to_half_basebit": "tfhe_to_half_basebit"}
        for model, runtime in mapping.items():
            if result["parameters"][model] != run[runtime]:
                parser.error(f"model/benchmark parameter mismatch: {model}")
        if run["bfv_plaintext_modulus"] != 2**result["parameters"]["plaintext_bits"]:
            parser.error("BFV plaintext modulus mismatch")
        if run["evaluation_key_bytes"] != result["evaluation_key_bytes"]:
            parser.error("evaluation-key size mismatch")
        if run["bfv_dd_levels"] != 8 or run["bfv_dd_basebit"] != 16:
            parser.error("analysis requires full 8x16-bit DD limb coverage")
        observed = {
            "forward": run["max_forward_integer_unit_error"],
            "product": run["max_product_integer_unit_error"],
        }
        for stage, error in observed.items():
            if not math.isfinite(error) or error < 0:
                parser.error(f"invalid observed {stage} error")
            if error > result["envelope"][f"{stage}_error"] * 2**params.plaintext_bits:
                parser.error(f"observed {stage} error exceeds the conditional envelope")
        result["observed_integer_unit_errors"] = observed
        result["observed_errors_within_envelope"] = True
        result["measured_source_sha256"] = measured["source_sha256"]
        result["benchmark_report"] = str(args.benchmark_report)
    rendered = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
    return int(not result["meets_target_conditionally"])


if __name__ == "__main__":
    raise SystemExit(main())
