#!/usr/bin/env python3
"""Prospective nine-level TFHEpp DD schedules within Table-2 HT key budgets.

Run from this directory's parent with ordinary Python. These are engineering
screens; encrypted measurements and structured/ring/KDM security are required.
"""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from noiseestimation.narrow_ckks import NarrowParams, estimate, retained_unit_product_chain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--lean-dd", action="store_true",
                        help="Screen primary18/secondary8 digits with the unchanged HT level schedules")
    args = parser.parse_args()
    if args.output is None:
        args.output = Path(__file__).with_name("narrow_ckks_throughput_lean_candidates.json") if args.lean_dd else Path(__file__).with_suffix(".json")
    cases = [
        ("ht14-inspired-dd", 438, 4, 58, 23, 0.19,
         NarrowParams(n=16384, selector_bits=383, scale_bits=25, plain_bits=23,
                      source_bits=28, input_scale_bits=26, s2c_bits=23)),
        ("ht15-inspired-dd", 881, 16, 212, 39, 0.13,
         NarrowParams(n=32768, selector_bits=669, scale_bits=46, plain_bits=39,
                      source_bits=50, input_scale_bits=43, s2c_bits=39)),
    ]
    reports = []
    for name, budget, precision, auxiliary, product_scale, bound, params in cases:
        if args.lean_dd:
            name = name.replace("inspired-dd", "lean-dd")
            params = replace(params, primary_bits=18, bbar_bits=8)
        refresh = estimate(params, streamed=True, relinearize_before_rescale=True,
                           auxiliary_switch_bits=auxiliary,
                           selector_chunk_size=40 if args.lean_dd else 0)
        chain = retained_unit_product_chain(params, 9, product_scale, bound, auxiliary)
        passes = (refresh["maximum_evaluation_key_bits"] <= budget and
                  chain["guarded_range_passed"] and
                  chain["eight_rms_plus_sine_screen"] < 2**-precision)
        reports.append({"profile": name, "parameters": asdict(params),
                        "auxiliary_switch_bits": auxiliary,
                        "key_modulus_budget": budget, "target_precision_bits": precision,
                        "refresh": refresh, "retained_unit_chain": chain,
                        "engineering_screen_passed": passes,
                        "encrypted_validation_scope": "Configured full ring size; N=512 diagnostic results are separate",
                        "encrypted_validation_passed": False, "security_validated": False})
    args.output.write_text(json.dumps({"scope": __doc__, "profiles": reports}, indent=2) + "\n")
    for report in reports:
        print(report["profile"], report["engineering_screen_passed"],
              report["retained_unit_chain"]["eight_rms_plus_sine_screen"])
    return 0 if all(r["engineering_screen_passed"] for r in reports) else 1


if __name__ == "__main__":
    raise SystemExit(main())
