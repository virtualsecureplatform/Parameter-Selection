#!/usr/bin/env python3
"""Screen/search average-case narrow CKKS DD refresh parameters; ordinary Python."""
import argparse
from dataclasses import asdict, replace
import json
from noiseestimation.narrow_ckks import NarrowParams, estimate, minimum_square_output_bits, retained_square_range


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    default = NarrowParams()
    for name, value in asdict(default).items():
        parser.add_argument("--" + name.replace("_", "-"), type=type(value), default=value)
    parser.add_argument("--model", choices=("independent", "correlated"), default="independent")
    parser.add_argument("--search", action="store_true")
    parser.add_argument("--target-precision", type=float, default=11)
    parser.add_argument("--max-key-modulus", type=int, default=218)
    parser.add_argument("--max-cache-gib", type=float, default=32)
    parser.add_argument("--streamed", action="store_true")
    parser.add_argument("--relinearize-before-rescale", action="store_true")
    parser.add_argument("--auxiliary-switch-bits", type=int, default=0)
    parser.add_argument("--selector-chunk-size", type=int, default=0,
                        help="Bound each transient selector NTT cache; zero reports a whole group")
    parser.add_argument("--batch-columns", type=int, default=1,
                        help="Input/frequency columns for chunked modular sum storage")
    parser.add_argument("--min-scale-bits", type=int, default=20)
    parser.add_argument("--max-scale-bits", type=int, default=58)
    parser.add_argument("--min-plain-bits", type=int, default=20)
    parser.add_argument("--max-plain-bits", type=int, default=40)
    parser.add_argument("--retained-component-bound", type=float)
    parser.add_argument("--retained-guard-bits", type=int, default=3)
    parser.add_argument("--range-only", action="store_true",
                        help="Report retained DD square range without the standard sine noise model")
    args = parser.parse_args()
    if not (0 < args.min_scale_bits <= args.max_scale_bits and
            0 < args.min_plain_bits <= args.max_plain_bits <= 50):
        parser.error("invalid scale/plain search bounds")
    if args.auxiliary_switch_bits < 0:
        parser.error("auxiliary switch bits must be nonnegative")
    if args.selector_chunk_size < 0 or args.batch_columns < 1:
        parser.error("selector chunk size must be nonnegative and batch columns positive")
    p = NarrowParams(**{name: getattr(args, name) for name in asdict(default)})
    if args.range_only:
        if args.retained_component_bound is None:
            parser.error("--range-only requires --retained-component-bound")
        p.validate()
        result = retained_square_range(p.output_bits, p.scale_bits,
            args.retained_component_bound, args.retained_guard_bits)
        selector = p.plain_bits + p.depth * p.scale_bits + result["minimum_output_bits"]
        print(json.dumps({"parameters": asdict(p), "range": result,
            "minimum_selector_bits": selector,
            "minimum_packed_bits": selector + p.plain_bits + p.n.bit_length() - 1,
            "security_validated": False}, indent=2))
        return 0
    if not args.search:
        result = estimate(p, args.model, args.streamed, args.relinearize_before_rescale, args.auxiliary_switch_bits,
                          args.selector_chunk_size, args.batch_columns)
        if args.retained_component_bound is not None:
            result["retained_square_range"] = retained_square_range(p.output_bits,
                p.scale_bits, args.retained_component_bound, args.retained_guard_bits)
        print(json.dumps({"parameters": asdict(p), "estimate": result}, indent=2))
        return 0
    candidates = []
    for scale in range(args.min_scale_bits, args.max_scale_bits + 1):
        for plain in range(args.min_plain_bits, args.max_plain_bits + 1):
            output = max(p.source_bits + p.s2c_bits + 1, 2 * scale + 2)
            if args.retained_component_bound is not None:
                output = max(output, minimum_square_output_bits(scale,
                    args.retained_component_bound, args.retained_guard_bits))
            candidate = replace(p, scale_bits=scale, plain_bits=plain,
                                selector_bits=plain + p.depth * scale + output)
            if candidate.packed_bits > args.max_key_modulus:
                continue
            result = estimate(candidate, args.model, args.streamed, args.relinearize_before_rescale, args.auxiliary_switch_bits,
                              args.selector_chunk_size, args.batch_columns)
            if result["maximum_evaluation_key_bits"] > args.max_key_modulus:
                continue
            if args.retained_component_bound is not None:
                result["retained_square_range"] = retained_square_range(candidate.output_bits,
                    candidate.scale_bits, args.retained_component_bound, args.retained_guard_bits)
            if (result["eight_rms_plus_sine_screen"] <= 2**-args.target_precision and
                (result["chunked_selector_memory"]["ntt_chunk_bytes"] if args.selector_chunk_size else
                 result["ntt_cache_peak_bytes"]) <= args.max_cache_gib * 2**30):
                candidates.append({"parameters": asdict(candidate), "estimate": result})
    candidates.sort(key=lambda c: (c["estimate"]["packed_bits"],
                                  c["estimate"]["ntt_cache_bytes_upper_bound"]))
    print(json.dumps({"target_precision": args.target_precision,
        "scope": "Average-case refresh screen; eight RMS is not a tail proof; FFT/security excluded",
        "candidates": candidates[:10], "candidate_count": len(candidates)}, indent=2))
    return 0 if candidates else 1


if __name__ == "__main__":
    raise SystemExit(main())
