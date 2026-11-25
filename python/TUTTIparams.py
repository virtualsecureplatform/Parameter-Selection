#!/usr/bin/env python3
"""Budget, noise screen and search for TFHEpp's TUTTI CKKS refresh; ordinary Python.

The lattice screens of the two published keys need Sage:
`sage -python python/estimates/tutti_security.py --profile ll13 --output report.json`.
"""
import argparse
from dataclasses import asdict, replace
import json

from noiseestimation.block_binary_meet_lwe import combinatorial_report
from noiseestimation.tutti import PROFILES, TuttiParams, estimate, published_caps, search


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", choices=sorted(PROFILES),
                        help="start from a TFHEpp schedule; other options override it")
    for name, value in asdict(TuttiParams()).items():
        option = "--" + name.replace("_", "-")
        if isinstance(value, bool):
            parser.add_argument(option, action=argparse.BooleanOptionalAction)
        else:
            parser.add_argument(option, type=type(value))
    parser.add_argument("--model", choices=("independent", "correlated"), default="independent")
    parser.add_argument("--search", action="store_true",
                        help="choose the scales under --max-dense-bits")
    parser.add_argument("--max-dense-bits", type=int,
                        help="largest modulus at which the dense key may be published")
    parser.add_argument("--application-levels", type=int, default=1)
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument("--combinatorial", action="store_true",
                        help="add Odlyzko and Meet-LWE costs of the sparse key (about 5 s)")
    args = parser.parse_args()
    p = PROFILES[args.profile] if args.profile else TuttiParams()
    overrides = {name: getattr(args, name) for name in asdict(p)
                 if getattr(args, name) is not None}
    p = replace(p, **overrides)
    if args.search:
        if args.max_dense_bits is None:
            parser.error("--search requires --max-dense-bits")
        found = search(p, args.max_dense_bits, args.application_levels, model=args.model)
        if not found:
            parser.error("no schedule fits the dense-key budget")
        result = {"maximum_dense_bits": args.max_dense_bits, "candidates": found[:args.top],
                  "security_validated": False}
    else:
        try:
            result = estimate(p, args.model)
        except ValueError as error:
            parser.error(str(error))
        result["parameters"] = asdict(p)
        result["published_caps"] = published_caps(p)
    if args.combinatorial:
        result["sparse_key_combinatorial"] = combinatorial_report(
            p.sparse_degree, p.blocks, p.base_bits, max(p.sigma, 1e-9))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
