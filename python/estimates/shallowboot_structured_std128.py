#!/usr/bin/env sage -python
"""Security and noise screen for TFHEpp's structured shallow bootstrap.

Run from Parameter-Selection/python with the repository's Sage environment:

    sage -python estimates/shallowboot_structured_std128.py
    sage -python estimates/shallowboot_structured_std128.py --source-n 512 --source-h 64

The source estimate uses a structured one-hot binary secret but models it as a
uniform fixed-weight sparse-binary secret.  The optional shorter setting
therefore screens the proposed post-key-switch blind-rotation secret; it does
not analyze the key switch itself.  The target estimate is an RLWE-as-LWE
proxy for TFHEpp's N=1024, Q=2^32 ternary accumulator.
"""

import argparse
import importlib
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
estimator = importlib.import_module(".estimator", "lattice-estimator")


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-n", type=int, default=1024)
    parser.add_argument("--source-h", type=int, default=64)
    parser.add_argument("--source-qbits", type=int, default=9)
    parser.add_argument("--source-sigma", type=float, default=3.2)
    parser.add_argument("--target-n", type=int, default=1024)
    parser.add_argument("--target-qbits", type=int, default=32)
    parser.add_argument("--target-alpha-bits", type=int, default=25)
    parser.add_argument("--gadget-basebit", type=int, default=6)
    parser.add_argument("--gadget-levels", type=int, default=3)
    parser.add_argument("--skip-security", action="store_true")
    return parser.parse_args()


def estimate(label, parameters):
    results = estimator.LWE.estimate(
        parameters, red_cost_model=estimator.RC.BDGL16
    )
    security = min(math.log2(result["rop"]) for result in results.values())
    print(f"{label}: estimated classical security = {security:.1f} bits")
    return security


def shallow_pbs_variance(source_n, source_h, target_n, alpha_bk,
                          gadget_basebit, gadget_levels):
    """Algorithm-2 structured-BR variance in normalized target-torus units.

    This is the TFHEnoise.py blind-rotation model adapted to h encrypted
    bucket products.  A bucket contains n/h RGSWs, so its key error variance
    grows by n/h while the expensive external-product and rounding terms occur
    only h times.  The model is an average-case variance screen and excludes
    FFT numerical error and a source-return key switch.
    """
    bucket_width = source_n / source_h
    base = 2**gadget_basebit
    beta = base / 2
    epsilon = 1 / (2 * base**gadget_levels)
    key_error = (source_h * 2 * gadget_levels * target_n * beta**2 *
                 (bucket_width * alpha_bk**2))
    decomposition = source_h * (target_n + 1) * epsilon**2
    return key_error + decomposition, key_error, decomposition


def log2_failure_probability(variance):
    # TFHEpp's conventional Boolean margin is 1/16 of the torus.  The erfc
    # call is stable for this range under Sage's MPFR-backed real numbers.
    from sage.all import RealField, erfc

    real = RealField(256)
    probability = erfc(real(1) / 16 / (2 * real(variance)).sqrt())
    if probability == 0:
        return float("-inf")
    return float(probability.log(2))


args = parse_arguments()
if args.source_n <= 0 or args.source_h <= 0:
    raise SystemExit("source dimension and Hamming weight must be positive")
if args.source_n % args.source_h:
    raise SystemExit("one-hot structured source requires source-n divisible by source-h")
if args.source_h > args.source_n:
    raise SystemExit("Hamming weight cannot exceed source dimension")

bucket_width = args.source_n // args.source_h
entropy_bits = args.source_h * math.log2(bucket_width)
target_alpha = 2**(-args.target_alpha_bits)
variance, key_error, decomposition = shallow_pbs_variance(
    args.source_n, args.source_h, args.target_n, target_alpha,
    args.gadget_basebit, args.gadget_levels)
print("structured shallow PBS noise screen")
print(f"  source n={args.source_n} h={args.source_h} bucket_width={bucket_width}")
print(f"  one-hot secret entropy={entropy_bits:.1f} bits")
print(f"  external products={args.source_h}")
print(f"  key-error log2(variance)={math.log2(key_error):.2f}")
print(f"  decomposition log2(variance)={math.log2(decomposition):.2f}")
print(f"  output log2(variance)={math.log2(variance):.2f}")
print(f"  Boolean failure log2 bound={log2_failure_probability(variance):.2f}")

if not args.skip_security:
    source = estimator.lwe_parameters.LWEParameters(
        n=args.source_n,
        q=2**args.source_qbits,
        Xs=estimator.nd.SparseBinary(args.source_h, n=args.source_n),
        Xe=estimator.nd.DiscreteGaussian(stddev=args.source_sigma),
        m=float("inf"),
        tag="shallowboot-structured-source-proxy",
    )
    target_proxy = estimator.lwe_parameters.LWEParameters(
        n=args.target_n,
        q=2**args.target_qbits,
        Xs=estimator.nd.Ternary,
        Xe=estimator.nd.DiscreteGaussian(
            stddev=2**(args.target_qbits - args.target_alpha_bits)),
        m=float("inf"),
        tag="shallowboot-tfhepp-accumulator-proxy",
    )

    source_bits = estimate("source (uniform sparse-binary proxy)", source)
    target_bits = estimate("accumulator (RLWE-as-LWE proxy)", target_proxy)
    if min(source_bits, target_bits) < 128:
        raise SystemExit("candidate does not meet the 128-bit screen")
