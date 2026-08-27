#!/usr/bin/env sage -python
"""Algorithm-3 parameter and source-security screen for shallow bootstrap.

This script records the concrete Section 7.2 schedule from ePrint 2026/1730:
the sparse LWE source has n=1450, h=29, q=512 and sigma=3.2; the Binary-NTT
RLWE product uses N=4096 and a 105-bit-to-50-bit modulus boundary after three
then two multiplication layers.

The LWE estimator can screen the sparse LWE source.  It cannot establish the
paper's Binary-NTT RLWE assumption, so that portion is intentionally reported
as an assumption rather than an estimated security claim.
"""

import argparse
import importlib
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
estimator = importlib.import_module(".estimator", "lattice-estimator")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lwe-n", type=int, default=1450)
    parser.add_argument("--lwe-h", type=int, default=29)
    parser.add_argument("--lwe-qbits", type=int, default=9)
    parser.add_argument("--lwe-sigma", type=float, default=3.2)
    parser.add_argument("--ring-n", type=int, default=4096)
    parser.add_argument("--initial-qbits", type=int, default=105)
    parser.add_argument("--post-switch-qbits", type=int, default=50)
    parser.add_argument("--windows", type=int, nargs="+", default=[8, 4])
    parser.add_argument("--pbc-copies", type=int, default=3)
    parser.add_argument("--pbc-slack", type=int, default=3)
    parser.add_argument("--rlwe-sigma", type=float, default=0.75)
    parser.add_argument("--key-switch-digits", type=int, default=4)
    parser.add_argument("--double-decomposition", action="store_true")
    parser.add_argument("--skip-security", action="store_true")
    parser.add_argument("--ring-proxy-security", action="store_true",
                        help="also run a ternary RLWE-as-LWE proxy at the initial modulus")
    return parser.parse_args()


def log2_product_depth(factor_count):
    return math.ceil(math.log2(max(1, factor_count)))


args = arguments()
if args.lwe_n <= 0 or args.lwe_h <= 0 or args.lwe_h > args.lwe_n:
    raise SystemExit("invalid sparse LWE dimensions")
if args.ring_n <= 0 or args.ring_n & (args.ring_n - 1):
    raise SystemExit("ring-n must be a power of two")
if any(window <= 1 or window & (window - 1) for window in args.windows):
    raise SystemExit("every multiplication window must be a power of two > 1")
if args.pbc_copies <= 0 or args.pbc_slack < 0:
    raise SystemExit("invalid PBC parameters")
if args.rlwe_sigma <= 0 or args.key_switch_digits <= 0:
    raise SystemExit("invalid RLWE noise or key-switch digits")

window_layers = sum(int(math.log2(window)) for window in args.windows)
pbc_leaves = args.lwe_h + args.pbc_slack
tree_layers = log2_product_depth(pbc_leaves)
if pbc_leaves > math.prod(args.windows):
    raise SystemExit("PBC leaves do not fit in the requested product windows")
print("Algorithm 3 low-depth schedule")
print(f"  LWE: n={args.lwe_n} h={args.lwe_h} q=2^{args.lwe_qbits} sigma={args.lwe_sigma}")
print(f"  Binary-NTT RLWE: N={args.ring_n} Q=2^{args.initial_qbits} -> 2^{args.post_switch_qbits}")
print(f"  PBC: c={args.pbc_copies} k=h+{args.pbc_slack}={pbc_leaves}")
print(f"  multiplication windows={args.windows} ({window_layers} product layers)")
print(f"  balanced tree depth for {pbc_leaves} PBC factors={tree_layers}")
print("  NTT/INTT inside relinearization-free tree=0")
print("  Binary-NTT RLWE security=paper assumption (not modeled by standard LWE estimator)")
if args.double_decomposition:
    torus_bits = 128
    embedding_bits = torus_bits - args.initial_qbits
    print("  prospective TFHEpp Double-Decomposition mapping:")
    print(f"    torus=2^{torus_bits}, native-Q embedding scale=2^{embedding_bits}")
    print("    DD limbs=8, limb bits=16 (the existing lvl3 128-bit layout)")
    print(f"    encoded RLWE sigma'=0.75*2^{embedding_bits}")
    print("    requires Algorithm-3 two-component DD multiplication; not executable")

# Linear boundary-cost diagnostic for Algorithm 3's intended native-Q
# boundary.  This deliberately excludes ciphertext-ciphertext multiplication
# noise and therefore is not an end-to-end correctness estimate.  Use
# shallowboot_noise.py for the PBC aggregation and the N*V1*V2 product term.
bucket_entries = math.ceil(args.pbc_copies * args.lwe_n / pbc_leaves)
factor_stddev = args.rlwe_sigma * math.sqrt(bucket_entries)
high_group_error = args.windows[0] * factor_stddev
ks_basebits = math.ceil(args.initial_qbits / args.key_switch_digits)
ks_error = (args.rlwe_sigma * (2**ks_basebits) *
            math.sqrt(args.ring_n * args.key_switch_digits / 3.0))
scaled_ks_error = ks_error * 2 ** (args.post_switch_qbits -
                                   args.initial_qbits)
rounding_error = 1.0 + args.ring_n / 2.0
low_inputs = math.ceil(pbc_leaves / args.windows[0])
low_error = low_inputs * (scaled_ks_error + rounding_error)
plaintext_scale = 2 ** (args.post_switch_qbits - 3)
margin_bits = math.log2(plaintext_scale / max(1.0, low_error))
print("  native-Q linear boundary diagnostic (not an end-to-end noise bound):")
print(f"    PBC entries/bucket <= {bucket_entries}; factor sigma <= {factor_stddev:.2f}")
print(f"    high-group factor error <= {high_group_error:.2f}")
print(f"    boundary KSK: d={args.key_switch_digits}, base≈2^{ks_basebits}, "
      f"scaled sigma <= 2^{math.log2(max(scaled_ks_error, 1e-300)):.1f}")
print(f"    low-stage error bound <= 2^{math.log2(low_error):.1f}; "
      f"plaintext scale=2^{args.post_switch_qbits - 3}; "
      f"margin≈{margin_bits:.1f} bits")
print("    WARNING: excludes product noise; run estimates/shallowboot_noise.py")

if not args.skip_security:
    source = estimator.lwe_parameters.LWEParameters(
        n=args.lwe_n,
        q=2**args.lwe_qbits,
        Xs=estimator.nd.SparseBinary(args.lwe_h, n=args.lwe_n),
        Xe=estimator.nd.DiscreteGaussian(stddev=args.lwe_sigma),
        m=float("inf"),
        tag="algorithm3-sparse-source-proxy",
    )
    # Some estimator attack implementations have a hard dimension bound.  Do
    # not turn that implementation limit into an opaque KeyError; retain every
    # successfully evaluated attack and label an incomplete screen below.
    results = estimator.LWE.estimate(
        source, red_cost_model=estimator.RC.BDGL16, catch_exceptions=True
    )
    security = min(math.log2(result["rop"]) for result in results.values())
    print(f"  sparse-LWE source proxy security={security:.1f} bits")
    if security < 128:
        raise SystemExit("sparse LWE source proxy does not meet 128 bits")
    if args.ring_proxy_security:
        ring_proxy = estimator.lwe_parameters.LWEParameters(
            n=args.ring_n,
            q=2**args.initial_qbits,
            Xs=estimator.nd.Ternary,
            Xe=estimator.nd.DiscreteGaussian(stddev=args.rlwe_sigma),
            m=float("inf"),
            tag="algorithm3-qh-ss-rlwe-as-lwe-proxy",
        )
        ring_results = estimator.LWE.estimate(
            ring_proxy, red_cost_model=estimator.RC.BDGL16,
            catch_exceptions=True
        )
        ring_security = min(
            math.log2(result["rop"]) for result in ring_results.values()
        )
        print(f"  QH-SS-RLWE-as-LWE proxy security={ring_security:.1f} bits")
        if ring_security < 128:
            raise SystemExit("ring proxy does not meet 128 bits")
