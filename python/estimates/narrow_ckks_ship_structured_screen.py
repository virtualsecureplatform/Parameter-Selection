#!/usr/bin/env python3
"""Compare public narrow block geometry with SHIP's Section 5.1 heuristic.

SHIP uses (4*w)^(0.4*h) scalar products and asks for more than 2^115
scalar products at N >= 8192. Its experimental w is max(175,N/(2*h)).
For disjoint ragged narrow blocks, use the exact signed support product
as a geometry-only analogue. Transferring the exponent to that distribution
is an assumption, not a proved attack cost or security qualification.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def screen(n, h=31):
    if n < 8192 or n & (n - 1) or not 1 < h < n:
        raise ValueError('expected power-of-two N>=8192 and 1<h<N')
    g = (n + h - 1) // h
    widths = [min(g, n - j*g) for j in range(h)]
    if min(widths) <= 0:
        raise ValueError('ragged schedule contains an empty block')
    support = math.prod(2*w for w in widths)
    entropy = math.log2(support)
    scalar_bits = 0.4 * entropy
    ship_w = max(175, n/(2*h))
    return dict(n=n, h=h, block_widths=widths,
                signed_key_support_count=str(support),
                signed_key_support_log2=entropy,
                geometry_analogue_scalar_product_log2=scalar_bits,
                ship_scalar_product_threshold_log2=115,
                meets_ship_geometry_heuristic=scalar_bits > 115,
                ship_experimental_w=ship_w,
                ship_experimental_scalar_product_log2=0.4*h*math.log2(4*ship_w),
                security_validated=False,
                inference='Apply SHIP exponent to exact disjoint-block support; not a proved attack or lower bound.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    sources = ['refs/978-3-031-91131-6_14.pdf',
               'TFHEpp/include/ckks/narrow_bootstrapping.hpp']
    cases = [screen(n) for n in [8192, 16384, 32768]]
    # Integer support products independently agree with log-of-factor sums.
    for case in cases:
        assert sum(case['block_widths']) == case['n']
        assert abs(sum(math.log2(2*w) for w in case['block_widths']) -
                   case['signed_key_support_log2']) < 1e-10
    args.output.write_text(json.dumps(dict(
        scope=__doc__, source_sha256={p: hashlib.sha256((repo/p).read_bytes()).hexdigest()
                                    for p in sources},
        cases=cases, security_validated=False,
        omitted=['lattice/hybrid attacks with public block constraints',
                 'ring structure and correlations', 'noise and evaluation-key moduli',
                 'circular/KDM assumptions', 'algorithmic constants and memory tradeoffs'],
        references=['Local SHIP Section 5.1, regularly-spaced nonzero coefficients',
                    'https://eprint.iacr.org/2021/216.pdf']), indent=2)+'\n')


if __name__ == '__main__':
    main()
