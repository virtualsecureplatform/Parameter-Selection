#!/usr/bin/env python3
"""Build a public guessing-mask descriptor for bounded ring-attack studies.

The legacy adapters call the descriptor's anchor 'best'; it contains no cost
estimate or optimality claim. Secret block widths, sigma and modulus remain
those of the baseline. Only the public guessing mask/weight changes.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import prod
from pathlib import Path
from narrow_ckks_block_hybrid_grid import saturated_layout
from narrow_ckks_block_rotation_probability import coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--zeta', type=int, required=True)
    parser.add_argument('--saturated-blocks', type=int, required=True)
    parser.add_argument('--weight', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    baseline = json.loads(args.baseline.read_text())
    widths = baseline['block_widths']
    if len(widths) != baseline['h'] or sum(widths) != baseline['n']:
        parser.error('baseline block geometry mismatch')
    if not 0 <= args.weight <= baseline['h']:
        parser.error('invalid exact guessed weight')
    guessed = saturated_layout(widths, args.zeta, args.saturated_blocks)
    mass, support = coefficients(widths, guessed, args.weight)
    if not mass or not support:
        parser.error('unsupported zero-rotation guessed-weight stratum')
    probability = Fraction(mass, prod(widths))
    anchor = baseline['attacks']['bdd_mitm_hybrid']['best']
    result = dict(scope=__doc__, geometry_only=True, security_validated=False,
        complete_attack_optimization=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        baseline_report=str(args.baseline),
        baseline_report_sha256=hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
        n=baseline['n'], h=baseline['h'], sigma=baseline['sigma'],
        log_q=baseline['log_q'], block_widths=widths,
        attacks={'bdd_mitm_hybrid': {'best': dict(
            grid_beta=anchor['grid_beta'], grid_d=anchor['grid_d'],
            grid_zeta=args.zeta, guessed_weight=args.weight,
            guessed_per_block=guessed,
            layout=f'saturated-{args.saturated_blocks}-prefix',
            candidate_count=str(support),
            hit_probability_numerator=str(probability.numerator),
            hit_probability_denominator=str(probability.denominator))}})
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Public guessing geometry:', args.zeta, args.saturated_blocks,
          args.weight, 'positions:', sum(guessed))


if __name__ == '__main__':
    main()
