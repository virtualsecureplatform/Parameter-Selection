#!/usr/bin/env python3
"""Independently check saved joint-decoder probabilities and grid minima.

Decimal arithmetic checks the posterior formulas, not the scalar GSA decoder
or the attack cost formula. Passing does not qualify structured-key security.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(path, posterior_path):
    grid = json.loads(path.read_text())
    posterior = json.loads(posterior_path.read_text())
    if (grid['complete'] is not True or grid['security_validated'] is not False
            or grid['complete_attack_optimization'] is not False):
        raise ValueError('incomplete grid or unsupported security claim')
    if grid['probability_report_sha256'] != sha(posterior_path):
        raise ValueError('posterior provenance mismatch')
    if (posterior['complete'] is not True or grid['n'] != posterior['n']
            or grid['weight'] != posterior['guessed_weight']):
        raise ValueError('posterior scope mismatch')
    histogram = {int(k): int(v) for k, v in
                 posterior['matching_rotation_histogram'].items()}
    samples = sum(histogram.values())
    if min(histogram) < 1 or min(histogram.values()) < 1 or samples != grid['histogram_samples']:
        raise ValueError('invalid posterior histogram')
    event_fraction = Fraction(posterior['expected_matching_rotations'])
    inverse_fraction = sum((Fraction(v, k) for k, v in histogram.items()), Fraction()) / samples
    checked = 0
    unavailable = []
    minima = {}
    for model, scenario in grid['scenarios'].items():
        points = scenario['points']
        coordinates = [(p['beta'], p['d']) for p in points]
        if len(set(coordinates)) != len(coordinates):
            raise ValueError('duplicate grid point')
        if len(points) != len(grid['beta_offsets']) * len(grid['d_offsets']):
            raise ValueError('missing grid point (validator requires a fully admissible grid)')
        costs = {}
        for point in points:
            if 'decode_marginal' not in point:
                if point.get('finite') is not False or not point.get('reason'):
                    raise ValueError('missing probability without an explicit unavailable record')
                unavailable.append(dict(model=model, **point))
                continue
            p = Decimal(point['decode_marginal'])
            if not p.is_finite() or not 0 < p <= 1:
                raise ValueError('invalid decoder probability')
            with localcontext() as context:
                # Extra precision makes direct subtraction independent of
                # the producer's binary recursive geometric-sum algorithm.
                context.prec = 200 + max(0, -p.adjusted())
                event = Decimal(event_fraction.numerator) / event_fraction.denominator
                inverse = Decimal(inverse_fraction.numerator) / inverse_fraction.denominator
                common = event * inverse * p
                independent = event * sum(
                    (Decimal(count) * (1 - (1-p)**k) / k
                     for k, count in histogram.items()), Decimal()) / samples
                upper = min(event * inverse, event * p)
                for field, expected in (
                        ('joint_success_common', common),
                        ('joint_success_independent', independent),
                        ('joint_success_upper_bound', upper)):
                    actual = Decimal(point[field])
                    if not actual.is_finite() or abs(actual / expected - 1) > Decimal('1e-130'):
                        raise ValueError('joint probability mismatch: ' + field)
                    checked += 1
            for name, cost in point['costs'].items():
                bits = cost['bits']
                if type(bits) not in (int, float) or not math.isfinite(bits) or bits <= 0:
                    raise ValueError('nonfinite or missing cost')
                candidate = (bits, point['beta'], point['d'])
                if name not in costs or bits < costs[name][0]:
                    costs[name] = candidate
        for name, best in scenario['best'].items():
            if (best['bits'], best['beta'], best['d']) != costs[name]:
                raise ValueError('saved minimum does not match grid')
        if set(costs) != set(scenario['best']):
            raise ValueError('missing saved minimum')
        minima[model] = costs
    return dict(report=str(path), report_sha256=sha(path),
                checked_probabilities=checked, minima=minima,
                unavailable_points=unavailable,
                unavailable_decoder_probabilities_validated=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid', type=Path, action='append', required=True)
    parser.add_argument('--posterior', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reports = [validate(path, args.posterior) for path in args.grid]
    result = dict(scope=__doc__, validator_source_sha256=sha(Path(__file__)),
                  posterior_sha256=sha(args.posterior), reports=reports,
                  passed=True, security_validated=False)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print('checked joint probabilities', sum(r['checked_probabilities'] for r in reports))


if __name__ == '__main__':
    main()
