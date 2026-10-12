#!/usr/bin/env python3
"""Bounded public-block hybrid search with exact guessed-weight strata.

Candidate counts and hit probabilities match one independent signed position
per public block. GSA uses scalar residual variance; heterogeneous coordinates,
ring structure and correlated evaluation keys remain outside this heuristic.
This is a bounded attack search, never a security certificate.
"""
import argparse
import datetime
from fractions import Fraction
import itertools
import json
from pathlib import Path
import sys

from narrow_ckks_block_hybrid import (cost_record,
    estimator_revision, guessing_counts, hash_file, layout, multiply_linear)


def saturated_layout(widths, zeta, count, short_first=False):
    """Fully guess count public blocks, then balance remaining coordinates."""
    if not 0 <= count < len(widths) or not 0 <= zeta < sum(widths):
        raise ValueError('invalid saturated guessing geometry')
    order = sorted(range(len(widths)), key=lambda j: (widths[j], j)) if short_first else list(range(len(widths)))
    full, remaining = order[:count], order[count:]
    used = sum(widths[j] for j in full)
    if used > zeta:
        raise ValueError('saturated blocks exceed guessing dimension')
    result = [0] * len(widths)
    for j in full:
        result[j] = widths[j]
    balanced = layout([widths[j] for j in remaining], zeta - used, 'balanced')
    for j, value in zip(remaining, balanced):
        result[j] = value
    assert sum(result) == zeta and all(0 <= g <= w for g, w in zip(result, widths))
    return result


def exact_strata(widths, guessed):
    cumulative = guessing_counts(widths, guessed)
    previous_support, previous_probability = 0, Fraction(0)
    for weight, support, probability in cumulative:
        count, mass = support - previous_support, probability - previous_probability
        previous_support, previous_probability = support, probability
        if not count or not mass:
            continue
        masses = [1]
        for width, selected in zip(widths, guessed):
            masses = multiply_linear(masses, width - selected, selected)
        accepted = masses[weight]
        variances = []
        for j, (width, selected) in enumerate(zip(widths, guessed)):
            if selected == width:
                variances.append(None)
                continue
            others = [1]
            for k, (w, g) in enumerate(zip(widths, guessed)):
                if k != j:
                    others = multiply_linear(others, w - g, g)
            variances.append(Fraction(others[weight] if weight < len(others) else 0, accepted))
        norm = sum((w - g) * v for w, g, v in zip(widths, guessed, variances) if v is not None)
        assert norm == len(widths) - weight
        yield weight, count, mass, norm, variances


def self_test():
    cases = 0
    for widths in ([2, 3, 1], [1, 1], [3]):
        secrets = list(itertools.product(*[[(p, s) for p in range(w) for s in (-1, 1)] for w in widths]))
        for guessed in itertools.product(*(range(w + 1) for w in widths)):
            total_probability = Fraction(0)
            for weight, count, probability, norm, variances in exact_strata(widths, guessed):
                accepted = [s for s in secrets if sum(p < g for (p, _), g in zip(s, guessed)) == weight]
                projections = {tuple((p, sign) if p < g else None for (p, sign), g in zip(s, guessed)) for s in accepted}
                assert len(projections) == count
                assert probability == Fraction(len(accepted), len(secrets))
                assert all(sum(p >= g for (p, _), g in zip(s, guessed)) == norm for s in accepted)
                for j, (w, g) in enumerate(zip(widths, guessed)):
                    for p in range(g, w):
                        assert variances[j] == Fraction(sum(s[j][0] == p for s in accepted), len(accepted))
                total_probability += probability
            assert total_probability == 1
            cases += 1
        for zeta in range(sum(widths)):
            for count in range(len(widths)):
                for short_first in (False, True):
                    order = sorted(range(len(widths)), key=lambda j: (widths[j], j)) if short_first else list(range(len(widths)))
                    if sum(widths[j] for j in order[:count]) <= zeta:
                        guessed = saturated_layout(widths, zeta, count, short_first)
                        assert sum(guessed) == zeta
                        assert all(guessed[j] == widths[j] for j in order[:count])
    print('exact-weight block counts, probabilities and residual moments passed:', cases, 'partitions', flush=True)


def offsets(value):
    try:
        result = sorted(set(int(v) for v in value.split(',')))
    except ValueError as error:
        raise argparse.ArgumentTypeError('comma-separated integer offsets required') from error
    if not result or len(result) > 11:
        raise argparse.ArgumentTypeError('between 1 and 11 distinct offsets required')
    return result


def run(args):
    sys.path.insert(0, str(args.estimator.resolve()))
    from estimator import LWE, ND, RC
    from estimator.lwe_primal import PrimalHybrid
    from sage.all import QQ, log, sqrt
    from narrow_ckks_security import configure_max_beta
    configure_max_beta(args.max_beta)
    baseline = json.loads(args.baseline_report.read_text())
    revision = estimator_revision(args.estimator.resolve())
    if baseline['estimator_revision'] != revision:
        raise ValueError('baseline estimator revision mismatch')
    n, h, q, sigma = (baseline[k] for k in ('n', 'h', 'log_q', 'sigma'))
    widths = baseline['block_widths']
    if sum(widths) != n or len(widths) != h or min(widths) < 1:
        raise ValueError('invalid baseline block geometry')
    width = (n + h - 1) // h
    if widths != [min(width, n - j * width) for j in range(h)]:
        raise ValueError('baseline differs from public ceil-width generator')
    report = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), scope=__doc__,
        source_sha256=hash_file(__file__), helper_source_sha256=hash_file(Path(__file__).with_name('narrow_ckks_block_hybrid.py')),
        estimator_revision=revision, baseline_report=str(args.baseline_report),
        baseline_report_sha256=hash_file(args.baseline_report), n=n, h=h, log_q=q, sigma=sigma,
        block_widths=widths, beta_offsets=args.beta_offsets, zeta_offsets=args.zeta_offsets,
        d_offsets=args.d_offsets, max_beta=args.max_beta, guessed_weight_model='exact-weight',
        saturated_blocks=args.saturated_blocks,
        geometry_model='exact residual norm with scalar average variance', self_test_passed=True,
        attacks={}, bounded_grid_complete=False, complete_attack_optimization=False, security_validated=False,
        limitations=['Only specified offsets/layouts searched; no global optimum claim',
            'Scalar GSA geometry omits heterogeneous variances and higher-order conditional correlations',
            'Estimator Babai/MITM probabilities and square-root speedup remain heuristic',
            'Zero-residual strata require separate exhaustive verification and are skipped',
            'Ring-specific and circular/correlated-key attacks remain unmodeled'])
    def save():
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    save()
    for attack in args.attacks:
        mitm = attack == 'bdd_mitm_hybrid'
        anchor = baseline['attacks'][attack]['uniform_baseline']
        if anchor['bits'] is None:
            raise ValueError('nonfinite baseline point')
        origin = anchor['parameters']
        entry = dict(uniform_baseline=anchor, points=[], best=None, skipped_invalid_points=0)
        report['attacks'][attack] = entry
        strategies = ['balanced', 'prefix'] + [f'saturated-{count}-{order}'
            for count in args.saturated_blocks for order in ('prefix', 'short-first')]
        for strategy in strategies:
            for bz, zz, dz in itertools.product(args.beta_offsets, args.zeta_offsets, args.d_offsets):
                beta, zeta, d = origin['beta'] + bz, origin['zeta'] + zz, origin['d'] + dz
                if not (40 <= beta <= args.max_beta and 0 <= zeta < n and d >= max(beta, n - zeta + 1)):
                    entry['skipped_invalid_points'] += 1
                    continue
                if strategy.startswith('saturated-'):
                    _, count, order = strategy.split('-', 2)
                    try:
                        guessed = saturated_layout(widths, zeta, int(count), order == 'short-first')
                    except ValueError:
                        entry['skipped_invalid_points'] += 1
                        continue
                else:
                    guessed = layout(widths, zeta, strategy)
                point = dict(layout=strategy, beta=beta, zeta=zeta, d=d, best=None, evaluated_strata=0, zero_residual_strata=0)
                for weight, count, probability, norm, variances in exact_strata(widths, guessed):
                    if not norm:
                        point['zero_residual_strata'] += 1
                        continue
                    variance = norm / (n - zeta)
                    scalar = QQ(variance.numerator) / variance.denominator
                    secret = ND.NoiseDistribution(n=n, mean=0, stddev=sqrt(scalar), bounds=(-1, 1), _density=scalar)
                    params = LWE.Parameters(n=n, q=2**q, Xs=secret, Xe=ND.DiscreteGaussian(sigma))
                    result = PrimalHybrid.cost(beta=beta, params=params, zeta=zeta, d=d,
                        babai=mitm, mitm=mitm, red_cost_model=RC.BDGL16,
                        search_space=QQ(count), hit_probability=QQ(probability.numerator) / probability.denominator)
                    record = cost_record(result, log)
                    point['evaluated_strata'] += 1
                    if record['bits'] is not None and (point['best'] is None or record['bits'] < point['best']['bits']):
                        point['best'] = dict(record, guessed_weight=weight, residual_norm_squared=str(norm),
                            residual_average_variance=str(variance), candidate_count=str(count),
                            hit_probability_numerator=str(probability.numerator), hit_probability_denominator=str(probability.denominator),
                            remaining_coordinate_variance_per_block=[str(v) if v is not None else None for v in variances])
                entry['points'].append(point)
                if point['best'] is not None and (entry['best'] is None or point['best']['bits'] < entry['best']['bits']):
                    entry['best'] = dict(point['best'], layout=strategy, grid_beta=beta, grid_zeta=zeta, grid_d=d, guessed_per_block=guessed)
                print(attack, strategy, beta, zeta, d, 'bits', None if point['best'] is None else point['best']['bits'], flush=True)
                save()
                # Bound memoized Sage GSA vectors; reuse within each weight sweep.
                PrimalHybrid.cost.clear_cache()
                PrimalHybrid.beta_params.clear_cache()
    report['bounded_grid_complete'] = True
    save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--estimator', type=Path, default=Path(__file__).resolve().parents[1] / 'lattice-estimator')
    parser.add_argument('--baseline-report', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--max-beta', type=int, default=4096)
    parser.add_argument('--beta-offsets', type=offsets, default=offsets('-8,0,8'))
    parser.add_argument('--zeta-offsets', type=offsets, default=offsets('-128,0,128'))
    parser.add_argument('--d-offsets', type=offsets, default=offsets('-128,0,128'))
    parser.add_argument('--attacks', nargs='+', choices=('bdd_hybrid','bdd_mitm_hybrid'), default=['bdd_mitm_hybrid'])
    parser.add_argument('--saturated-blocks', nargs='+', type=int, default=[],
                        help='also test fully guessed block counts with prefix and shortest-block-first orders')
    args = parser.parse_args()
    if args.max_beta < 40:
        parser.error('--max-beta must be at least 40')
    if any(count < 1 for count in args.saturated_blocks) or len(set(args.saturated_blocks)) != len(args.saturated_blocks):
        parser.error('--saturated-blocks must contain distinct positive counts')
    self_test()
    if args.self_test:
        return 0
    if args.baseline_report is None or args.output is None:
        parser.error('--baseline-report and --output are required')
    run(args)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
