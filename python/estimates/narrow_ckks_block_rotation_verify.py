#!/usr/bin/env python3
"""Verify actual-block rotation counts with bit masks and grouped binomials.

This independent oracle does not use the sampler's NumPy correlations or
linear-factor dynamic program. It verifies exact event masses and labelled
signed projection counts, not decoding probabilities or security.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from math import comb, prod
from pathlib import Path


def polynomial(groups, maximum, signed):
    result = [1] + [0] * maximum
    for (width, selected), multiplicity in groups:
        a = int(selected < width) if signed else width - selected
        b = 2 * selected if signed else selected
        factor = [comb(multiplicity, k) * a**(multiplicity-k) * b**k
                  for k in range(min(multiplicity, maximum) + 1)]
        result = [sum(result[j] * factor[k-j]
                      for j in range(max(0, k-len(factor)+1), k+1))
                  for k in range(maximum+1)]
    return result


def rotation_counts(widths, guessed, maximum):
    mask = 0
    start = 0
    windows = []
    for width, selected in zip(widths, guessed):
        if not 0 <= selected <= width or width < 1:
            raise ValueError('invalid block geometry')
        mask |= ((1 << selected)-1) << start
        windows.append((start, (1 << width)-1, width))
        start += width
    cache = {}
    for _ in range(start):
        groups = tuple(sorted(Counter(
            (width, ((mask >> offset) & window).bit_count())
            for offset, window, width in windows).items()))
        if groups not in cache:
            cache[groups] = (polynomial(groups, maximum, False),
                             polynomial(groups, maximum, True))
        yield cache[groups]
        mask = (mask >> 1) | ((mask & 1) << (start-1))


def self_test():
    widths = [2, 3, 1]
    blocks = []
    start = 0
    for width in widths:
        blocks.append(list(range(start, start+width)))
        start += width
    cases = 0
    for guessed in itertools.product(*(range(w+1) for w in widths)):
        selected = {p for block, count in zip(blocks, guessed) for p in block[:count]}
        for r, (masses, supports) in enumerate(rotation_counts(widths, guessed, len(widths))):
            accepted = [[] for _ in range(len(widths)+1)]
            for secret in itertools.product(*blocks):
                projection = tuple((p+r) % start for p in secret if (p+r) % start in selected)
                accepted[len(projection)].append(projection)
            for w, projections in enumerate(accepted):
                signed = {tuple(zip(p, signs)) for p in projections
                          for signs in itertools.product((-1, 1), repeat=w)}
                assert masses[w] == len(projections)
                assert supports[w] == len(signed)
                cases += 1
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reports = [json.loads(p.read_text()) for p in args.report]
    first = reports[0]
    widths, guessed = first['widths'], first['guessed_per_block']
    weights = sorted({d['guessed_weight'] for d in reports})
    if len(widths) != len(guessed) or not all(type(w) is int and 0 <= w <= len(widths) for w in weights):
        raise ValueError('invalid widths or guessed weights')
    mass_sum = [0] * (max(weights)+1)
    support_sum = [0] * len(mass_sum)
    for masses, supports in rotation_counts(widths, guessed, max(weights)):
        for w in weights:
            mass_sum[w] += masses[w]
            support_sum[w] += supports[w]
    for d in reports:
        if not d['complete'] or d['n'] != sum(widths) or d['h'] != len(widths) or d['widths'] != widths or d['guessed_per_block'] != guessed:
            raise ValueError('report geometry mismatch')
        w = d['guessed_weight']
        if (int(d['rotation_event_mass_sum']) != mass_sum[w] or
                int(d['rotation_labelled_candidate_count']) != support_sum[w] or
                int(d['position_vector_denominator']) != prod(widths) or
                Fraction(d['expected_matching_rotations']) != Fraction(mass_sum[w], prod(widths))):
            raise ValueError('exact rotation count mismatch')
    result = dict(scope=__doc__, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        exhaustive_toy_cases=self_test(), n=sum(widths), weights=weights, rotations_checked=sum(widths),
        exact_mass_sums={str(w): str(mass_sum[w]) for w in weights},
        exact_signed_labelled_support_sums={str(w): str(support_sum[w]) for w in weights},
        reports={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in args.report},
        passed=True, security_validated=False)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Grouped-binomial and cyclic bit-mask verification passed:', weights, 'toy cases:', result['exhaustive_toy_cases'])


if __name__ == '__main__':
    main()
