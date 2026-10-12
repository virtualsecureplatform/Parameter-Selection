#!/usr/bin/env python3
"""Pinned-author ring-rotation hybrid fixed-point sensitivity, uniform proxy only.

This does not model the port's public one-per-block support distribution or
certify security. It evaluates the author's rotation-amplification formula with
exact uniform-support counts and 512-bit probabilities to avoid cancellation.
"""
import argparse
import datetime
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import sys
from narrow_ckks_block_hybrid import cost_record, estimator_revision, hash_file


def counts(n, h, zeta, maximum):
    low = max(0, h - (n - zeta))
    support = sum(math.comb(zeta, k) * 2**k for k in range(low, maximum + 1))
    probability = Fraction(sum(math.comb(zeta, k) * math.comb(n - zeta, h - k)
                              for k in range(low, maximum + 1)), math.comb(n, h))
    return support, probability


def self_test():
    cases = 0
    for n, h in ((4, 2), (5, 3), (3, 1)):
        secrets = []
        for positions in itertools.combinations(range(n), h):
            for signs in itertools.product((-1, 1), repeat=h):
                vector = [0] * n
                for position, sign in zip(positions, signs):
                    vector[position] = sign
                secrets.append(tuple(vector))
        for zeta in range(n + 1):
            for maximum in range(max(0, h - (n - zeta)), min(h, zeta) + 1):
                accepted = [s for s in secrets if sum(x != 0 for x in s[:zeta]) <= maximum]
                support, probability = counts(n, h, zeta, maximum)
                assert support == len({s[:zeta] for s in accepted})
                assert probability == Fraction(len(accepted), len(secrets))
                cases += 1
    print('uniform signed-support counts/probabilities passed exhaustive cases:', cases, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--implementation', type=Path)
    parser.add_argument('--baseline-report', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--probability-bits', type=int, default=512)
    args = parser.parse_args()
    self_test()
    if args.self_test:
        return
    if any(value is None for value in (args.implementation, args.baseline_report, args.output)):
        parser.error('--implementation, --baseline-report and --output are required')
    if not 256 <= args.probability_bits <= 2048:
        parser.error('--probability-bits must be between 256 and 2048')
    implementation = args.implementation.resolve()
    sys.path.insert(0, str(implementation))
    from lattice_estimator.estimator import LWE, ND, RC
    from lattice_estimator.estimator import conf, lwe_primal
    from lwe_rot_primal import RotPrimalHybrid
    from sage.all import RealField, log
    conf.max_beta = lwe_primal.max_beta_global = 4096
    # Rotated optimizers have their own imported global; this probe uses fixed
    # parameters, but retain the same limit explicitly for provenance.
    import lwe_rot_primal
    lwe_rot_primal.max_beta_global = 4096
    saved = json.loads(args.baseline_report.read_text())
    n, h, q, sigma = (saved[k] for k in ('n', 'h', 'log_q', 'sigma'))
    anchor = saved['attacks']['bdd_mitm_hybrid']['uniform_baseline']
    point = anchor['parameters']
    params = LWE.Parameters(n=n, q=2**q, Xs=ND.SparseTernary((h+1)//2, h//2, n=n), Xe=ND.DiscreteGaussian(sigma))
    real = RealField(args.probability_bits)
    report = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), scope=__doc__,
        source_sha256=hash_file(__file__), author_repository='https://github.com/TabOg/mlwe-hybrids',
        author_revision=estimator_revision(implementation.parent),
        author_source_sha256=hash_file(implementation/'lwe_rot_primal.py'),
        estimator_revision=estimator_revision(implementation/'lattice_estimator'),
        paper='https://eprint.iacr.org/2026/279', baseline_report_sha256=hash_file(args.baseline_report),
        n=n,h=h,log_q=q,sigma=sigma,point=point,probability_precision_bits=args.probability_bits,
        exact_count_self_test_passed=True,
        models={},fixed_point_complete=False,complete_attack_optimization=False,security_validated=False,
        limitations=['Uniform-support proxy; actual public blocks are not rotation invariant',
            'Author rotation-event independence and square-root MITM heuristics retained',
            'Only saved beta/zeta/d point evaluated, no ring-parameter optimization',
            'Correlated evaluation keys and other ring-specific attacks unqualified'])
    def save():
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save()
    for cost_name in ('BDGL16','MATZOV'):
        entries={}
        report['models'][cost_name]=entries
        for name,degree,heuristic in (('lwe_control',1,'square root'),('ring_estimator_mitm',n,'estimator'),('ring_square_root_mitm',n,'square root')):
            entry=dict(best=None,nonfinite_weights=[],weights_evaluated=0)
            entries[name]=entry
            for weight in range(max(0,h-(n-point['zeta'])), min(h,point['zeta'])+1):
                support,probability=counts(n,h,point['zeta'],weight)
                p=real(probability.numerator)/probability.denominator
                amplified=1-(1-p)**degree
                if amplified<=0 or amplified>1:
                    raise ValueError('rotation probability lost precision or range')
                result=RotPrimalHybrid.cost(beta=point['beta'],params=params,zeta=point['zeta'],d=point['d'],
                    babai=True,mitm=True,poly_degree=degree,mitm_heuristic=heuristic,
                    red_cost_model=getattr(RC,cost_name),search_space=support*degree,hit_probability=amplified)
                record=cost_record(result,log);entry['weights_evaluated']+=1
                if record['bits'] is None:
                    entry['nonfinite_weights'].append(dict(weight=weight,raw_cost=record['cost']))
                elif entry['best'] is None or record['bits']<entry['best']['bits']:
                    entry['best']=dict(record,maximum_guessed_weight=weight,
                        single_rotation_support=str(support),single_rotation_probability=str(probability),
                        amplified_probability=str(amplified),poly_degree=degree,mitm_heuristic=heuristic)
            print(cost_name,name,None if entry['best'] is None else entry['best']['bits'],flush=True)
            save()
    report['fixed_point_complete']=True
    save()


if __name__=='__main__':
    main()
