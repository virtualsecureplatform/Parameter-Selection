#!/usr/bin/env python3
"""Actual public-block rotation events: exact marginals and posterior sampling.

For K matching rotations, Pr[K>0] = sum_r Pr[event_r] E[1/K | event_r].
We sample r proportional to its exact event mass, then draw the secret exactly
conditional on that rotation's guessed weight. This avoids sampling a rare
unconditional event and does not assume independence between rotations.
Confidence intervals are probability diagnostics, not security certificates.
"""
import argparse
from bisect import bisect_right
from collections import Counter
import datetime
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import statistics


def coefficients(widths, guessed, weight):
    mass = [1] + [0] * weight
    support = [1] + [0] * weight
    for width, selected in zip(widths, guessed):
        if not 0 <= selected <= width or width < 1:
            raise ValueError('invalid block geometry')
        for k in range(weight, -1, -1):
            mass[k] = (width-selected)*mass[k] + (selected*mass[k-1] if k else 0)
            support[k] = int(selected<width)*support[k] + (2*selected*support[k-1] if k else 0)
    return mass[weight], support[weight]


def suffix_masses(widths, guessed, weight):
    suffix = [[0]*(weight+1) for _ in range(len(widths)+1)]
    suffix[-1][0] = 1
    for j in range(len(widths)-1, -1, -1):
        for k in range(weight+1):
            suffix[j][k] = (widths[j]-guessed[j])*suffix[j+1][k] + (guessed[j]*suffix[j+1][k-1] if k else 0)
    return suffix


def self_test():
    widths = [2,3,1]
    n = sum(widths)
    blocks = [list(range(sum(widths[:j]),sum(widths[:j+1]))) for j in range(len(widths))]
    secrets = list(itertools.product(*blocks))
    cases = 0
    for mask in itertools.product((0,1),repeat=n):
        for weight in range(len(widths)+1):
            matching = [[sum(mask[(p+r)%n] for p in secret)==weight for r in range(n)] for secret in secrets]
            union = Fraction(sum(any(events) for events in matching),len(secrets))
            identity = Fraction(0)
            decode_probability = Fraction(1,2)
            joint_union = sum((1-(1-decode_probability)**sum(events) for events in matching),Fraction(0))/len(secrets)
            joint_identity = Fraction(0)
            for r in range(n):
                guessed = [sum(mask[(p+r)%n] for p in block) for block in blocks]
                mass,support = coefficients(widths,guessed,weight)
                accepted = [s for s,events in zip(secrets,matching) if events[r]]
                assert mass == len(accepted)
                suffix = suffix_masses(widths,guessed,weight)
                assert suffix[0][weight] == mass
                projections = set()
                for secret in accepted:
                    probability = Fraction(1)
                    left = weight
                    for j,p in enumerate(secret):
                        yes = bool(mask[(p+r)%n])
                        choices = guessed[j] if yes else widths[j]-guessed[j]
                        next_mass = suffix[j+1][left-int(yes)]
                        probability *= Fraction(choices*next_mass,suffix[j][left])*Fraction(1,choices)
                        left -= int(yes)
                    assert left==0 and probability==Fraction(1,mass)
                    for signs in itertools.product((-1,1),repeat=len(widths)):
                        projected = [0]*n
                        for p,sign in zip(secret,signs):
                            q=(p+r)%n
                            if mask[q]: projected[q]=sign if p+r<n else -sign
                        projections.add(tuple(projected))
                assert support==len(projections)
                if mass:
                    inverse_mean=sum((Fraction(1,sum(events)) for events in matching if events[r]),Fraction(0))/mass
                    identity += Fraction(mass,len(secrets))*inverse_mean
                    joint_mean = sum(((1-(1-decode_probability)**sum(events))/sum(events)
                        for events in matching if events[r]),Fraction(0))/mass
                    joint_identity += Fraction(mass,len(secrets))*joint_mean
            assert identity==union
            assert joint_identity==joint_union
            cases += 1
    print('exact rotation marginals, signed projections, posterior paths and secret/joint-decoding union identities passed:',cases,'cases',flush=True)


def run(args):
    import numpy as np
    profile=json.loads(args.profile_report.read_text())
    best=profile['attacks']['bdd_mitm_hybrid']['best']
    widths=profile['block_widths']
    guessed=best['guessed_per_block']
    n=sum(widths);h=len(widths);weight=best['guessed_weight'] if args.weight is None else args.weight
    if n!=profile['n'] or not 0<=weight<=h or sum(guessed)!=best['grid_zeta']:
        raise ValueError('invalid profile geometry')
    starts=np.array([sum(widths[:j]) for j in range(h)],dtype=np.int64)
    mask=np.zeros(n,dtype=np.int64)
    for start,count in zip(starts,guessed): mask[start:start+count]=1
    prefix=np.concatenate(([0],np.cumsum(np.tile(mask,2))))
    indices=(np.arange(n,dtype=np.int64)[:,None]+starts[None,:])%n
    all_guessed=(prefix[indices+np.array(widths)[None,:]]-prefix[indices]).tolist()
    masses=[];supports=[]
    for geometry in all_guessed:
        mass,support=coefficients(widths,geometry,weight)
        masses.append(mass);supports.append(support)
    denominator=math.prod(widths)
    if args.weight is None:
        assert supports[0]==int(best['candidate_count'])
        assert Fraction(masses[0],denominator)==Fraction(int(best['hit_probability_numerator']),int(best['hit_probability_denominator']))
    total_mass=sum(masses)
    if not total_mass:raise ValueError('no supported rotation event')
    cumulative=[];total=0
    for mass in masses:total+=mass;cumulative.append(total)
    expected_matches=Fraction(total_mass,denominator)
    probabilities=[mass/denominator for mass in masses]
    independent=1.0 if any(p==1 for p in probabilities) else -math.expm1(sum(math.log1p(-p) for p in probabilities))
    report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=__doc__,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        profile_report=str(args.profile_report),profile_report_sha256=hashlib.sha256(args.profile_report.read_bytes()).hexdigest(),
        n=n,h=h,widths=widths,guessed_per_block=guessed,guessed_weight=weight,layout=best['layout'],
        seed=args.seed,samples_requested=args.samples,numpy_version=np.__version__,exact_self_test_passed=True,
        zero_rotation_control_matches_profile=args.weight is None,
        rotation_event_mass_sum=str(total_mass),position_vector_denominator=str(denominator),
        expected_matching_rotations=str(expected_matches),rotation_labelled_candidate_count=str(sum(supports)),
        minimum_single_rotation_probability=min(probabilities),maximum_single_rotation_probability=max(probabilities),
        independent_rotation_probability=independent,union_probability_upper_bound=min(1.0,float(expected_matches)),
        complete=False,security_validated=False,
        limitations=['Exact-weight event at one public guessing mask; no attack optimization',
            'Rotation-labelled guessing count is an upper bound; feasible isometric MITM dictionaries not constructed',
            'Monte Carlo intervals quantify this event model, not security or lattice geometry',
            'Ring lattice cost and correlated evaluation keys are not modeled'])
    def save():args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save()
    rng=random.Random(args.seed);mask_spectrum=np.fft.rfft(mask);inverse_counts=[];matching_counts=[];max_fft_error=0
    for trial in range(args.samples):
        r=bisect_right(cumulative,rng.randrange(total_mass));geometry=all_guessed[r]
        suffix=suffix_masses(widths,geometry,weight);left=weight;positions=[]
        for j,(start,width,selected) in enumerate(zip(starts,widths,geometry)):
            yes_mass=selected*suffix[j+1][left-1] if left else 0
            no_mass=(width-selected)*suffix[j+1][left]
            yes=rng.randrange(yes_mass+no_mass)<yes_mass
            choices=np.flatnonzero(mask[(int(start)+np.arange(width)+r)%n]==int(yes))
            if len(choices)!=(selected if yes else width-selected):raise ValueError('conditional coordinate count mismatch')
            positions.append(int(start)+int(choices[rng.randrange(len(choices))]));left-=int(yes)
        if left:raise ValueError('conditional weight mismatch')
        secret=np.zeros(n,dtype=np.int64);secret[positions]=1
        raw=np.fft.irfft(np.conj(np.fft.rfft(secret))*mask_spectrum,n=n)
        rounded=np.rint(raw).astype(np.int64);error=float(np.max(np.abs(raw-rounded)));max_fft_error=max(max_fft_error,error)
        if error>1e-6 or rounded.min()<0 or rounded.max()>h:raise ValueError('FFT integer correlation check failed')
        for check in (0,r,n-1):
            if int(rounded[check])!=sum(int(mask[(p+check)%n]) for p in positions):raise ValueError('direct rotation correlation mismatch')
        k=int(np.count_nonzero(rounded==weight))
        if k<1 or rounded[r]!=weight:raise ValueError('conditioned rotation absent')
        inverse_counts.append(1/k);matching_counts.append(k)
        if (trial+1)%500==0:print('posterior samples',trial+1,flush=True)
    inverse_mean=statistics.mean(inverse_counts)
    se=statistics.stdev(inverse_counts)/math.sqrt(args.samples)
    scale=float(expected_matches);estimate=scale*inverse_mean
    hoeffding=scale*math.sqrt(math.log(40)/(2*args.samples))
    report.update(complete=True,samples_completed=args.samples,mean_inverse_matching_rotations=inverse_mean,
        matching_rotations_min=min(matching_counts),matching_rotations_median=statistics.median(matching_counts),matching_rotations_max=max(matching_counts),
        matching_rotation_histogram={str(k):count for k,count in sorted(Counter(matching_counts).items())},
        union_probability_estimate=estimate,mc_standard_error=scale*se,
        approximate_mc95_interval=[max(0,estimate-1.96*scale*se),min(1,estimate+1.96*scale*se)],
        iid_hoeffding95_interval=[max(0,estimate-hoeffding),min(1,estimate+hoeffding)],
        independence_overestimate_ratio=independent/estimate,max_fft_integer_residual=max_fft_error)
    save();print('union probability',estimate,'independence overestimate',independent/estimate,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--profile-report',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--samples',type=int,default=4000)
    parser.add_argument('--seed',type=int,default=2119)
    parser.add_argument('--weight',type=int)
    args=parser.parse_args();self_test()
    if args.self_test:return
    if args.profile_report is None or args.output is None or args.samples<2:
        parser.error('--profile-report, --output and at least 2 samples required')
    run(args)


if __name__=='__main__':main()
