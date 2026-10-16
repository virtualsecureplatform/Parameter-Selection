#!/usr/bin/env python3
"""Structured ring joint-decoding sensitivity from posterior rotation histograms.

Compares common decoding outcomes, conditional decoding independence, and an
optimistic union bound. None establishes actual correlations or security.
The full square-root speedup on labelled dictionaries remains a heuristic.
"""
import argparse
import datetime
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
from narrow_ckks_block_hybrid import cost_record,estimator_revision
from narrow_ckks_block_hybrid_grid import offsets
from narrow_ckks_security import configure_max_beta


def geometric_sum(x,k):
    """Return sum(x**j,j=0..k-1) without subtracting nearly equal numbers."""
    if k==0:return x*0
    power,total=geometric_pair(x,k)
    return total


def geometric_pair(x,k):
    if k==0:return x*0+1,x*0
    power,total=geometric_pair(x,k//2)
    total=total*(1+power);power=power*power
    if k%2:return power*x,total+power
    return power,total


def joint_factor(histogram,p):
    size=sum(histogram.values())
    return sum((count*p*geometric_sum(1-p,k)/k for k,count in histogram.items()),p*0)/size



def normalize_joint_success(common, independent, upper, relative_tolerance=0):
    """Correct only arithmetic roundoff against proven probability bounds."""
    if not 0<=common<=upper<=1:
        raise ValueError('invalid joint probability bounds')
    tolerance=relative_tolerance*max(common,independent,upper)
    if common-independent>tolerance or independent-upper>tolerance:
        raise ValueError('joint decoder probability bounds failed')
    adjusted=independent<common or independent>upper
    return min(upper,max(common,independent)),adjusted


def stable_probability_cost(geometry, support, joint, beta, zeta):
    """Apply the pinned hybrid formula to an actual joint probability.

    Also supports rotation unions whose success exceeds one decoder marginal;
    their joint/marginal ratio is not a valid single-secret hit probability.
    Keep the estimator's BKZ/SVP geometry and its 53-bit square-root search
    factor. Only success multiplication/amplification use MPFR and log1p.
    This does not strengthen the underlying attack model.
    """
    from estimator.cost import Cost
    from sage.all import RR, RealField, ceil, sqrt
    real=RealField(512)
    joint=real(joint)
    if not 0<joint<=1:
        raise ValueError('stable amplification requires a probability in (0,1]')
    target=real(0.99)  # Match the pinned backend's binary-float target.
    trials=1 if joint>=target else ceil((-target).log1p()/(-joint).log1p())
    svp=geometry['svp_cost'].repeat(RR(sqrt(support)))
    red=geometry['bkz_cost']['rop']
    result=Cost(rop=red+svp['rop'],red=red,svp=svp['rop'],beta=beta,
                eta=geometry['eta'],zeta=zeta,d=geometry['d'],prob=joint)
    result['|S|']=support
    result.register_impermanent({'|S|':False},rop=True,red=True,svp=True,
                                eta=False,zeta=False,prob=False)
    return result.repeat(trials)


def self_test():
    histogram={1:3,2:2,5:1}
    for p in (Fraction(0),Fraction(1,4),Fraction(1,2),Fraction(1)):
        direct=Fraction(0)
        for k,count in histogram.items():
            success=Fraction(0)
            for outcomes in itertools.product((0,1),repeat=k):
                if any(outcomes):success+=p**sum(outcomes)*(1-p)**(k-sum(outcomes))
            direct+=count*success/k
        direct/=sum(histogram.values())
        assert joint_factor(histogram,p)==direct
        inverse_mean=sum((Fraction(count,k) for k,count in histogram.items()),Fraction(0))/sum(histogram.values())
        assert p*inverse_mean<=direct<=min(inverse_mean,p)
    # Summation by parts: a uniform empirical-CDF error bounds every
    # decreasing decoder function at once, including a data-selected p.
    left=[Fraction(x,16) for x in (1,3,2,4,1,2,1,2)]
    right=[Fraction(x,16) for x in (2,1,4,1,3,1,2,2)]
    cdf=[];running=Fraction(0)
    for a,b in zip(left,right):
        running+=a-b;cdf.append(running)
    assert running==0
    epsilon=max(map(abs,cdf))
    for p in (Fraction(0),Fraction(1,100),Fraction(1,4),Fraction(1,2),Fraction(1)):
        values=[p*geometric_sum(1-p,k)/k for k in range(1,9)]
        assert all(a>=b for a,b in zip(values,values[1:]))
        error=sum((a-b)*v for a,b,v in zip(left,right,values))
        by_parts=sum(c*(a-b) for c,a,b in zip(cdf,values,values[1:]))
        assert error==by_parts and abs(error)<=epsilon*p
    values=[Fraction(1,k) for k in range(1,9)]
    assert abs(sum((a-b)*v for a,b,v in zip(left,right,values)))<=epsilon
    # A rotation union can legitimately succeed more often than one decoder.
    # Two matching rotations, with secret-event mass two, have union mass one.
    p=Fraction(1,4)
    joint=2*joint_factor({2:100},p)
    assert joint==Fraction(7,16) and joint/p>1 and joint<=1
    tiny=Fraction(1,2**500)
    assert normalize_joint_success(Fraction(1,2),Fraction(3,4)+tiny,Fraction(3,4),Fraction(1,2**450))==(Fraction(3,4),True)
    try:normalize_joint_success(Fraction(1,2),Fraction(4,5),Fraction(3,4),Fraction(1,2**450))
    except ValueError:pass
    else:raise AssertionError('material probability-bound violation was hidden')
    print('joint decoding formula, stable geometric sum and probability bounds passed exact Bernoulli enumeration',flush=True)
    print('monotone decoder functions and uniform CDF expectation bounds passed exact summation by parts',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--estimator',type=Path)
    parser.add_argument('--profile-report',type=Path)
    parser.add_argument('--probability-report',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--beta-offsets',type=offsets,default=[0])
    parser.add_argument('--d-offsets',type=offsets,default=[0])
    parser.add_argument('--log-q',type=int,
                        help='candidate sparse-key modulus bits; default is the saved profile')
    parser.add_argument('--noise-stddev',type=float,
                        help='integer-noise RMS for the Gaussian proxy; default is the saved profile')
    args=parser.parse_args();self_test()
    if args.self_test:return
    if any(value is None for value in (args.estimator,args.profile_report,args.probability_report,args.output)):
        parser.error('--estimator, --profile-report, --probability-report and --output are required')
    if args.noise_stddev is not None and (not math.isfinite(args.noise_stddev) or args.noise_stddev<=0):
        parser.error('--noise-stddev must be finite and positive')
    sys.path.insert(0,str(args.estimator.resolve()))
    from estimator import LWE,ND,RC
    from estimator.lwe_primal import PrimalHybrid
    from sage.all import QQ,RealField,sqrt,log,oo
    configure_max_beta(4096);real=RealField(512)
    profile=json.loads(args.profile_report.read_text());prob=json.loads(args.probability_report.read_text())
    log_q=profile['log_q'] if args.log_q is None else args.log_q
    sigma=profile['sigma'] if args.noise_stddev is None else args.noise_stddev
    if not math.isfinite(sigma) or sigma<=0:
        raise ValueError('noise standard deviation must be finite and positive')
    if type(log_q) is not int or not 1<=log_q<=4096:
        raise ValueError('sparse-key modulus bits must be between 1 and 4096')
    b=profile['attacks']['bdd_mitm_hybrid']['best'];n=profile['n'];zeta=b['grid_zeta']
    # The mask fixes zeta and the layout, but the posterior sampler can
    # condition on a different exact guessed-weight stratum of that mask.
    # Its marginals and labelled support count already use this weight.
    weight=prob['guessed_weight']
    if type(weight) is not int or not 0<=weight<=profile['h']:
        raise ValueError('invalid posterior guessed weight')
    if (not prob['complete'] or prob['n']!=n or prob['widths']!=profile['block_widths'] or
        prob['h']!=profile['h'] or prob['layout']!=b['layout'] or
        prob['guessed_per_block']!=b['guessed_per_block'] or
        prob['profile_report_sha256']!=hashlib.sha256(args.profile_report.read_bytes()).hexdigest()):
        raise ValueError('posterior histogram/profile mismatch')
    histogram={int(k):int(count) for k,count in prob['matching_rotation_histogram'].items()}
    if not histogram or min(histogram)<1 or max(histogram)>n or min(histogram.values())<1 or sum(histogram.values())!=prob['samples_completed']:
        raise ValueError('invalid posterior matching histogram')
    inverse_mean=sum((real(count)/k for k,count in histogram.items()),real(0))/sum(histogram.values())
    event_sum=real(QQ(prob['expected_matching_rotations']))
    empirical_union=event_sum*inverse_mean
    union=min(real(1),empirical_union)
    # One-sample DKW-Massart bounds the entire posterior K CDF with 95%
    # coverage. Monotone h_p(K)=p*sum((1-p)**j)/K has variation <=p;
    # 1/K has variation <=1. These bounds hold simultaneously for all p,
    # hence survive beta/d selection on this same fixed-layout histogram.
    iid_radius=(real(40).log()/(2*sum(histogram.values()))).sqrt()
    union_iid_upper=min(real(1),union+event_sum*iid_radius)
    if abs(float(empirical_union)/prob['union_probability_estimate']-1)>1e-12:
        raise ValueError('histogram does not reproduce secret union estimate')
    variance=QQ(profile['h']-weight)/(n-zeta)
    if variance<=0:raise ValueError('zero residual requires separate exhaustive model')
    params=LWE.Parameters(n=n,q=2**log_q,Xs=ND.NoiseDistribution(n=n,mean=0,stddev=sqrt(variance),bounds=(-1,1),_density=variance),Xe=ND.DiscreteGaussian(sigma))
    support=QQ(prob['rotation_labelled_candidate_count'])
    report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=__doc__,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),estimator_revision=estimator_revision(args.estimator),
        profile_sha256=hashlib.sha256(args.profile_report.read_bytes()).hexdigest(),probability_report=str(args.probability_report),
        probability_report_sha256=hashlib.sha256(args.probability_report.read_bytes()).hexdigest(),probability_seed=prob['seed'],
        n=n,log_q=log_q,profile_log_q=profile['log_q'],candidate_modulus=log_q!=profile['log_q'],
        sigma=sigma,profile_sigma=profile['sigma'],candidate_noise=sigma!=profile['sigma'],
        noise_proxy='DiscreteGaussian with supplied integer RMS; exact sampler distribution is not modeled',
        probability_reuse_scope='Posterior matching histogram depends on the unchanged secret/guessing geometry; decoder probabilities are recomputed at the selected modulus and noise',
        empirical_secret_union=str(empirical_union),central_secret_union=str(union),
        central_union_clipped_to_one=empirical_union>1,
        relative_probability_roundoff_tolerance=str(real(2)**-450),
        zeta=zeta,weight=weight,layout=b['layout'],histogram_samples=sum(histogram.values()),beta_offsets=args.beta_offsets,d_offsets=args.d_offsets,
        scenarios={},complete=False,security_validated=False,complete_attack_optimization=False,
        posterior_uncertainty=dict(method='one-sample DKW-Massart empirical-CDF bound',
            confidence=0.95,cdf_radius=str(iid_radius),
            theorem_source='https://arxiv.org/abs/1107.5356',
            formula='epsilon=sqrt(log(2/0.05)/(2*samples))',
            simultaneous_scope='All decoder marginals and beta/d choices using this same fixed-layout posterior histogram',
            expectation_bound='epsilon for 1/K; p*epsilon for h_p(K)=(1-(1-p)^K)/K by monotone summation by parts',
            excludes='Other layouts, weights, zeta, seeds, decoding-correlation hypotheses, and security certification'),
        limitations=['Decoder marginal uses scalar GSA Babai and MITM probabilities; actual cross-rotation correlations unmodeled',
            'Conditional decoder independence and optimistic union bound are sensitivities, not established attacks',
            'Only one zeta/layout/weight and supplied beta/d offsets; no global optimum',
            'Full square-root work on labelled dictionaries assumed; compatible decompositions not constructed',
            'Posterior finite-sample uncertainty is not a security confidence interval',
            'Candidate modulus overrides do not establish matching C++ accuracy or exposed-key security',
            'Correlated evaluation keys remain unqualified'])
    for cost_name in ('BDGL16','MATZOV'):
        entry=dict(points=[],best={});report['scenarios'][cost_name]=entry
        for bo,do in itertools.product(args.beta_offsets,args.d_offsets):
            beta=b['grid_beta']+bo;d=b['grid_d']+do
            if not 40<=beta<=4096 or d<max(beta,n-zeta+1):continue
            geometry=PrimalHybrid.beta_params(beta=beta,params=params,zeta=zeta,babai=True,mitm=True,m=oo,d=d,red_cost_model=getattr(RC,cost_name))
            if 'babai_probability' not in geometry:
                entry['points'].append(dict(beta=beta,d=d,finite=False,reason='decoder geometry unavailable'));continue
            p=real(geometry['babai_probability'])*real(geometry['mitm_probability'])
            if p<=0 or p>1:
                entry['points'].append(dict(beta=beta,d=d,finite=False,reason='zero/out-of-range decoder marginal'));continue
            independent=min(real(1),event_sum*joint_factor(histogram,p))
            common=union*p;upper=min(union,event_sum*p)
            independent,roundoff_adjusted=normalize_joint_success(common,independent,upper,real(2)**-450)
            common_iid_upper=union_iid_upper*p
            independent_iid_upper=min(real(1),event_sum*p,independent+event_sum*p*iid_radius)
            optimistic_iid_upper=min(union_iid_upper,event_sum*p)
            point=dict(beta=beta,d=d,probability_roundoff_adjusted=roundoff_adjusted,decode_marginal=str(p),joint_success_common=str(common),joint_success_independent=str(independent),joint_success_upper_bound=str(upper),costs={})
            for name,joint in (('common_decoder',common),('independent_decoder',independent),('optimistic_joint_bound',upper),
                               ('common_decoder_iid95',common_iid_upper),('independent_decoder_iid95',independent_iid_upper),
                               ('optimistic_joint_bound_iid95',optimistic_iid_upper)):
                effective=joint/p
                if not 0<joint<=1:raise ValueError('joint success probability out of range')
                if effective>1:
                    result=stable_probability_cost(geometry,support,joint,beta,zeta)
                    record=cost_record(result,log)
                    record['joint_probability_method']='pinned geometry/cost formula with rotation-union success; joint/marginal exceeds one, so it is not passed as a single-secret hit probability'
                else:
                    result=PrimalHybrid.cost(beta=beta,params=params,zeta=zeta,d=d,babai=True,mitm=True,
                        red_cost_model=getattr(RC,cost_name),search_space=support,hit_probability=effective)
                    record=cost_record(result,log)
                backend_probability=effective*geometry['babai_probability']*geometry['mitm_probability']
                if record['bits'] is None and not backend_probability:
                    original=record
                    result=stable_probability_cost(geometry,support,joint,beta,zeta)
                    record=cost_record(result,log)
                    record['backend_result']=original
                    record['backend_probability']=str(backend_probability)
                    record['backend_probability_parent']=str(backend_probability.parent())
                    record['numerical_method']='same pinned geometry/cost formula; 512-bit joint probability and log1p amplification'
                point['costs'][name]=record
                if record['bits'] is not None and (name not in entry['best'] or record['bits']<entry['best'][name]['bits']):
                    entry['best'][name]=dict(record,beta=beta,d=d,decode_marginal=str(p))
            entry['points'].append(point)
            args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(cost_name,{name:value['bits'] for name,value in entry['best'].items()},flush=True)
    report['complete']=True;args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
