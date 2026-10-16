#!/usr/bin/env python3
"""Bounded ring/block hybrid sensitivities with exact rotation event masses.

Use the optimistic joint success bound min(1, sum_r Pr[event_r]*p_decode).
Compare a rotation-labelled dictionary with an optimistic deduplicated size:
the largest individual rotation's support is a lower bound on its union.
The latter gives an optimistic cost sensitivity, not a realizable attack.
Scalar GSA, decoder heuristics and square-root dictionary search are assumed.
No result is a security certificate or a global attack bound.
"""
import argparse
from collections import Counter
import datetime
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
from narrow_ckks_block_hybrid import layout, multiply_linear, estimator_revision
from narrow_ckks_block_hybrid_grid import saturated_layout
from narrow_ckks_block_ring_joint_decode import stable_probability_cost


def rotation_counts(widths, guessed):
    n, h = sum(widths), len(widths)
    if len(guessed) != h or any(not 0 <= g <= w for w,g in zip(widths,guessed)):
        raise ValueError('invalid guessing geometry')
    starts = list(itertools.accumulate([0]+widths))[:-1]
    mask = [int(j < g) for w,g in zip(widths,guessed) for j in range(w)]
    prefix = [0]+list(itertools.accumulate(mask+mask))
    # Polynomial products commute: group identical multisets of (width,count).
    geometries = Counter(tuple(sorted((w,prefix[(start+r)%n+w]-prefix[(start+r)%n])
                                     for start,w in zip(starts,widths))) for r in range(n))
    masses, labelled, largest = [0]*(h+1), [0]*(h+1), [0]*(h+1)
    for geometry, multiplicity in geometries.items():
        mass, support = [1], [1]
        for width,g in geometry:
            mass = multiply_linear(mass,width-g,g)
            support = multiply_linear(support,int(g<width),2*g)
        for weight in range(h+1):
            masses[weight] += multiplicity*mass[weight]
            labelled[weight] += multiplicity*support[weight]
            largest[weight] = max(largest[weight],support[weight])
    denominator = math.prod(widths)
    return dict(event_sums=[str(Fraction(m,denominator)) for m in masses],
                labelled_support=labelled, largest_rotation_support=largest,
                distinct_rotation_geometries=len(geometries))


def self_test():
    widths, guessed = [2,3,1], [1,2,0]
    n, h = sum(widths), len(widths)
    starts = list(itertools.accumulate([0]+widths))[:-1]
    mask = [int(j<g) for w,g in zip(widths,guessed) for j in range(w)]
    secrets = list(itertools.product(*(range(start,start+w) for start,w in zip(starts,widths))))
    labelled, largest, event_counts = [0]*(h+1), [0]*(h+1), [0]*(h+1)
    for r in range(n):
        projections = [set() for _ in range(h+1)]
        for secret in secrets:
            weight = sum(mask[(p+r)%n] for p in secret)
            event_counts[weight] += 1
            for signs in itertools.product((-1,1),repeat=h):
                projection = [0]*n
                for p,sign in zip(secret,signs):
                    q = (p+r)%n
                    if mask[q]:projection[q] = sign if p+r<n else -sign
                projections[weight].add(tuple(projection))
        for weight in range(h+1):
            labelled[weight] += len(projections[weight])
            largest[weight] = max(largest[weight],len(projections[weight]))
    result = rotation_counts(widths,guessed)
    assert result['labelled_support'] == labelled
    assert result['largest_rotation_support'] == largest
    assert list(map(Fraction,result['event_sums'])) == [Fraction(c,len(secrets)) for c in event_counts]
    print('Exact rotation masses and labelled/largest signed supports match exhaustive secrets',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--estimator',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--n',type=int,default=16384)
    parser.add_argument('--h',type=int,default=31)
    parser.add_argument('--log-q',type=int,default=56)
    parser.add_argument('--noise-stddev',type=float,default=2.828301687834186)
    parser.add_argument('--zetas',type=int,nargs='+',default=[13824,14336,14848,15360,15872])
    parser.add_argument('--weights',type=int,nargs='+',default=list(range(11,30,2)))
    parser.add_argument('--betas',type=int,nargs='+',default=[200,280,360,440,520])
    parser.add_argument('--dimension-offsets',type=int,nargs='+',default=[500,1000,1500,2000])
    parser.add_argument('--cost-model',choices=['MATZOV','BDGL16'],default='MATZOV')
    args = parser.parse_args();self_test()
    if args.self_test:return
    if args.estimator is None or args.output is None:parser.error('--estimator and --output required')
    if not (1<args.h<args.n and 1<=args.log_q<=448 and math.isfinite(args.noise_stddev) and args.noise_stddev>0):parser.error('invalid ring/noise parameters')
    if any(not 0<z<args.n for z in args.zetas) or any(not 0<=w<args.h for w in args.weights):parser.error('invalid guessing dimensions/weights; zero residual needs a separate model')
    if any(not 40<=b<=4096 for b in args.betas) or any(d<=0 for d in args.dimension_offsets):parser.error('invalid beta/dimension offsets')
    sys.path.insert(0,str(args.estimator.resolve()))
    from estimator import LWE,ND,RC
    from estimator.lwe_primal import PrimalHybrid
    from narrow_ckks_security import configure_max_beta
    from sage.all import QQ,RealField,sqrt,log,oo
    configure_max_beta(4096);real=RealField(512);model=getattr(RC,args.cost_model)
    width=(args.n+args.h-1)//args.h
    widths=[min(width,args.n-j*width) for j in range(args.h)]
    if min(widths)<1:raise ValueError('invalid ragged public blocks')
    report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=__doc__,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                estimator_revision=estimator_revision(args.estimator),
                settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
                widths=widths,profiles=[],points=[],best={},complete=False,security_validated=False)
    def save():args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    save()
    for zeta in args.zetas:
        profiles=[]
        for strategy in ('balanced','prefix','saturated-5-prefix'):
            guessed = saturated_layout(widths,zeta,5) if strategy=='saturated-5-prefix' else layout(widths,zeta,strategy)
            x=dict(zeta=zeta,layout=strategy,guessed_per_block=guessed,**rotation_counts(widths,guessed));profiles.append(x);report['profiles'].append(x)
        for weight in args.weights:
            variance=QQ(args.h-weight)/(args.n-zeta)
            params=LWE.Parameters(n=args.n,q=2**args.log_q,Xs=ND.NoiseDistribution(n=args.n,mean=0,stddev=sqrt(variance),bounds=(-1,1),_density=variance),Xe=ND.DiscreteGaussian(args.noise_stddev))
            for beta,offset in itertools.product(args.betas,args.dimension_offsets):
                d=args.n-zeta+offset
                if d<max(beta,args.n-zeta+1):continue
                geometry=PrimalHybrid.beta_params(beta=beta,params=params,zeta=zeta,babai=True,mitm=True,m=oo,d=d,red_cost_model=model)
                if 'babai_probability' not in geometry:
                    report['points'].append(dict(zeta=zeta,weight=weight,beta=beta,d=d,finite=False,reason='decoder geometry unavailable'));continue
                marginal=real(geometry['babai_probability'])*real(geometry['mitm_probability'])
                if not 0<marginal<=1:
                    report['points'].append(dict(zeta=zeta,weight=weight,beta=beta,d=d,finite=False,reason='zero or invalid decoder marginal'));continue
                for profile in profiles:
                    event_sum=real(QQ(profile['event_sums'][weight]));joint=min(real(1),event_sum*marginal)
                    if not joint:continue
                    for kind in ('labelled_support','largest_rotation_support'):
                        support=profile[kind][weight]
                        if not support:raise ValueError('nonzero event mass without signed support')
                        result=stable_probability_cost(geometry,QQ(support),joint,beta,zeta);bits=float(log(result['rop'],2))
                        if not math.isfinite(bits):raise ValueError('nonfinite reconstructed cost')
                        point=dict(zeta=zeta,weight=weight,beta=beta,d=d,layout=profile['layout'],dictionary=kind,bits=bits,decode_marginal=str(marginal),joint_success_upper=str(joint))
                        report['points'].append(point)
                        if kind not in report['best'] or bits<report['best'][kind]['bits']:report['best'][kind]=point
            save();print('zeta',zeta,'weight',weight,'best', {k:v['bits'] for k,v in report['best'].items()},flush=True)
    report['complete']=True;save()


if __name__=='__main__':main()
