#!/usr/bin/env python3
"""Fixed-point structured ring cost sensitivity using measured joint rotations.

Uses the attacker's full square-root speedup on the rotation-labelled guessing
set as a heuristic. Actual compatible MITM dictionaries are not constructed.
This is neither a global attack optimization nor a security certificate.
"""
import argparse
import datetime
from fractions import Fraction
import hashlib
import json
import math
import itertools
from pathlib import Path
import sys
from narrow_ckks_block_hybrid import cost_record, estimator_revision
from narrow_ckks_security import configure_max_beta
from narrow_ckks_block_hybrid_grid import offsets


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--estimator',type=Path,required=True)
    parser.add_argument('--profile-report',type=Path,required=True)
    parser.add_argument('--probability-report',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--beta-offsets',type=offsets,default=[0])
    parser.add_argument('--d-offsets',type=offsets,default=[0])
    args=parser.parse_args()
    sys.path.insert(0,str(args.estimator.resolve()))
    from estimator import LWE,ND,RC
    from estimator.lwe_primal import PrimalHybrid
    from sage.all import QQ,RealField,sqrt,log
    configure_max_beta(4096)
    profile=json.loads(args.profile_report.read_text());prob=json.loads(args.probability_report.read_text())
    b=profile['attacks']['bdd_mitm_hybrid']['best']
    if (not prob['complete'] or prob['n']!=profile['n'] or prob['widths']!=profile['block_widths'] or
        prob['guessed_per_block']!=b['guessed_per_block'] or prob['guessed_weight']!=b['guessed_weight'] or
        prob['profile_report_sha256']!=hashlib.sha256(args.profile_report.read_bytes()).hexdigest()):
        raise ValueError('joint probability report does not match the saved fixed point')
    n=profile['n'];weight=b['guessed_weight'];norm=profile['h']-weight
    variance=QQ(norm)/(n-b['grid_zeta'])
    if not norm:raise ValueError('zero residual requires separate exhaustive model')
    params=LWE.Parameters(n=n,q=2**profile['log_q'],Xs=ND.NoiseDistribution(n=n,mean=0,stddev=sqrt(variance),bounds=(-1,1),_density=variance),Xe=ND.DiscreteGaussian(profile['sigma']))
    report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=__doc__,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),estimator_revision=estimator_revision(args.estimator),
        profile_sha256=hashlib.sha256(args.profile_report.read_bytes()).hexdigest(),
        probability_report=str(args.probability_report),probability_report_sha256=hashlib.sha256(args.probability_report.read_bytes()).hexdigest(),
        probability_seed=prob['seed'],point=dict(beta=b['grid_beta'],zeta=b['grid_zeta'],d=b['grid_d'],weight=weight,layout=b['layout']),
        models={},grid={},beta_offsets=args.beta_offsets,d_offsets=args.d_offsets,
        complete=False,security_validated=False,complete_attack_optimization=False,
        limitations=['Full square-root speedup assumed on labelled dictionaries, without constructing compatible isometric decompositions',
            'Only specified beta/d offsets at one zeta/layout/weight; other choices could be cheaper',
            'Scalar GSA residual variance omits heterogeneous coordinates and conditioned noise correlations',
            'Monte Carlo/Hoeffding probability intervals assume iid posterior draws; not security confidence intervals',
            'Correlated evaluation-key and other ring attack models remain unqualified'])
    scenarios=[('unrotated_control',int(b['candidate_count']),QQ(int(b['hit_probability_numerator']))/int(b['hit_probability_denominator']))]
    for name,p in (('joint_rotation_central',prob['union_probability_estimate']),
                   ('joint_rotation_iid95_lower',prob['iid_hoeffding95_interval'][0]),
                   ('joint_rotation_iid95_upper',prob['iid_hoeffding95_interval'][1]),
                   ('independent_rotation_sensitivity',prob['independent_rotation_probability'])):
        scenarios.append((name,int(prob['rotation_labelled_candidate_count']),RealField(256)(p)))
    for cost_name in ('BDGL16','MATZOV'):
        report['models'][cost_name]={}
        for name,support,p in scenarios:
            result=PrimalHybrid.cost(beta=b['grid_beta'],params=params,zeta=b['grid_zeta'],d=b['grid_d'],
                babai=True,mitm=True,red_cost_model=getattr(RC,cost_name),search_space=QQ(support),hit_probability=p)
            record=cost_record(result,log)
            report['models'][cost_name][name]=dict(record,support=str(support),hit_probability=str(p))
            print(cost_name,name,record['bits'],flush=True)
    report['legacy_control_identical']=report['models']['BDGL16']['unrotated_control']['bits']==b['bits']
    if not report['legacy_control_identical']:raise ValueError('legacy cost control changed')
    for cost_name in ('BDGL16','MATZOV'):
        entry=dict(points=[],best_central=None,best_iid95_attacker_favorable=None,skipped_invalid_points=0)
        report['grid'][cost_name]=entry
        for beta_offset,d_offset in itertools.product(args.beta_offsets,args.d_offsets):
            beta=b['grid_beta']+beta_offset;d=b['grid_d']+d_offset
            if not 40<=beta<=4096 or d<max(beta,n-b['grid_zeta']+1):
                entry['skipped_invalid_points']+=1
                continue
            point=dict(beta=beta,d=d,costs={})
            for name,p in (('central',prob['union_probability_estimate']),('iid95_attacker_favorable',prob['iid_hoeffding95_interval'][1])):
                result=PrimalHybrid.cost(beta=beta,params=params,zeta=b['grid_zeta'],d=d,
                    babai=True,mitm=True,red_cost_model=getattr(RC,cost_name),
                    search_space=QQ(int(prob['rotation_labelled_candidate_count'])),hit_probability=RealField(256)(p))
                record=cost_record(result,log);point['costs'][name]=record
                best_key='best_'+name
                if record['bits'] is not None and (entry[best_key] is None or record['bits']<entry[best_key]['bits']):
                    entry[best_key]=dict(record,beta=beta,d=d)
            entry['points'].append(point)
            args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(cost_name,'grid best',entry['best_central']['bits'] if entry['best_central'] else None,flush=True)
    report['complete']=True
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
