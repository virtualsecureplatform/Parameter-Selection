#!/usr/bin/env python3
"""Validate MPFR repair of pinned hybrid probability underflow, not security."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import sys

from narrow_ckks_block_ring_joint_decode import stable_probability_cost
from narrow_ckks_block_hybrid import estimator_revision
from narrow_ckks_security import configure_max_beta


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--estimator',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    sys.path.insert(0,str(args.estimator.resolve()))
    from estimator import LWE,ND,RC
    from estimator.lwe_primal import PrimalHybrid
    from sage.all import QQ,RR,RealField,ceil,log,oo,sqrt
    configure_max_beta(4096)
    real=RealField(512)
    high=RealField(2048)
    # Independent subtraction at sufficient precision checks log1p at rare p.
    amplification=[]
    for text in ('0.25','0.5','0.99','1','1e-332'):
        p=real(text)
        low=1 if p>=real(0.99) else ceil((-real(0.99)).log1p()/(-p).log1p())
        q=high(p)
        reference=1 if q>=high(0.99) else ceil(log(1-high(0.99))/log(1-q))
        # Huge integer trial counts have rounding uncertainty; compare logs.
        error=abs(float(log(high(low)/reference,2)))
        assert error<1e-140
        if text in ('0.25','0.5','0.99','1'):
            assert low==reference
        amplification.append(dict(probability=text,trial_bits=float(log(real(low),2)),
                                  relative_log2_error=error))
    base=Path(__file__).resolve().parent
    posterior_path=base/'narrow_ckks_ll14_rotation_z14560-s5.json'
    posterior=json.loads(posterior_path.read_text())
    hist={int(k):int(v) for k,v in posterior['matching_rotation_histogram'].items()}
    union=real(QQ(posterior['expected_matching_rotations']))*sum(
        (real(v)/k for k,v in hist.items()),real(0))/sum(hist.values())
    support=QQ(posterior['rotation_labelled_candidate_count'])
    variance=QQ(19)/(16384-14560)
    params=LWE.Parameters(n=16384,q=2**64,
        Xs=ND.NoiseDistribution(n=16384,mean=0,stddev=sqrt(variance),
                               bounds=(-1,1),_density=variance),
        Xe=ND.DiscreteGaussian(3.2))
    points=[]
    for model in ('BDGL16','MATZOV'):
        for beta,d in ((174,3671),(182,3415),(206,3351),(230,3223)):
            geometry=PrimalHybrid.beta_params(beta=beta,params=params,zeta=14560,
                babai=True,mitm=True,m=oo,d=d,red_cost_model=getattr(RC,model))
            p=real(geometry['babai_probability'])*real(geometry['mitm_probability'])
            joint=union*p
            raw_probability=union*geometry['babai_probability']*geometry['mitm_probability']
            backend=PrimalHybrid.cost(beta=beta,params=params,zeta=14560,d=d,
                babai=True,mitm=True,red_cost_model=getattr(RC,model),
                search_space=support,hit_probability=union)
            repaired=stable_probability_cost(geometry,support,joint,beta,14560)
            repaired_bits=float(log(repaired['rop'],2))
            if beta==174:
                assert joint>0 and not raw_probability and backend['rop']==oo
                difference=None
            else:
                assert backend['rop']!=oo
                difference=abs(repaired_bits-float(log(backend['rop'],2)))
                assert difference<1e-10
            # Compute the same per-attempt work independently, without Cost.repeat.
            q=high(joint)
            trials=1 if q>=high(0.99) else ceil(log(1-high(0.99))/log(1-q))
            work=high(geometry['bkz_cost']['rop'])+high(geometry['svp_cost']['rop'])*high(RR(sqrt(support)))
            reference_bits=float(log(work*trials,2))
            assert abs(reference_bits-repaired_bits)<1e-10
            points.append(dict(model=model,beta=beta,d=d,joint=str(joint),
                backend_probability=str(raw_probability),
                backend_probability_parent=str(raw_probability.parent()),
                backend_missing=backend['rop']==oo,repaired_bits=repaired_bits,
                finite_control_log2_difference=difference,
                independent_reference_bits=reference_bits))
    # The repair rejects invalid probabilities rather than inventing finite costs.
    for invalid in (0,-1,2):
        try:
            stable_probability_cost(geometry,support,invalid,beta,14560)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid probability accepted')
    report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scope=__doc__,estimator_revision=estimator_revision(args.estimator),
        sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                 (Path(__file__),base/'narrow_ckks_block_ring_joint_decode.py',posterior_path)},
        amplification=amplification,points=points,passed=True,
        security_validated=False,
        cause='Sage multiplication with Real Double Field coerces the positive MPFR joint probability to double zero',
        repair='Preserve backend result; use the same geometry and cost formula with MPFR joint probability and log1p only when backend multiplication is zero')
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('PASS: rare amplification, two reproduced underflows, six finite controls, independent 2048-bit references',flush=True)


if __name__=='__main__':
    main()
