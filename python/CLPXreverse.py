#!/usr/bin/env python3
"""Compose conditional whole-run CLPX reverse error envelopes for 16 output bits."""
import argparse
import json
import hashlib
import subprocess
from pathlib import Path
from noiseestimation.clpx_reverse_bound import estimate

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-log2-variance',type=float,default=109.20114828158572,
                   help='CLPX coefficient variance in 64-bit integer units; a tail assumption, not a certificate')
    p.add_argument('--basebit',type=int,choices=(2,4),default=4)
    p.add_argument('--half-n',type=int,default=760)
    p.add_argument('--half-alpha-bits',type=int,default=17)
    p.add_argument('--ks20-levels',type=int,default=7)
    p.add_argument('--ks20-basebit',type=int,default=2)
    p.add_argument('--multiplications',type=int,default=256)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--tfhepp',type=Path,help='Record hashes of the analyzed C++ source')
    a=p.parse_args()
    if a.half_n<=0 or not 1<=a.half_alpha_bits<=31:
        p.error('positive half dimension and alpha bits in 1..31 required')
    r=estimate(input_sigma=2**(a.input_log2_variance/2-64),basebit=a.basebit,
               half_n=a.half_n,half_alpha_bits=a.half_alpha_bits,multiplications=a.multiplications,
               ks20_levels=a.ks20_levels,ks20_basebit=a.ks20_basebit)
    r['input_log2_integer_variance']=a.input_log2_variance
    root=Path(__file__).resolve().parent.parent
    r['parameter_selection_base_revision']=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
    r['analysis_sha256']={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__).resolve(),root/'python/noiseestimation/clpx_reverse_bound.py']}
    r['is_measured_parameter_configuration']=(a.basebit,a.half_n,a.half_alpha_bits,a.ks20_levels,a.ks20_basebit)==(4,760,17,7,2)
    if a.tfhepp:
        names=['include/clpx/bfv-clpx.hpp','include/tfhe/homdecomp.hpp','include/tfhe/gatebootstrapping.hpp','include/tfhe/keyswitch.hpp','include/params/128bit.hpp','include/clpx/params/SS2CLPX.hpp']
        r['tfhepp_source_sha256']={f:hashlib.sha256((a.tfhepp/f).read_bytes()).hexdigest() for f in names}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ('parameters','maximum_uniform_cutoff_sigma',
        'cutoff_required_for_target','conditional_log2_failure_bound','conditional_target_passes',
        'bottleneck','implementation_certified')},indent=2))
    return 0
if __name__=='__main__':raise SystemExit(main())
