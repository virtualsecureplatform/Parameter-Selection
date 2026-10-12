#!/usr/bin/env python3
"""Validate a small-ring comparison of the running old and current CPU builds."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re


EXPECTED=dict(profile='ll14-inspired-dd-streamed',function='identity',n=512,h=31,
    source_log_q=48,selector_log_q=379,packed_log_q=422,input_log_q=85,input_log_delta=41,
    output_log_q=100,output_log_delta=49,root_plain_log_delta=34,s2c_plain_log_delta=37,
    storage_bits=448,gaussian_sigma=3.2,gadget_bits=4,secondary_digit_bits=4,
    packed_selectors=True,file_backed_selectors=True,packed_fft_levels=2,
    compact_dd=True,compact_s2c=True,bounded_rotation_fft=True,staged_fft_caches=True,
    relinearize_before_rescale=False,auxiliary_switch_bits=0,factor_key_log_q=345,
    structured_switch_key_log_q=48,s2c_input_key_log_q=85,cache_mode='streamed',
    tree_depth=5,selector_ntt_chunk_size=17,message_bound=0.125,input_pattern='random-with-endpoints')


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda:source.read(8*1024**2),b''):h.update(chunk)
    return h.hexdigest()


def checked_report(report,n=512):
    if type(n) is not int or n not in (512,16384):
        raise ValueError('unsupported packed report ring')
    profile=dict(EXPECTED,n=n,packed_log_q=413+(n.bit_length()-1),
                 selector_ntt_chunk_size=min(40,(n+30)//31))
    for name,expected in profile.items():
        value=report.get(name)
        if value!=expected or (type(expected) in (bool,int) and type(value) is not type(expected)):
            raise ValueError('incompatible profile: '+name)
    for name,value in report.items():
        if type(value) is float and not math.isfinite(value):
            raise ValueError('nonfinite report field: '+name)
    tolerance=2**-16
    if report['evaluation_tolerance']!=tolerance or report['target_error_tolerance']!=tolerance:
        raise ValueError('changed fixed accuracy gate')
    for name in ('max_error','rms_error','fresh_input_max_error','evaluation_max_error'):
        if not 0<=report[name]<tolerance:raise ValueError('accuracy gate failed: '+name)
    if (report['selector_reconstruction_tolerance']!=2**-20 or
        not 0<=report['selector_reconstruction_max_error']<2**-20):
        raise ValueError('selector reconstruction gate failed')
    # The identity benchmark uses its fixed evaluation tolerance for both
    # square gates, rather than the propagated functional-approximation bound.
    product_tolerance=tolerance
    for name in ('post_refresh_product_max_error','post_refresh_product_evaluation_max_error'):
        if not 0<=report[name]<product_tolerance:raise ValueError('square gate failed')
    if report['accuracy_passed'] is not True:raise ValueError('reported accuracy failed')
    if report.get('accuracy_failure_injected',False) is not False:
        raise ValueError('injected accuracy failure')
    if (report.get('degree')!=0 or report.get('approximation_method')!='direct-remez' or
        report.get('retained_unit_product_levels')!=0 or report.get('retained_unit_products')!=[] or
        report.get('retained_unit_product_max_error')!=0):
        raise ValueError('packed identity report contains a different evaluation path')
    for name in ('selector_reconstruction_seconds','key_setup_seconds','refresh_seconds'):
        if not report[name]>0:raise ValueError('invalid phase timing')
    return profile


def checked(path,prefix):
    report=json.loads(path.read_text())
    checked_report(report)
    time_path=Path(str(prefix)+'-time.log');time=time_path.read_text()
    if ('Command terminated by signal' in time or
        not re.search(r'Exit status:\s*0\s*$',time) or Path(str(prefix)+'.stderr').stat().st_size):
        raise ValueError('run did not complete cleanly')
    cpu=sum(float(re.search(label+r':\s*([0-9.]+)',time)[1])
            for label in (r'User time \(seconds\)',r'System time \(seconds\)'))
    return dict(report=str(path),report_sha256=digest(path),phase_seconds={name:report[name] for name in
                ('selector_reconstruction_seconds','key_setup_seconds','refresh_seconds')},
        max_error=report['max_error'],selector_max_error=report['selector_reconstruction_max_error'],
        post_refresh_product_max_error=report['post_refresh_product_max_error'],
        cpu_seconds=cpu,peak_rss_kib=int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',time)[1]),
        time_log=str(time_path),time_log_sha256=digest(time_path),actual_exit_code=0,stderr_bytes=0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old-report',type=Path,required=True)
    parser.add_argument('--current-report',type=Path,required=True)
    parser.add_argument('--budget',type=Path,required=True)
    parser.add_argument('--old-prefix',type=Path,default=Path('/tmp/ckks-packed-cpu-old'))
    parser.add_argument('--current-prefix',type=Path,default=Path('/tmp/ckks-packed-cpu-current'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    budget=json.loads(args.budget.read_text())
    old=checked(args.old_report,args.old_prefix);current=checked(args.current_report,args.current_prefix)
    for entry in budget['binaries'].values():
        if digest(Path(entry['path']))!=entry['sha256']:raise ValueError('captured executable changed')
    root=Path(budget['private_root'])
    for name in ('selectors-old','selectors-current'):
        if (root/name).exists():raise ValueError('owned selectors were not cleaned up')
    result=dict(scope=__doc__,passed=True,old=old,current=current,profile=EXPECTED,
        launch_budget=str(args.budget),launch_budget_sha256=digest(args.budget),binaries=budget['binaries'],
        validator_sha256=digest(Path(__file__)),
        observed_selector_reconstruction_ratio=old['phase_seconds']['selector_reconstruction_seconds']/current['phase_seconds']['selector_reconstruction_seconds'],
        observed_total_cpu_ratio=old['cpu_seconds']/current['cpu_seconds'],
        limitations=['N512 identity profile only; no N16384 speed or completion claim',
            'One sample per build, independent keys/inputs, concurrent full CPU jobs',
            'Builds differ beyond packed workspace changes; no isolated causal attribution',
            'Old workspace source was not recovered; executable hash identifies the baseline'],
        security_validated=False,full_packed_validation_complete=False)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('PASS: identical fixed CPU profile, finite accuracy/square gates, normal exits and cleanup; observed reconstruction ratio',result['observed_selector_reconstruction_ratio'])


if __name__=='__main__':main()
