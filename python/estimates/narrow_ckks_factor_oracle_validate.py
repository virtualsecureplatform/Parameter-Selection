#!/usr/bin/env python3
"""Check complete Q364 candidate and Q376 baseline factor-oracle evidence."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import re


def checked_run(path, log_q, encrypted):
    timing=path.with_name(path.stem+'-time.log')
    stderr=path.with_suffix('.stderr')
    time=timing.read_text()
    if ('Command terminated by signal' in time or
        not re.search(r'Exit status:\s*0\s*$',time) or stderr.stat().st_size):
        raise ValueError('run did not complete cleanly: '+str(path))
    records=[json.loads(line) for line in path.read_text().splitlines()
             if line.startswith('{')]
    phase=[r for r in records if 'phase_passed' in r]
    evaluation=[r for r in records if 'accuracy_passed' in r]
    if len(records)!=len(phase)+len(evaluation) or len(phase)!=15:
        raise ValueError('incomplete or unknown factor diagnostics')
    if len({r['factor_key_index'] for r in phase})!=15:
        raise ValueError('duplicate factor key')
    for r in phase:
        index=r['factor_key_index']
        if type(index) is not int or not 0<=index<15:
            raise ValueError('invalid factor key index')
        exponent=32767 if index==14 else pow(5,2**index,32768)
        if (r['automorphism_exponent']!=exponent or r['key_log_q']!=log_q or
            r['phase_coefficients_checked']!=(log_q//4)*16384 or
            r['phase_passed'] is not True or type(r['maximum_noise']) is not int or
            not 0<=r['maximum_noise']<=39):
            raise ValueError('factor phase diagnostic failed')
    expected={(index,kind) for index in range(15) for kind in range(3)} if encrypted else set()
    actual={(r['factor_key_index'],r['evaluation_kind']) for r in evaluation}
    if actual!=expected or len(evaluation)!=len(expected):
        raise ValueError('missing or duplicate encrypted factor cases')
    for r in evaluation:
        if r['accuracy_passed'] is not True or r['slots_checked']!=8192:
            raise ValueError('encrypted factor gate failed')
        for name in ('maximum_error','rms_error'):
            value=r[name]
            if type(value) not in (int,float) or not math.isfinite(value) or not 0<=value<2**-23:
                raise ValueError('invalid or excessive encrypted factor error')
    return dict(log=str(path),log_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        time_log=str(timing),time_log_sha256=hashlib.sha256(timing.read_bytes()).hexdigest(),
        actual_exit_code=0,stderr_bytes=0,key_log_q=log_q,phase=phase,evaluation=evaluation,
        phase_coefficients=sum(r['phase_coefficients_checked'] for r in phase),
        encrypted_slots=sum(r['slots_checked'] for r in evaluation),
        maximum_noise=max(r['maximum_noise'] for r in phase),
        maximum_error=max((r['maximum_error'] for r in evaluation),default=None),
        peak_rss_kib=int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',time)[1]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-log',type=Path,required=True)
    parser.add_argument('--baseline-log',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    candidate=checked_run(args.candidate_log,364,True)
    baseline=checked_run(args.baseline_log,376,False)
    root=Path(__file__).resolve().parents[3]
    source=root/'TFHEpp/test/ckks/ckks_narrow_dd_exact.cpp'
    executable=root/'TFHEpp/build/test/ckks/ckks_narrow_dd_exact'
    report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=__doc__,
        candidate=candidate,baseline=baseline,passed=True,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
        validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        fixed_precision_bits=23,noise_limit=39,noise_limit_is_diagnostic_12_sigma=True,
        phase_oracle='Independent exact single-prime integer NTT, explicit negacyclic substitution and low-Q key-message subtraction',
        evaluation_oracle='Direct substitution of encoded plaintext coefficients; shared CKKS encoder/decoder, no production automorphism in reference',
        candidate_keys_are_fresh=True,candidate_live_refresh_keys_checked=False,
        full_refresh_validated=False,security_validated=False)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('PASS: all15 Q364 exact key phases and45 encrypted cases; all15 retained Q376 phases')


if __name__=='__main__':
    main()
