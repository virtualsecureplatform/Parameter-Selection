#!/usr/bin/env python3
"""Run the retained one-bit selector experiment after a checked baseline.

Waits for the original process to terminate and for its authoritative completion
report. By default requires archived-output identity. The explicit matching
baseline option instead requires distinct byte-identical current traced and
untraced prefixes plus independent metrics, preserving any archive mismatch.
Rechecks frozen executable/source,
retained state, corrected selector hash and available memory before launching.
No keys or inputs are generated. Diagnostic accuracy is separate from exit code.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import time

from narrow_ckks_auxiliary_probe_prefix import inspect as inspect_prefix

GIB = 2**30


def process(pid, start, cap):
    try:
        stat = (Path('/proc')/str(pid)/'stat').read_text().rsplit(') ',1)[1].split()
    except FileNotFoundError:
        return None
    if int(stat[19]) != start or stat[0] == 'Z':
        return None
    return dict(pid=pid,start_ticks=start,state=stat[0],rss_bytes=int(stat[21])*4096,
                cpu_ticks=int(stat[11])+int(stat[12]),cap_bytes=cap)


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def save(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def baseline_gate(primary, matching=None):
    measured = primary.get('report', {})
    expected_mode = '--replay' if matching is not None else '--trace'
    def valid(record, mode):
        report = record.get('report', {})
        return (record.get('terminal') is True and record.get('actual_exit_code') == 0
                and report.get('n') == 16384 and report.get('function') == 'sigmoid'
                and report.get('auxiliary_bits') == 8 and report.get('mode') == mode
                and report.get('fresh_keys_generated') is False)
    if not valid(primary, expected_mode):
        raise ValueError('baseline terminal status or profile mismatch')
    if matching is None:
        if measured.get('saved_output_equal') is not True:
            raise ValueError('baseline did not complete with identical saved ciphertext')
        return dict(archived_output_reproduced=True)
    # A reproducible current baseline can support an intervention comparison
    # even when the archived output differs. Preserve that distinction.
    if not valid(matching, '--trace'):
        raise ValueError('matching baseline terminal status or profile mismatch')
    independent = primary.get('independent_prefix', {})
    for key in ['fresh_input_max_error', 'max_error', 'rms_error']:
        values = [measured.get(key), matching['report'].get(key), independent.get(key)]
        if any(type(x) not in (int, float) or not math.isfinite(x) for x in values):
            raise ValueError('baseline independent metric missing or nonfinite')
        if not all(math.isclose(values[0], x, rel_tol=1e-11, abs_tol=1e-11) for x in values[1:]):
            raise ValueError('baseline independent metrics differ')
    left = Path(primary['destination'])/'replayed-prefix.bin'
    right = Path(matching['directory'])/'replayed-prefix.bin'
    if left == right or left.stat().st_size != 5636193 or right.stat().st_size != 5636193:
        raise ValueError('baseline prefixes must be distinct full-ring files')
    hashed = digest(left)
    if hashed != digest(right):
        raise ValueError('current baseline replay ciphertexts differ')
    return dict(archived_output_reproduced=measured.get('saved_output_equal') is True,
                current_traced_untraced_prefixes_identical=True, prefix_sha256=hashed)


def run(a):
    os.umask(0o077)
    # Keep the parent small while allowing its child a separately checked cap.
    resource.setrlimit(resource.RLIMIT_AS,(128*2**20,resource.getrlimit(resource.RLIMIT_AS)[1]))
    result=dict(scope=__doc__,terminal=False,stage='waiting for matching baseline',
                baseline_pid=a.baseline_pid,baseline_start_ticks=a.baseline_start,
                fresh_keys_generated=False,failure_cause_established=False,security_validated=False)
    save(a.output,result)
    announced=None
    while True:
        baseline=process(a.baseline_pid,a.baseline_start,16*GIB)
        if baseline is not None:
            result.update(stage='waiting for live baseline',baseline_observation=baseline)
        elif not a.baseline_validation.exists():
            result.update(stage='waiting for baseline completion report')
        else:
            baseline_report=json.loads(a.baseline_validation.read_text())
            matching = (json.loads(a.matching_baseline_validation.read_text())
                        if a.matching_baseline_validation is not None else None)
            result['baseline_comparison'] = baseline_gate(baseline_report, matching)
            provenance=json.loads(a.provenance.read_text())
            correction=json.loads(a.correction.read_text())
            if (correction.get('phase_validation_actual_exit_code') != 0
                    or correction.get('corrected_group_phase_validation',{}).get('phase_accuracy_passed') is not True
                    or correction.get('modified_bits') != 1):
                raise ValueError('corrected selector phase evidence incomplete')
            jobs=[]
            for reservation in a.reserved_job:
                pid,start,cap=map(int,reservation.split(':'))
                observation=process(pid,start,cap*GIB)
                if observation is not None:
                    jobs.append(observation)
            available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                               if line.startswith('MemAvailable:')))*1024
            required=sum(max(0,j['cap_bytes']-j['rss_bytes']) for j in jobs)+10*GIB+16*GIB+128*2**20
            budget=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),jobs=jobs,
                        available_bytes=available,required_bytes=required,host_reserve_bytes=10*GIB,
                        child_limit_bytes=16*GIB,parent_limit_bytes=128*2**20,eligible=available>=required)
            save(a.budget,budget)
            if budget['eligible']:
                break
            result.update(stage='waiting for eligible memory budget',budget=budget)
        save(a.output,result)
        if result['stage'] != announced:
            print(result['stage'],flush=True);announced=result['stage']
        time.sleep(30)
    exe=a.repo/'build/test/ckks/ckks_narrow_auxiliary_replay'
    source=a.repo/'test/ckks/ckks_narrow_auxiliary_replay.cpp'
    if digest(exe) != provenance['executable_sha256'] or digest(source) != provenance['source_sha256']:
        raise ValueError('frozen replay source or executable changed')
    preparation=json.loads(a.preparation.read_text())
    if digest(a.state) != preparation['state_sha256']:
        raise ValueError('retained state hash mismatch')
    if a.matching_baseline_validation is not None:
        for key in ['state_sha256', 'source_sha256', 'executable_sha256']:
            expected = preparation[key] if key == 'state_sha256' else provenance[key]
            if baseline_report.get(key) != expected:
                raise ValueError('current baseline provenance mismatch '+key)
    selectors=Path(correction['directory'])
    if digest(selectors/correction['modified_group']) != correction['corrected_group_sha256']:
        raise ValueError('corrected selector group hash mismatch')
    if a.destination.exists():
        raise ValueError('replay destination must be new')
    # Recheck memory after streaming the provenance hashes.
    available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                       if line.startswith('MemAvailable:')))*1024
    if available < required:
        raise ValueError('memory budget changed before replay launch')
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(16*GIB,16*GIB))
        resource.setrlimit(resource.RLIMIT_STACK,(64*2**20,64*2**20))
        resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    log=a.destination.parent/(a.destination.name+'.log')
    stderr=log.with_suffix('.stderr');timing=log.with_suffix('.time')
    with log.open('x') as out,stderr.open('x') as err:
        mode = '--replay' if a.matching_baseline_validation is not None else '--trace'
        child=subprocess.Popen(['/usr/bin/time','-v','-o',str(timing),str(exe),mode,
                                str(a.state),str(selectors),str(a.destination)],stdout=out,stderr=err,
                               env=dict(os.environ,OMP_NUM_THREADS='4',OMP_STACKSIZE='64M'),preexec_fn=limits)
        result.update(stage='corrected replay running',mode=mode,child_pid=child.pid,
                      destination=str(a.destination),log=str(log),stderr=str(stderr),budget=budget,
                      source_sha256=provenance['source_sha256'],executable_sha256=provenance['executable_sha256'],
                      state_sha256=preparation['state_sha256'])
        save(a.output,result);print('corrected replay launched',child.pid,flush=True)
        code=child.wait()
    result.update(terminal=True,stage='corrected replay terminal',actual_exit_code=code,
                  utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    if code == 0:
        report=json.loads((a.destination/'report.json').read_text())
        phase=inspect_prefix(a.destination/'replayed-prefix.bin')
        import math
        for key in ['fresh_input_max_error','max_error','rms_error']:
            if not math.isfinite(report[key]) or not math.isclose(report[key],phase[key],rel_tol=1e-11,abs_tol=1e-11):
                raise ValueError('corrected replay independent phase mismatch '+key)
        result.update(report=report,independent_prefix=phase,
                      fixed_23bit_refresh_passed=max(phase[k] for k in ['fresh_input_max_error','max_error','rms_error'])<2**-23)
    save(a.output,result);print('corrected replay exit',code,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['repo','state','destination','baseline-validation','provenance','correction','preparation','budget','output']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--baseline-pid',type=int,required=True)
    p.add_argument('--baseline-start',type=int,required=True)
    p.add_argument('--reserved-job',action='append',default=[])
    p.add_argument('--matching-baseline-validation',type=Path,
                   help='Require distinct byte-identical traced/untraced current baselines instead of archived-output identity')
    a=p.parse_args()
    try:
        run(a)
    except Exception as error:
        current=json.loads(a.output.read_text()) if a.output.exists() else {}
        current.update(terminal=True,stage='diagnostic launcher failed',error=str(error),failure_cause_established=False)
        save(a.output,current)
        raise


if __name__ == '__main__':
    main()
