#!/usr/bin/env python3
"""Replay the original retained aux8 state without root instrumentation.

Checks provenance and memory before loading keys. Reports ciphertext identity
and independently decoded accuracy separately from diagnostic process success.
"""
import argparse
import datetime
import json
import math
import os
from pathlib import Path
import resource
import subprocess

from narrow_ckks_corrected_selector_replay import digest, process, save, GIB
from narrow_ckks_auxiliary_probe_prefix import inspect


def run(a):
    os.umask(0o077)
    resource.setrlimit(resource.RLIMIT_AS, (128*2**20, resource.getrlimit(resource.RLIMIT_AS)[1]))
    provenance = json.loads(a.provenance.read_text())
    retained = json.loads(a.preparation.read_text())
    exe = a.repo/'build/test/ckks/ckks_narrow_auxiliary_replay'
    source = a.repo/'test/ckks/ckks_narrow_auxiliary_replay.cpp'
    if digest(exe) != provenance['executable_sha256'] or digest(source) != provenance['source_sha256']:
        raise ValueError('frozen replay source or executable changed')
    if digest(a.state) != retained['state_sha256']:
        raise ValueError('retained state changed')
    if a.destination.exists():
        raise ValueError('baseline destination must be new')
    jobs = []
    for reservation in a.reserved_job:
        pid, start, cap = map(int, reservation.split(':'))
        live = process(pid, start, cap*GIB)
        if live is not None:
            jobs.append(live)
    available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                         if line.startswith('MemAvailable:')))*1024
    required = sum(max(0, j['cap_bytes']-j['rss_bytes']) for j in jobs)+26*GIB+128*2**20
    budget = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), jobs=jobs,
                  available_bytes=available, required_bytes=required, host_reserve_bytes=10*GIB,
                  child_limit_bytes=16*GIB, parent_limit_bytes=128*2**20, eligible=available >= required)
    save(a.budget, budget)
    if not budget['eligible']:
        raise ValueError('insufficient memory budget')
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (16*GIB, 16*GIB))
        resource.setrlimit(resource.RLIMIT_STACK, (64*2**20, 64*2**20))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    log = a.destination.parent/(a.destination.name+'.log')
    stderr = log.with_suffix('.stderr')
    timing = log.with_suffix('.time')
    result = dict(scope=__doc__, terminal=False, mode='--replay', fresh_keys_generated=False,
                  failure_cause_established=False, security_validated=False,
                  state_sha256=retained['state_sha256'], source_sha256=provenance['source_sha256'],
                  executable_sha256=provenance['executable_sha256'], destination=str(a.destination),
                  selectors=str(a.selectors), budget=budget, log=str(log), stderr=str(stderr), time_log=str(timing))
    with log.open('x') as out, stderr.open('x') as err:
        child = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(timing), str(exe), '--replay',
                                  str(a.state), str(a.selectors), str(a.destination)], stdout=out, stderr=err,
                                 env=dict(os.environ, OMP_NUM_THREADS='4', OMP_STACKSIZE='64M'), preexec_fn=limits)
        result.update(child_pid=child.pid, stage='original untraced replay running')
        save(a.output, result)
        print(result['stage'], child.pid, flush=True)
        code = child.wait()
    result.update(terminal=True, actual_exit_code=code, stage='original untraced replay terminal',
                  utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    save(a.output, result)
    if code == 0:
        report = json.loads((a.destination/'report.json').read_text())
        phase = inspect(a.destination/'replayed-prefix.bin')
        for key in ['fresh_input_max_error', 'max_error', 'rms_error']:
            if not math.isfinite(report[key]) or not math.isclose(report[key], phase[key], rel_tol=1e-11, abs_tol=1e-11):
                raise ValueError('independent replay phase mismatch '+key)
        result.update(report=report, independent_prefix=phase,
                      saved_output_equal=report['saved_output_equal'],
                      fixed_23bit_refresh_passed=max(phase[k] for k in ['fresh_input_max_error', 'max_error', 'rms_error']) < 2**-23)
        save(a.output, result)
    print('original untraced replay exit', code, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['repo', 'state', 'selectors', 'destination', 'provenance', 'preparation', 'budget', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--reserved-job', action='append', default=[])
    args = parser.parse_args()
    try:
        run(args)
    except Exception as error:
        current = json.loads(args.output.read_text()) if args.output.exists() else {}
        current.update(terminal=True, stage='baseline diagnostic failed', error=str(error))
        save(args.output, current)
        raise


if __name__ == '__main__':
    main()
