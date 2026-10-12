#!/usr/bin/env python3
"""Square a completed corrected aux8 replay using its original retained key.

Requires independent replay metrics, checks retained-state and square-program
provenance, and budgets memory before evaluation. Process success and accuracy
are reported separately. This does not establish the archived failure's cause.
"""
import argparse
import datetime
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import time

from narrow_ckks_corrected_selector_replay import GIB, digest, process, save
from narrow_ckks_auxiliary_square_validate import inspect


def replay_gate(record):
    report = record.get('report', {})
    if (record.get('terminal') is not True or record.get('actual_exit_code') != 0
            or record.get('stage') != 'corrected replay terminal'
            or report.get('n') != 16384 or report.get('function') != 'sigmoid'
            or report.get('auxiliary_bits') != 8 or report.get('mode') != '--replay'
            or report.get('fresh_keys_generated') is not False):
        raise ValueError('corrected replay is incomplete or has the wrong profile')
    phase = record.get('independent_prefix', {})
    for key in ['fresh_input_max_error', 'max_error', 'rms_error']:
        values = [report.get(key), phase.get(key)]
        if any(type(x) not in (int, float) or not math.isfinite(x) or x < 0 for x in values):
            raise ValueError('independent replay metrics missing or invalid')
        if not math.isclose(*values, rel_tol=1e-11, abs_tol=1e-11):
            raise ValueError('independent replay metrics disagree')
    passed = max(phase[k] for k in ['fresh_input_max_error', 'max_error', 'rms_error']) < 2**-23
    if record.get('fixed_23bit_refresh_passed') is not passed:
        raise ValueError('corrected replay accuracy gate disagrees with metrics')
    return passed


def run(a):
    os.umask(0o077)
    resource.setrlimit(resource.RLIMIT_AS, (128*2**20, resource.getrlimit(resource.RLIMIT_AS)[1]))
    if a.wait_for_replay is not None:
        pid, start = map(int, a.wait_for_replay.split(':'))
        while True:
            live = process(pid, start, 128*2**20)
            if live is None:
                break
            save(a.output, dict(scope=__doc__, terminal=False,
                                stage='waiting for corrected replay launcher',
                                replay_launcher=live, fresh_keys_generated=False,
                                failure_cause_established=False, security_validated=False))
            time.sleep(30)
    replay = json.loads(a.replay_validation.read_text())
    refresh_passed = replay_gate(replay)
    provenance = json.loads(a.square_provenance.read_text())
    exe = a.repo/'build/test/ckks/ckks_narrow_auxiliary_square'
    source = a.repo/'test/ckks/ckks_narrow_auxiliary_square.cpp'
    if digest(source) != provenance['source_sha256'] or digest(exe) != provenance['executable_sha256']:
        raise ValueError('retained square source or executable changed')
    if digest(a.state) != replay['state_sha256']:
        raise ValueError('retained state changed')
    prefix = Path(replay['destination'])/'replayed-prefix.bin'
    if prefix.stat().st_size != 5636193:
        raise ValueError('corrected replay prefix is not full-ring')
    if a.destination.exists():
        raise ValueError('square destination must be new')
    jobs = []
    for reservation in a.reserved_job:
        pid, start, cap = map(int, reservation.split(':'))
        live = process(pid, start, cap*GIB)
        if live is not None:
            jobs.append(live)
    available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                         if line.startswith('MemAvailable:')))*1024
    required = sum(max(0, j['cap_bytes']-j['rss_bytes']) for j in jobs)+12*GIB+128*2**20
    budget = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), jobs=jobs,
                  available_bytes=available, required_bytes=required, host_reserve_bytes=10*GIB,
                  child_limit_bytes=2*GIB, parent_limit_bytes=128*2**20, eligible=available >= required)
    save(a.budget, budget)
    if not budget['eligible']:
        raise ValueError('insufficient memory budget')
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (2*GIB, 2*GIB))
        resource.setrlimit(resource.RLIMIT_STACK, (64*2**20, 64*2**20))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    log = a.destination.parent/(a.destination.name+'.log')
    stderr = log.with_suffix('.stderr')
    result = dict(scope=__doc__, terminal=False, stage='corrected replay square running',
                  destination=str(a.destination), state_sha256=replay['state_sha256'],
                  prefix_sha256=digest(prefix), source_sha256=provenance['source_sha256'],
                  executable_sha256=provenance['executable_sha256'], budget=budget,
                  fixed_23bit_refresh_passed=refresh_passed, fresh_keys_generated=False,
                  failure_cause_established=False, security_validated=False)
    with log.open('x') as out, stderr.open('x') as err:
        child = subprocess.Popen([str(exe), str(a.state), str(prefix), str(a.destination)],
                                 stdout=out, stderr=err, preexec_fn=limits,
                                 env=dict(os.environ, OMP_NUM_THREADS='1', OMP_STACKSIZE='64M'))
        result['child_pid'] = child.pid
        save(a.output, result)
        code = child.wait()
    result.update(terminal=True, actual_exit_code=code, stage='corrected replay square terminal')
    save(a.output, result)
    if code == 0:
        report = json.loads((a.destination/'report.json').read_text())
        phase = inspect(a.state, a.destination/'square.bin')
        for key in ['post_refresh_product_max_error', 'rms_error', 'target_amplitude', 'product_tolerance']:
            if not math.isfinite(report[key]) or not math.isclose(report[key], phase[key], abs_tol=1e-11, rel_tol=1e-11):
                raise ValueError('independent square metric mismatch '+key)
        for key in ['n', 'function', 'log_q', 'log_delta', 'product_accuracy_passed']:
            if report[key] != phase[key]:
                raise ValueError('independent square profile or gate mismatch '+key)
        result.update(report=report, independent_phase=phase,
                      square_sha256=digest(a.destination/'square.bin'),
                      corrected_refresh_and_square_passed=refresh_passed and phase['product_accuracy_passed'])
        save(a.output, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['repo', 'state', 'destination', 'replay-validation', 'square-provenance', 'budget', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--reserved-job', action='append', default=[])
    parser.add_argument('--wait-for-replay', help='Wait for matching replay launcher PID:START_TICKS before checking completed validation')
    a = parser.parse_args()
    try:
        run(a)
    except Exception as error:
        record = json.loads(a.output.read_text()) if a.output.exists() else {}
        record.update(terminal=True, stage='corrected square diagnostic failed', error=str(error))
        save(a.output, record)
        raise


if __name__ == '__main__':
    main()
