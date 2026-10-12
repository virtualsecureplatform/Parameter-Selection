#!/usr/bin/env python3
"""Run full retained aux8 GeLU/tanh with checked corrected selectors.

Reuses the independently validated sigmoid fixture's original keys and input.
Each function needs fresh provenance and memory checks, and its refresh and
square outputs are independently decoded before completion is reported.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import resource
import subprocess

from narrow_ckks_corrected_selector_replay import GIB, digest, process, save
from narrow_ckks_auxiliary_functions_validate import validate
from narrow_ckks_selector_integrity import snapshot, verify


def provenance_gate(build, retained, correction):
    if (build.get('terminal') is not True or build.get('build_actual_exit_code') != 0
            or build.get('all_small_functions_passed') is not True
            or set(build.get('functions', {})) != {'sigmoid', 'gelu', 'tanh'}):
        raise ValueError('retained-function build and small validation incomplete')
    if (retained.get('terminal') is not True or retained.get('actual_exit_code') != 0
            or retained.get('fixed_23bit_refresh_passed') is not True
            or retained.get('report', {}).get('n') != 16384
            or retained.get('report', {}).get('function') != 'sigmoid'
            or retained.get('report', {}).get('auxiliary_bits') != 8
            or retained.get('fresh_keys_generated') is not False):
        raise ValueError('validated full retained sigmoid fixture missing')
    phase = correction.get('corrected_group_phase_validation', {})
    if (correction.get('modified_bits') != 1 or correction.get('modified_ciphertexts') != 1
            or correction.get('phase_validation_actual_exit_code') != 0
            or phase.get('n') != 16384 or phase.get('phase_accuracy_passed') is not True):
        raise ValueError('full corrected selector phase evidence missing')


def selector_integrity_gate(selectors, manifest_path, manifest_digest):
    if digest(manifest_path) != manifest_digest:
        raise ValueError('selector integrity baseline manifest changed')
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('n') != 16384:
        raise ValueError('full retained-function selector manifest profile mismatch')
    return verify(selectors, manifest)


def run(a):
    os.umask(0o077)
    resource.setrlimit(resource.RLIMIT_AS, (128*2**20, resource.getrlimit(resource.RLIMIT_AS)[1]))
    build = json.loads(a.build_validation.read_text())
    retained = json.loads(a.retained_validation.read_text())
    correction = json.loads(a.correction.read_text())
    provenance_gate(build, retained, correction)
    selectors = Path(correction['directory'])
    source = a.repo/'test/ckks/ckks_narrow_auxiliary_functions.cpp'
    exe = a.repo/'build/test/ckks/ckks_narrow_auxiliary_functions'
    if a.destination.exists():
        raise ValueError('retained-function root must be new')
    a.destination.mkdir(mode=0o700)
    result = dict(scope=__doc__, terminal=False, destination=str(a.destination), functions={},
                  source_sha256=build['source_sha256'], executable_sha256=build['executable_sha256'],
                  state_sha256=retained['state_sha256'], fresh_keys_generated=False,
                  selector_bit_error_origin_established=False, security_validated=False)
    save(a.output, result)
    manifest_path = a.selector_manifest or a.destination/'selector-integrity.json'
    if a.selector_manifest is None:
        manifest = snapshot(selectors, 16384)
        with manifest_path.open('x') as file:
            json.dump(manifest, file, indent=2, allow_nan=False)
            file.write('\n')
    manifest_digest = digest(manifest_path)
    result.update(selector_integrity_manifest=str(manifest_path),
                  selector_integrity_manifest_sha256=manifest_digest,
                  selector_integrity_baseline_scope='Observation time only; not proof of generation-time correctness')
    save(a.output, result)
    for function in a.functions:
        if digest(source) != build['source_sha256'] or digest(exe) != build['executable_sha256']:
            raise ValueError('retained-function source or executable changed')
        if digest(a.state) != retained['state_sha256']:
            raise ValueError('original retained state changed')
        if digest(selectors/correction['modified_group']) != correction['corrected_group_sha256']:
            raise ValueError('corrected selector group changed')
        integrity_before = selector_integrity_gate(selectors, manifest_path, manifest_digest)
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
            raise ValueError('insufficient full retained-function memory budget')
        def limits():
            resource.setrlimit(resource.RLIMIT_AS, (16*GIB, 16*GIB))
            resource.setrlimit(resource.RLIMIT_STACK, (64*2**20, 64*2**20))
            resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        destination = a.destination/function
        with (a.destination/(function+'.log')).open('x') as log, (a.destination/(function+'.stderr')).open('x') as err:
            child = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(a.destination/(function+'.time')),
                                      str(exe), function, str(a.state), str(selectors), str(destination)],
                                     stdout=log, stderr=err, preexec_fn=limits,
                                     env=dict(os.environ, OMP_NUM_THREADS='4', OMP_STACKSIZE='64M'))
            result.update(stage='retained '+function+' running', current_function=function, child_pid=child.pid)
            result['functions'][function] = dict(terminal=False, budget=budget, directory=str(destination),
                                                selector_integrity_before=integrity_before)
            save(a.output, result)
            print(result['stage'], child.pid, flush=True)
            code = child.wait()
        entry = result['functions'][function]
        entry.update(terminal=True, actual_exit_code=code)
        save(a.output, result)
        entry['selector_integrity_after'] = selector_integrity_gate(selectors, manifest_path, manifest_digest)
        save(a.output, result)
        if code != 0:
            raise ValueError('retained-function evaluator failed '+function)
        entry.update(report=json.loads((destination/'report.json').read_text()),
                     independent_validation=validate(a.state, destination),
                     prefix_sha256=digest(destination/'replayed-prefix.bin'), square_sha256=digest(destination/'square.bin'))
        save(a.output, result)
        print(function, 'independent accuracy', entry['independent_validation']['accuracy_passed'], flush=True)
    result.update(terminal=True, stage='retained functions terminal',
                  all_requested_functions_passed=all(x['independent_validation']['accuracy_passed'] for x in result['functions'].values()))
    save(a.output, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['repo', 'state', 'destination', 'build-validation', 'retained-validation', 'correction', 'budget', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--functions', nargs='+', choices=['sigmoid', 'gelu', 'tanh'], default=['gelu', 'tanh'])
    parser.add_argument('--reserved-job', action='append', default=[])
    parser.add_argument('--selector-manifest', type=Path,
                        help='Existing selector hash baseline; default records a new baseline before this run')
    a = parser.parse_args()
    try:
        if len(set(a.functions)) != len(a.functions):
            raise ValueError('duplicate requested functions')
        run(a)
    except Exception as error:
        result = json.loads(a.output.read_text()) if a.output.exists() else {}
        result.update(terminal=True, stage='retained-function launcher failed', error=str(error))
        save(a.output, result)
        raise


if __name__ == '__main__':
    main()
