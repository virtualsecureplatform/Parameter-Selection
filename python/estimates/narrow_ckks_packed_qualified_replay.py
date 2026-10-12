#!/usr/bin/env python3
"""Replay a completed packed checkpoint with a small-qualified frozen executable.

An optional PID/start-time pair delays evaluation until the baseline controller
has exited. No key generation or selector reconstruction is requested. This
checks standard refresh and multiplication, not Table 7 functional latency or
structured-key security.
"""
import argparse
import datetime
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time

from narrow_ckks_packed_checkpoint_validate import validate
from narrow_ckks_packed_cpu_compare import digest


def identity(pid):
    try:
        fields = (Path('/proc') / str(pid) / 'stat').read_text().rsplit(')', 1)[1].split()
        return dict(pid=pid, start_ticks=int(fields[19]), state=fields[0])
    except FileNotFoundError:
        return dict(pid=pid, state='absent')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--baseline-validation', type=Path, required=True)
    parser.add_argument('--small-qualification', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wait-runner', help='Baseline controller PID:start_ticks')
    args = parser.parse_args()
    expected = None
    if args.wait_runner:
        try:
            pid, start = map(int, args.wait_runner.split(':'))
            if pid <= 0 or start < 0:
                raise ValueError()
            expected = dict(pid=pid, start_ticks=start)
        except ValueError:
            parser.error('--wait-runner must be a positive PID and nonnegative start time')
    if args.output.exists():
        parser.error('output must be a new file')
    qualification = json.loads(args.small_qualification.read_text())
    if (qualification.get('terminal') is not True or
            qualification.get('all_checks_passed') is not True or
            qualification.get('controller_actual_exit_code') != 0):
        raise ValueError('small qualification is incomplete')
    if digest(args.exe) != qualification['exe_sha256']:
        raise ValueError('executable differs from small-qualified binary')
    source_hashes = dict(qualification['source_dependencies_sha256'])
    source_hashes['TFHEpp/test/ckks/ckks_narrow_benchmark.cpp'] = qualification['source_sha256']
    if any(digest(Path(path)) != value for path, value in source_hashes.items()):
        raise ValueError('small-qualified sources changed')
    args.root.mkdir(mode=0o700)
    root = args.root.resolve()
    executable = root / 'ckks_narrow_benchmark'
    shutil.copy2(args.exe, executable)
    if digest(executable) != qualification['exe_sha256']:
        raise ValueError('pinned executable copy mismatch')
    frozen_controller = root / 'controller-at-launch.py'
    shutil.copy2(__file__, frozen_controller)
    helper_hashes = {name: digest(Path(__file__).with_name(name)) for name in (
        'narrow_ckks_packed_checkpoint_validate.py', 'narrow_ckks_prefix_arithmetic.py',
        'narrow_ckks_packed_cpu_compare.py')}
    state = dict(scope=__doc__, terminal=False, all_checks_passed=False,
                 stage='waiting for baseline controller' if expected else 'checking baseline',
                 runner=identity(os.getpid()), wait_runner=expected,
                 baseline_validation=str(args.baseline_validation.resolve()),
                 small_qualification_sha256=digest(args.small_qualification),
                 exe_sha256=digest(executable), source_sha256=source_hashes,
                 controller_sha256=digest(frozen_controller), helper_sha256=helper_hashes,
                 omp_threads=4,
                 evaluator_cap_bytes=32 * 1024**3, fresh_keys_generated=False,
                 selector_reconstruction_requested=False,
                 security_validated=False, paper_performance_reproduced=False)

    def save():
        state['updated_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        pending = args.output.with_suffix(args.output.suffix + '.pending')
        pending.write_text(json.dumps(state, indent=2, allow_nan=False) + '\n')
        pending.replace(args.output)

    def limits():
        cap = state['evaluator_cap_bytes']
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        resource.setrlimit(resource.RLIMIT_STACK, (64 * 1024**2, 64 * 1024**2))

    save()
    try:
        while expected:
            observed = identity(expected['pid'])
            state['wait_runner_observation'] = observed
            save()
            if (observed.get('start_ticks') != expected['start_ticks'] or
                    observed['state'] == 'Z'):
                break
            time.sleep(30)
        baseline = json.loads(args.baseline_validation.read_text())
        if expected and any(baseline.get('runner', {}).get(k) != v for k, v in expected.items()):
            raise ValueError('baseline controller identity mismatch')
        for field in ('terminal', 'all_passed', 'checkpoint_files_unchanged',
                      'refresh_and_square_ciphertexts_identical',
                      'refresh_and_square_keys_reused'):
            if baseline.get(field) is not True:
                raise ValueError('baseline qualification failed: ' + field)
        if baseline.get('n') not in (512, 16384):
            raise ValueError('unsupported baseline ring')
        runs = baseline['runs']
        if [run.get('label') for run in runs] != ['retained', 'replay']:
            raise ValueError('baseline is missing creation/replay evidence')
        for run in runs:
            phase = run.get('independent_phase_validation', {})
            if (run.get('actual_exit_code') != 0 or phase.get('accuracy_passed') is not True or
                    phase.get('independent_phase_validated') is not True):
                raise ValueError('baseline evaluator or independent phase checks failed')
            for artifact in phase['artifacts'].values():
                if digest(Path(artifact['path'])) != artifact['sha256']:
                    raise ValueError('baseline output artifact changed')
        checkpoint = Path(baseline['checkpoint'])
        hashes = baseline['checkpoint_hashes']
        if len(hashes) != baseline['checkpoint_file_count'] or not hashes:
            raise ValueError('incomplete baseline file manifest')
        reference = Path(next(run['report'] for run in baseline['runs'] if run['label'] == 'replay'))
        state.update(stage='evaluating qualified retained replay', n=baseline['n'],
                     checkpoint=str(checkpoint), checkpoint_file_count=len(hashes),
                     checkpoint_hashes=hashes, baseline_validation_sha256=digest(args.baseline_validation),
                     checkpoint_rehashed_before_evaluation=False)
        report = root / 'replay.json'
        log = root / 'replay-time.log'
        mode = ('--small-high-packed-bounded-replay' if baseline['n'] == 512 else
                '--ll14-packed-bounded-replay')
        command = ['/usr/bin/time', '-v', str(executable), mode, str(checkpoint), str(report)]
        with log.open('wb') as output:
            process = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT,
                env=dict(os.environ, OMP_NUM_THREADS='4', OMP_STACKSIZE='64M'), preexec_fn=limits)
            state['active_evaluator'] = identity(process.pid)
            state['command'] = command
            save()
            code = process.wait()
        state.pop('active_evaluator', None)
        state['evaluator_actual_exit_code'] = code
        save()
        if code != 0:
            raise ValueError('qualified replay evaluator failed')
        result = json.loads(report.read_text())
        for field in ('refresh_keys_reused', 'post_refresh_square_key_reused',
                      'diagnostic_preparation_reused'):
            if result.get(field) is not True:
                raise ValueError('replay reuse gate failed: ' + field)
        if result['packed_setup_seconds'] != 0 or result['selector_reconstruction_seconds'] != 0:
            raise ValueError('replay repeated setup or reconstruction')
        if result['n'] != baseline['n']:
            raise ValueError('replay ring differs from baseline')
        for field in ('coefficient_preparation_seconds', 'prepared_evaluation_seconds',
                      'selector_validation_seconds'):
            if not math.isfinite(result[field]) or result[field] <= 0:
                raise ValueError('invalid phase timing: ' + field)
        if not math.isclose(result['refresh_seconds'], result['coefficient_preparation_seconds'] +
                            result['prepared_evaluation_seconds'], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError('refresh timing excludes preparation')
        width = (baseline['n'] + 30) // 31
        tail_blocks = (baseline['n'] // 2) // (2 * (width - 1))
        groups = 31 + (31 + tail_blocks - 1) // tail_blocks
        progress = [line for line in log.read_text().splitlines()
                    if line.startswith('selector validation groups=')]
        if (len(progress) != groups or any(not line.startswith(
                f'selector validation groups={index}/{groups} ')
                for index, line in enumerate(progress, 1))):
            raise ValueError('selector validation progress is incomplete')
        state['selector_slot_comparisons'] = width * groups * (baseline['n'] // 2)
        state['phase_seconds'] = {field: result[field] for field in (
            'selector_validation_seconds', 'coefficient_preparation_seconds',
            'prepared_evaluation_seconds', 'refresh_seconds', 'key_setup_seconds')}
        phase = validate(checkpoint / 'checkpoint.bin', Path(str(report) + '.ciphertext'),
                         Path(str(report) + '.square.ciphertext'), report)
        if not phase['accuracy_passed'] or not phase['independent_phase_validated']:
            raise ValueError('independent phase checks failed')
        state.update(stage='checking retained integrity', independent_phase_validation=phase)
        save()
        current = {str(p.relative_to(checkpoint)): digest(p)
                   for p in sorted(checkpoint.rglob('*')) if p.is_file()}
        if current != hashes:
            raise ValueError('checkpoint changed since qualified baseline')
        identities = {suffix: digest(Path(str(report) + suffix)) == digest(Path(str(reference) + suffix))
                      for suffix in ('.ciphertext', '.square.ciphertext')}
        if not all(identities.values()):
            raise ValueError('earlier executable ciphertext identity failed')
        if any(digest(Path(path)) != value for path, value in source_hashes.items()):
            raise ValueError('qualified sources changed during replay')
        state.update(terminal=True, all_checks_passed=True, stage='completed',
                     checkpoint_files_unchanged=True, earlier_ciphertext_identities=identities,
                     report=str(report), report_sha256=digest(report))
    except BaseException as error:
        state.update(terminal=True, all_checks_passed=False, stage='failed', error=repr(error))
        save()
        raise
    save()
    print('Qualified retained packed replay, independent phase and integrity checks passed')


if __name__ == '__main__':
    main()
