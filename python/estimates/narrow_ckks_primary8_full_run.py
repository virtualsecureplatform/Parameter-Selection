#!/usr/bin/env python3
"""Queue retained N16384 primary8 refresh behind a completed validation.

Freezes the small-qualified executable and independent validators before waiting.
Reuses the original state and diagnostically corrected selectors; verifies original
input identity, independent refresh/square accuracy and unchanged retained files.
This is not a security qualification or a reproduction of paper performance.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

GIB = 1 << 30


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as file:
        for block in iter(lambda: file.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def process(pid, start):
    try:
        fields = (Path('/proc') / str(pid) / 'stat').read_text().rsplit(')', 1)[1].split()
        if int(fields[19]) != start or fields[0] == 'Z':
            return None
        return dict(pid=pid, start_ticks=start, state=fields[0],
                    cpu_ticks=int(fields[11])+int(fields[12]),
                    rss_bytes=int(fields[21])*os.sysconf('SC_PAGE_SIZE'))
    except FileNotFoundError:
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--state', type=Path, required=True)
    parser.add_argument('--selectors', type=Path, required=True)
    parser.add_argument('--basis-evidence', type=Path, required=True)
    parser.add_argument('--mode', choices=['--replay-primary8', '--replay-dd-primary8'],
                        default='--replay-primary8')
    parser.add_argument('--packed-evidence', type=Path, required=True)
    parser.add_argument('--packed-process', required=True, help='PID:start_ticks')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume-waiting', action='store_true',
                        help='Resume an inactive waiting controller using its existing frozen artifacts')
    args = parser.parse_args()
    if not args.resume_waiting and (args.destination.exists() or args.output.exists()):
        parser.error('destination and output must be new; preserve existing evidence')
    repo = Path(__file__).resolve().parents[3]
    root = args.destination.resolve()
    if not args.resume_waiting:
        root.mkdir(mode=0o700)
    frozen = root / 'frozen'
    if not args.resume_waiting:
        frozen.mkdir(mode=0o700)
    pid, start = map(int, args.packed_process.split(':'))
    evidence = dict(scope=__doc__, terminal=False, all_checks_passed=False,
                    stage='freezing qualified artifacts', fresh_keys_generated=False,
                    security_validated=False, paper_performance_reproduced=False,
                    selector_bit_error_origin_established=False,
                    state=str(args.state), selectors=str(args.selectors), destination=str(root),
                    child_cap_bytes=20*GIB, parent_cap_bytes=128 << 20,
                    host_reserve_bytes=10*GIB)
    evidence['replay_mode'] = args.mode
    evidence['prerequisite_evidence'] = str(args.basis_evidence)
    if args.resume_waiting:
        evidence = json.loads(args.output.read_text())
        old = evidence['controller']
        if (evidence.get('terminal') is not False or evidence.get('active_child')
                or evidence.get('stage') != 'waiting for completed basis validation and memory budget'
                or evidence['destination'] != str(root) or evidence['state'] != str(args.state)
                or evidence['selectors'] != str(args.selectors)
                or evidence.get('replay_mode', '--replay-primary8') != args.mode
                or process(old['pid'], old['start_ticks'])):
            parser.error('only an inactive waiting controller with matching inputs may resume')
        evidence.setdefault('prior_controllers', []).append(old)

    def save():
        evidence['updated_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        pending = args.output.with_suffix('.pending')
        pending.write_text(json.dumps(evidence, indent=2, allow_nan=False)+'\n')
        pending.replace(args.output)

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (20*GIB, 20*GIB))
        resource.setrlimit(resource.RLIMIT_STACK, (64 << 20, 64 << 20))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))

    def verify_frozen():
        if any(sha(Path(path)) != digest for path, digest in evidence['frozen_sha256'].items()):
            raise ValueError('frozen executable or validation helper changed')

    def run(label, command):
        verify_frozen()
        evidence['stage'] = label
        with (root / (label+'.log')).open('w') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                     env=dict(os.environ, OMP_NUM_THREADS='4', OMP_STACKSIZE='64M'),
                                     preexec_fn=limits)
            fields = (Path('/proc') / str(child.pid) / 'stat').read_text().rsplit(')', 1)[1].split()
            evidence['active_child'] = dict(pid=child.pid, start_ticks=int(fields[19]))
            save()
            code = child.wait()
        evidence[label+'_actual_exit_code'] = code
        evidence.pop('active_child', None)
        save()
        if code != 0:
            raise ValueError(label+' failed; see '+str(root / (label+'.log')))
        print(label+' exit0', flush=True)

    originals = None
    try:
        if not args.resume_waiting:
            small_name = ('narrow-dd-auto-primary8-small-refresh-validation.json'
                          if args.mode == '--replay-primary8' else
                          'narrow-dd-primary8-dd-small-refresh-validation.json')
            small_path = repo / 'TFHEpp/docs' / small_name
            small = json.loads(small_path.read_text())
            if small.get('terminal') is not True or small.get('all_checks_passed') is not True:
                raise ValueError('small refresh qualification missing')
            if args.mode == '--replay-dd-primary8' and small.get('default_regression_exact') is not True:
                raise ValueError('default refresh regression qualification missing')
            evidence['small_qualification_sha256'] = sha(small_path)
            correction = repo / 'TFHEpp/docs/narrow-dd-aux8-selector-correction-preparation-validation.json'
            evidence['selector_correction_provenance_sha256'] = sha(correction)
            evidence['source_sha256'] = small['source_sha256']
            for name, digest in small['source_sha256'].items():
                if sha(repo / name) != digest:
                    raise ValueError('qualified source changed: '+name)
            artifacts = [(repo/'TFHEpp/build/test/ckks/ckks_narrow_auxiliary_replay', small['exe_sha256']),
                         (repo/'TFHEpp/build/test/ckks/ckks_narrow_auxiliary_square', small['square_exe_sha256'])]
            helpers = ['narrow_ckks_prefix_arithmetic.py', 'narrow_ckks_auxiliary_probe_prefix.py',
                       'narrow_ckks_auxiliary_square_validate.py', 'narrow_ckks_auxiliary_functions_validate.py']
            artifacts += [(Path(__file__).parent/name, sha(Path(__file__).parent/name)) for name in helpers]
            evidence['frozen_sha256'] = {}
            for source, digest in artifacts:
                target = frozen/source.name
                shutil.copy2(source, target)
                if sha(target) != digest:
                    raise ValueError('artifact differs from qualification: '+str(source))
                evidence['frozen_sha256'][str(target)] = digest
        else:
            verify_frozen()
        evidence['controller_sha256'] = sha(Path(__file__))
        own_fields = (Path('/proc')/str(os.getpid())/'stat').read_text().rsplit(')', 1)[1].split()
        evidence['controller'] = dict(pid=os.getpid(), start_ticks=int(own_fields[19]))
        # Keep the parent's soft cap small, allowing children to install their
        # larger cap without attempting to raise an inherited hard limit.
        resource.setrlimit(resource.RLIMIT_AS, (128 << 20, 20*GIB))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        while True:
            basis = json.loads(args.basis_evidence.read_text())
            if basis.get('terminal') is True and basis.get('all_checks_passed') is not True:
                raise ValueError('basis prerequisite terminated without validation success')
            packed = json.loads(args.packed_evidence.read_text())
            job = process(pid, start)
            remaining = (max(0, 64*GIB-job['rss_bytes']) if job else
                         0 if packed.get('terminal') is True and packed.get('all_passed') is True else 64*GIB)
            available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                                 if line.startswith('MemAvailable:')))*1024
            required = remaining + 30*GIB + (128 << 20)
            disk_free = shutil.disk_usage(root).free
            gate = basis.get('terminal') is True and basis.get('all_checks_passed') is True
            evidence.update(stage='waiting for completed basis validation and memory budget',
                            basis_gate_passed=gate, packed_process_observation=job,
                            memory_budget=dict(available_bytes=available, required_bytes=required,
                                               packed_remaining_reservation_bytes=remaining),
                            disk_budget=dict(available_bytes=disk_free, required_bytes=4*GIB))
            save()
            if gate and available >= required and disk_free >= 4*GIB:
                break
            print(evidence['stage'], flush=True)
            time.sleep(30)
        evidence['basis_qualification_sha256'] = sha(args.basis_evidence)
        evidence['stage'] = 'hashing original retained files'
        save()
        originals = [args.state]+sorted(p for p in args.selectors.rglob('*') if p.is_file())
        evidence['original_sha256'] = {str(p): sha(p) for p in originals}
        save()
        run('refresh', [str(frozen/'ckks_narrow_auxiliary_replay'), args.mode,
                        str(args.state), str(args.selectors), str(root/'refresh')])
        report = root/'refresh/report.json'
        evidence['report'] = json.loads(report.read_text())
        r = evidence['report']
        projected_flag = ('projected_primary8_conjugation' if args.mode == '--replay-primary8'
                          else 'projected_primary8_dd')
        if not (r['n'] == 16384 and r[projected_flag] is True
                and r['fresh_keys_generated'] is False and r['mode'] == args.mode
                and r['accuracy_passed'] is True):
            raise ValueError('full refresh profile or accuracy gate failed')
        prefix = root/'refresh/replayed-prefix.bin'
        phase = root/'refresh-independent.json'
        run('independent_refresh', [sys.executable, str(frozen/'narrow_ckks_auxiliary_probe_prefix.py'),
                                   '--state', str(prefix), '--report', str(report), '--output', str(phase)])
        evidence['independent_refresh'] = json.loads(phase.read_text())
        if evidence['independent_refresh']['fixed_23bit_refresh_passed'] is not True:
            raise ValueError('independent 23-bit refresh failed')
        # Identity comparison reads only the original prefix, never logs private material.
        sys.path.insert(0, str(frozen))
        from narrow_ckks_auxiliary_functions_validate import same_retained_input
        same_retained_input(args.state, prefix, 'sigmoid', 16384)
        evidence['same_retained_secrets_slots_and_input'] = True
        run('square', [str(frozen/'ckks_narrow_auxiliary_square'), str(args.state),
                       str(prefix), str(root/'square')])
        square_phase = root/'square-independent.json'
        run('independent_square', [sys.executable, str(frozen/'narrow_ckks_auxiliary_square_validate.py'),
                                  '--state', str(args.state), '--square', str(root/'square/square.bin'),
                                  '--report', str(root/'square/report.json'),
                                  '--output', str(square_phase)])
        evidence['square_report'] = json.loads((root/'square/report.json').read_text())
        evidence['independent_square'] = json.loads(square_phase.read_text())
        for metrics in (evidence['square_report'], evidence['independent_square']):
            if metrics['product_accuracy_passed'] is not True:
                raise ValueError('post-refresh product accuracy failed')
        for name in ('post_refresh_product_max_error', 'rms_error', 'target_amplitude', 'product_tolerance'):
            if not math.isclose(evidence['square_report'][name], evidence['independent_square'][name],
                                rel_tol=1e-11, abs_tol=1e-11):
                raise ValueError('independent square metric mismatch: '+name)
        evidence['prefix_sha256'] = sha(prefix)
        evidence['square_sha256'] = sha(root/'square/square.bin')
        verify_frozen()
        evidence['all_checks_passed'] = True
    except Exception as error:
        evidence['error'] = str(error)
        evidence['all_checks_passed'] = False
    finally:
        if originals is not None and 'original_sha256' in evidence:
            try:
                evidence['original_files_unchanged'] = evidence['original_sha256'] == {str(p): sha(p) for p in originals}
                evidence['original_files_checked'] = len(originals)
            except Exception as error:
                evidence['integrity_check_error'] = str(error)
                evidence['original_files_unchanged'] = False
            if not evidence['original_files_unchanged']:
                evidence['all_checks_passed'] = False
        evidence['terminal'] = True
        evidence['stage'] = 'validation passed' if evidence['all_checks_passed'] else 'validation failed'
        save()
    return 0 if evidence['all_checks_passed'] else 2


if __name__ == '__main__':
    sys.exit(main())
