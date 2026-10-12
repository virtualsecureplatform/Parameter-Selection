#!/usr/bin/env python3
"""Validate current full-ring root reuse against a completed retained checkpoint.

Freeze the small-qualified executable and independent decoder before waiting.
Wait for the original basis validation and termination of preceding large tests.
Reuse the checkpoint, compare all output bytes, decode phases independently and
check retained payload integrity. No fresh key generation or security claim.
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


def validate_root_tile_timings(row):
    stages = row.get('root_tile_stage_seconds', {})
    expected = {'input_and_factor_cache', 'leaf_factors_and_inline_alignment',
                'bounded_tail_alignment', 'dd_cache_preparation', 'dd_trees_and_beta'}
    if (set(stages) != expected or not all(
            isinstance(value, (float, int)) and math.isfinite(value) and value >= 0
            for value in stages.values())):
        raise ValueError('invalid exclusive root-tile timings: '+row['function'])
    if sum(stages.values()) > row['root_generation_seconds'] + 1e-3:
        raise ValueError('root-tile timings exceed root generation: '+row['function'])
    if (row['newly_evaluated_roots'] == 0 and any(stages.values())):
        raise ValueError('reused roots report new tile work: '+row['function'])
    if (sum(row['factor_stage_seconds'].values()) >
            stages['leaf_factors_and_inline_alignment'] + 1e-3):
        raise ValueError('factor timings exceed leaf scope: '+row['function'])


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as file:
        for block in iter(lambda: file.read(1 << 20), b''):
            result.update(block)
    return result.hexdigest()


def process(pid, start):
    try:
        fields = (Path('/proc')/str(pid)/'stat').read_text().rsplit(')', 1)[1].split()
        if int(fields[19]) != start or fields[0] == 'Z':
            return None
        return dict(pid=pid, start_ticks=start, state=fields[0],
                    rss_bytes=int(fields[21])*os.sysconf('SC_PAGE_SIZE'),
                    cpu_ticks=int(fields[11])+int(fields[12]))
    except FileNotFoundError:
        return None


def validated_payload_baseline(reference, checkpoint, baseline, originals):
    """Reuse a completed prerequisite's hashes; post-evaluation hashing stays mandatory."""
    if (reference.get('terminal') is not True or reference.get('all_checks_passed') is not True
            or reference.get('checkpoint_files_unchanged') is not True
            or reference.get('checkpoint_files_checked') != len(baseline)
            or Path(reference.get('checkpoint', '')).resolve() != checkpoint.resolve()
            or reference.get('original_sha256') != baseline
            or {str(p) for p in originals} != set(baseline)
            or any(p.stat().st_size != baseline[str(p)]['bytes'] for p in originals)):
        raise ValueError('completed prerequisite does not qualify this retained payload baseline')
    return baseline


def validated_frozen_resume(path, checkpoint, replay_mode, artifacts):
    """Recover a qualified queue that stopped at the baseline gate, before replay."""
    previous = json.loads(path.read_text())
    controller = previous.get('controller', {})
    active = previous.get('active_child')
    if (previous.get('terminal') is not True or previous.get('all_checks_passed') is not False
            or previous.get('error') !=
               'completed prerequisite does not qualify this retained payload baseline'
            or 'replay_actual_exit_code' in previous
            or Path(previous.get('checkpoint', '')).resolve() != checkpoint.resolve()
            or previous.get('replay_mode') != replay_mode
            or not controller or process(controller['pid'], controller['start_ticks'])
            or (active and process(active['pid'], active['start_ticks']))):
        raise ValueError('frozen resume requires a stopped pre-replay baseline-gate failure')
    frozen = previous.get('frozen_sha256', {})
    if len(frozen) != 7 or any(digest(Path(p)) != value for p, value in frozen.items()):
        raise ValueError('previous frozen qualification artifacts changed')
    for source, receipt in artifacts:
        expected = previous.get(receipt)
        if expected is None or frozen.get(str(source.resolve())) != expected:
            raise ValueError('resume artifact differs from previous qualification: '+str(source))
    return previous


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('destination', 'reference-root', 'reference-evidence',
                 'packed-evidence', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--wait-terminal', type=Path, action='append', required=True)
    parser.add_argument('--packed-process', required=True, help='PID:start_ticks')
    parser.add_argument('--small-validation', type=Path,
                        help='Qualified retained small replay evidence; defaults to original root reuse')
    parser.add_argument('--exe', type=Path,
                        help='Qualified replay executable; defaults to the all-profile benchmark')
    parser.add_argument('--factor-validation', type=Path,
                        help='Completed full22 factor qualification, required for frequency-union replay')
    parser.add_argument('--replay-mode', default='--lean-basis-fbs-replay',
                        choices=['--lean-basis-fbs-replay', '--lean-basis-fbs-replay-secondary11',
                                 '--lean-basis-fbs-replay-secondary11-selector56'])
    parser.add_argument('--dd-cache-validation', type=Path,
                        help='Retained all-level secondary11 DD qualification, required for that replay mode')
    parser.add_argument('--reuse-validated-baseline', action='store_true',
                        help='Reuse the completed prerequisite payload hashes; always rehash after evaluation')
    parser.add_argument('--resume-frozen-validation', type=Path,
                        help='Reuse verified artifacts from a stopped pre-replay baseline-gate failure')
    parser.add_argument('--omp-threads', type=int, default=4,
                        help='OpenMP worker count for the retained replay (default: 4)')
    args = parser.parse_args()
    if args.omp_threads < 1:
        parser.error('--omp-threads must be positive')
    selector56 = args.replay_mode.endswith('-secondary11-selector56')
    secondary11 = args.replay_mode.endswith('-secondary11') or selector56
    if secondary11 and args.dd_cache_validation is None:
        parser.error('secondary11 replay requires --dd-cache-validation')
    if args.destination.exists() or args.output.exists():
        parser.error('destination and evidence must be new')
    os.umask(0o077)
    repo = Path(__file__).resolve().parents[3]
    root = args.destination.resolve()
    root.mkdir(mode=0o700)
    frozen = root/'frozen'
    frozen.mkdir(mode=0o700)
    checkpoint = args.reference_root/'checkpoint'
    pid, start = map(int, args.packed_process.split(':'))
    e = dict(scope=__doc__, terminal=False, all_checks_passed=False,
             stage='freezing small-qualified artifacts', fresh_keys_generated=False,
             checkpoint=str(checkpoint), reference_root=str(args.reference_root),
             child_cap_bytes=40*GIB, parent_soft_cap_bytes=128 << 20,
             host_reserve_bytes=10*GIB, full_encrypted_root_reuse_validated=False,
             security_validated=False, paper_performance_reproduced=False)
    e['replay_mode'] = args.replay_mode
    e['omp_num_threads'] = args.omp_threads

    def save():
        e['updated_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        pending = args.output.with_suffix('.pending')
        pending.write_text(json.dumps(e, indent=2, allow_nan=False)+'\n')
        pending.replace(args.output)

    def verify_frozen():
        if any(digest(Path(path)) != value for path, value in e['frozen_sha256'].items()):
            raise ValueError('frozen executable or decoder changed')

    def run(label, command, cap):
        def limits():
            resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
            resource.setrlimit(resource.RLIMIT_STACK, (64 << 20, 64 << 20))
            resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        verify_frozen()
        e['stage'] = label
        with (root/(label+'.log')).open('x') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                     env=dict(os.environ, OMP_NUM_THREADS=str(args.omp_threads), OMP_STACKSIZE='64M'),
                                     preexec_fn=limits)
            fields = (Path('/proc')/str(child.pid)/'stat').read_text().rsplit(')', 1)[1].split()
            e['active_child'] = dict(pid=child.pid, start_ticks=int(fields[19]))
            save()
            code = child.wait()
        e.pop('active_child', None)
        e[label+'_actual_exit_code'] = code
        save()
        if code != 0:
            raise ValueError(label+' failed; see '+str(root/(label+'.log')))

    originals = None
    try:
        small_path = args.small_validation or repo/'TFHEpp/docs/narrow-dd-basis-root-reuse-validation.json'
        resumed = None
        if args.resume_frozen_validation:
            if not secondary11 or any(value is None for value in
                    (args.exe, args.small_validation, args.dd_cache_validation, args.factor_validation)):
                raise ValueError('frozen resume requires explicit secondary11 qualification artifacts')
            receipts = [(small_path, 'small_qualification_sha256'),
                        (args.dd_cache_validation, 'dd_cache_qualification_sha256'),
                        (args.factor_validation, 'factor_qualification_sha256')]
            resumed = validated_frozen_resume(args.resume_frozen_validation, checkpoint,
                                               args.replay_mode, receipts)
            executable_path = str(args.exe.resolve())
            if resumed['frozen_sha256'].get(executable_path) != json.loads(
                    small_path.read_text()).get('exe_sha256'):
                raise ValueError('resume executable differs from frozen small qualification')
            e['resumed_frozen_validation'] = str(args.resume_frozen_validation)
            e['resumed_frozen_validation_sha256'] = digest(args.resume_frozen_validation)
            e['qualification_source_scope'] = 'previous verified frozen executable; current source not requalified'
        small = json.loads(small_path.read_text())
        if (small.get('terminal') is not True or small.get('all_checks_passed') is not True
                or small.get('fresh_keys_generated') is not False or small.get('root_evaluations_after') != 22
                or small.get('replay_actual_exit_code') != 0 or small.get('build_actual_exit_code') != 0
                or set(small.get('ciphertext_identity', {})) != {'sigmoid', 'gelu', 'tanh'}
                or not all(small['ciphertext_identity'].values())
                or not all(small.get('independent_phase_checks', {}).get(f, {}).get(
                    'independent_phase_accuracy_passed') is True for f in ('sigmoid', 'gelu', 'tanh'))):
            raise ValueError('small root-reuse qualification incomplete')
        e['small_qualification_sha256'] = digest(small_path)
        if selector56 and (not small.get('results') or not all(
                row.get('selector56_replay') is True and
                row.get('selector_matrix_digit_bits') == 56 and
                row.get('selector_matrix_digits') == 10 for row in small['results'])):
            raise ValueError('small qualification did not exercise configured selector56 replay')
        if secondary11:
            if not small.get('results') or not all(
                    row.get('secondary11_replay') is True and
                    row.get('derived_tree_dd_secondary_bits') == 11
                    for row in small['results']):
                raise ValueError('small qualification did not exercise secondary11 replay')
            dd = json.loads(args.dd_cache_validation.read_text())
            if (dd.get('terminal') is not True or dd.get('all_component_checks_passed') is not True
                    or dd.get('validation_controller_actual_exit_code') != 0
                    or dd.get('extracted_public_key_files_unchanged') is not True
                    or {level.get('q') for level in dd.get('levels', [])} !=
                       {479, 432, 385, 338, 291, 244, 197}
                    or not all(level.get('run_actual_exit_code') == 0 and
                               level.get('report', {}).get('all_coefficients_equal') is True
                               for level in dd.get('levels', []))):
                raise ValueError('retained secondary11 DD qualification incomplete')
            for path, expected in dd['source_sha256'].items():
                if resumed is None and digest(repo/path) != expected:
                    raise ValueError('source changed after retained DD qualification: '+path)
            e['dd_cache_qualification_sha256'] = digest(args.dd_cache_validation)
        e['shared_input_tiles'] = small.get('shared_input_tiles') is True
        e['shared_basis_union'] = small.get('shared_basis_union') is True
        if e['shared_basis_union']:
            if (not e['shared_input_tiles'] or small.get('new_roots_by_function') != [22, 0, 0]
                    or small.get('reused_roots_by_function') != [0, 8, 18]
                    or small.get('small_selector_groups_after') != 35):
                raise ValueError('small frequency-union qualification incomplete')
            if args.factor_validation is None:
                raise ValueError('frequency-union replay requires completed full22 factor qualification')
            factor = json.loads(args.factor_validation.read_text())
            if (factor.get('terminal') is not True or factor.get('all_checks_passed') is not True
                    or factor.get('actual_exit_code') != 0 or factor.get('columns') != 22
                    or factor.get('slots_checked') != 720896
                    or factor.get('fresh_keys_generated') is not False):
                raise ValueError('full22 factor qualification incomplete')
            for path, expected in factor['qualified_source_sha256'].items():
                if resumed is None and digest(repo/path) != expected:
                    raise ValueError('source changed after full22 factor qualification: '+path)
            if selector56:
                trials = [row for row in factor.get('completed_rows', [])
                          if 'ciphertexts_exact' in row]
                if (len(trials) != 4 or not all(
                        row.get('ciphertexts_exact') is True and row.get('n') == 32768 and
                        row.get('columns') == 22 and row.get('selector_chunk_size') == 40 and
                        row.get('schedule_configured_width') == (row.get('matrix_digit_bits') == 56)
                        for row in trials) or
                        {(row.get('tail'), row.get('matrix_digit_bits'), row.get('selector_digits'))
                         for row in trials} !=
                        {(False, 48, 11), (False, 56, 10), (True, 48, 11), (True, 56, 10)}):
                    raise ValueError('full22 factor qualification did not exercise configured selector56')
            e['factor_qualification_sha256'] = digest(args.factor_validation)
        e['source_sha256'] = small['source_sha256']
        if resumed is not None and e['source_sha256'] != resumed.get('source_sha256'):
            raise ValueError('resume source provenance differs from previous qualification')
        for path, expected in small['source_sha256'].items():
            if resumed is None and digest(repo/path) != expected:
                raise ValueError('source changed after small qualification: '+path)
        executable = (args.exe or repo/'TFHEpp/build/test/ckks/ckks_narrow_benchmark').resolve()
        artifacts = [(executable, small['exe_sha256'])]
        artifacts.append((small_path.resolve(), e['small_qualification_sha256']))
        if secondary11:
            artifacts.append((args.dd_cache_validation.resolve(), e['dd_cache_qualification_sha256']))
        if e['shared_basis_union']:
            artifacts.append((args.factor_validation.resolve(), e['factor_qualification_sha256']))
        replay_executable = frozen/executable.name
        e['executable_source_path'] = str(executable)
        e['small_qualification_frozen_path'] = str(frozen/small_path.name)
        e['factor_stage_profiled'] = bool(small.get('results')) and all(
            isinstance(row.get('factor_stage_seconds'), dict) for row in small['results'])
        e['root_tile_stage_profiled'] = bool(small.get('results')) and all(
            isinstance(row.get('root_tile_stage_seconds'), dict) for row in small['results'])
        for name in ('narrow_ckks_basis_checkpoint_validate.py', 'narrow_ckks_prefix_arithmetic.py'):
            if resumed is None:
                source = Path(__file__).parent/name
                expected = digest(source)
            else:
                matches = [(Path(path), value) for path, value in resumed['frozen_sha256'].items()
                           if Path(path).name == name]
                if len(matches) != 1:
                    raise ValueError('resume decoder artifact missing: '+name)
                source, expected = matches[0]
            artifacts.append((source, expected))
        controller_source = Path(__file__).resolve()
        artifacts.append((controller_source, digest(controller_source)))
        e['frozen_sha256'] = {}
        for source, expected in artifacts:
            target = frozen/source.name
            shutil.copy2(source, target)
            if digest(target) != expected:
                raise ValueError('artifact differs from qualification: '+str(source))
            e['frozen_sha256'][str(target)] = expected
        own_fields = (Path('/proc')/str(os.getpid())/'stat').read_text().rsplit(')', 1)[1].split()
        e['controller'] = dict(pid=os.getpid(), start_ticks=int(own_fields[19]))
        e['controller_sha256'] = digest(Path(__file__))
        e['controller_frozen_path'] = str(frozen/controller_source.name)
        resource.setrlimit(resource.RLIMIT_AS, (128 << 20, 40*GIB))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        while True:
            reference = json.loads(args.reference_evidence.read_text())
            if reference.get('terminal') is True and reference.get('all_checks_passed') is not True:
                raise ValueError('original basis validation did not pass')
            preceding = {}
            for path in args.wait_terminal:
                previous = json.loads(path.read_text())
                active = previous.get('active_child')
                preceding[str(path)] = dict(
                    idle=previous.get('terminal') is True and
                         (not active or not process(active['pid'], active['start_ticks'])),
                    accuracy_passed=previous.get('all_checks_passed'))
            preceding_idle = all(job['idle'] for job in preceding.values())
            # The preceding test need only finish to release its reservation;
            # its optional-cache accuracy result does not gate this independent test.
            packed = json.loads(args.packed_evidence.read_text())
            packed_job = process(pid, start)
            remaining = (max(0, 64*GIB-packed_job['rss_bytes']) if packed_job else
                         0 if packed.get('terminal') is True and packed.get('all_passed') is True else 64*GIB)
            available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                                 if line.startswith('MemAvailable:')))*1024
            required = remaining+50*GIB+(128 << 20)
            disk_free = shutil.disk_usage(root).free
            basis_gate = reference.get('terminal') is True and reference.get('all_checks_passed') is True
            e.update(stage='waiting for original validation, preceding tests and memory headroom',
                     original_basis_validated=basis_gate, preceding_test_terminal=preceding_idle,
                     preceding_tests=preceding, packed_process=packed_job,
                     budget=dict(available_bytes=available, required_bytes=required,
                                 packed_remaining_reservation_bytes=remaining,
                                 free_disk_bytes=disk_free, required_disk_bytes=2*GIB))
            save()
            if basis_gate and preceding_idle and available >= required and disk_free >= 2*GIB:
                break
            print(e['stage'], flush=True)
            time.sleep(30)
        e['reference_qualification_sha256'] = digest(args.reference_evidence)
        baseline_path = args.reference_root/'retained-baseline.json'
        baseline = json.loads(baseline_path.read_text())
        e['reference_baseline_sha256'] = digest(baseline_path)
        e['stage'] = 'checking retained payload baseline'
        save()
        originals = sorted(p for p in checkpoint.rglob('*') if p.is_file())
        if args.reuse_validated_baseline:
            e['original_sha256'] = validated_payload_baseline(reference, checkpoint, baseline, originals)
            e['baseline_source'] = 'completed prerequisite full payload validation'
            e['payload_rehashed_before_evaluation'] = False
        else:
            e['original_sha256'] = {str(p): dict(bytes=p.stat().st_size, sha256=digest(p)) for p in originals}
            e['baseline_source'] = 'fresh full payload scan'
            e['payload_rehashed_before_evaluation'] = True
        if e['original_sha256'] != baseline:
            raise ValueError('checkpoint changed since original full validation')
        functions = ('sigmoid', 'gelu', 'tanh')
        archived = {f: digest(args.reference_root/(f+'-replay.bin')) for f in functions}
        e['reference_result_sha256'] = archived
        run('replay', ['/usr/bin/time', '-v', str(replay_executable),
                       args.replay_mode, 'all', str(checkpoint), str(root/'replay.json')], 40*GIB)
        reported = json.loads((root/'replay.json').read_text())
        rows = reported.get('results', [])
        if (reported.get('accuracy_passed') is not True or reported.get('setup_generation_count') != 0
                or len(rows) != 3 or {r.get('function') for r in rows} != set(functions)):
            raise ValueError('full root-reuse report gates failed')
        e['report'] = reported
        e['ciphertext_identity'] = {}
        e['independent_phase_checks'] = {}
        root_counts = (('sigmoid', 22, 0), ('gelu', 0, 8), ('tanh', 0, 18)) if e['shared_basis_union'] else (
            ('sigmoid', 8, 0), ('gelu', 0, 8), ('tanh', 14, 4))
        for function, new, reused in root_counts:
            row = next(r for r in rows if r['function'] == function)
            if selector56 and (row.get('selector56_replay') is not True or
                               row.get('selector_matrix_digit_bits') != 56 or
                               row.get('selector_matrix_digits') != 10):
                raise ValueError('full configured selector56 policy mismatch: '+function)
            if secondary11 and (row.get('secondary11_replay') is not True or
                                row.get('derived_tree_dd_secondary_bits') != 11 or
                                row.get('secondary_digit_bits') != 8):
                raise ValueError('full secondary11 cache policy mismatch: '+function)
            if (row.get('n') != 32768 or row.get('fresh_keys_generated') is not False
                    or row.get('input_reused') is not True or row.get('shared_frequency_roots') is not True
                    or row.get('retained_square_key_shared') is not True or row.get('accuracy_passed') is not True
                    or row.get('newly_evaluated_roots') != new or row.get('reused_roots') != reused):
                raise ValueError('full root-reuse profile or count mismatch: '+function)
            if e['shared_input_tiles'] and row.get('factor_groups_completed') != (0 if new == 0 else 36):
                raise ValueError('full joint-input selector group count mismatch: '+function)
            if e['shared_basis_union']:
                scope = ('shared all-function root generation and this function reconstruction' if new else
                         'this function reconstruction from shared roots')
                if (row.get('shared_basis_union') is not True or row.get('series_tile_size') != 22
                        or row.get('refresh_timing_scope') != scope):
                    raise ValueError('full frequency-union profile or timing scope mismatch: '+function)
            if e['factor_stage_profiled']:
                stages = row.get('factor_stage_seconds', {})
                expected_stages = {'selector_load', 'table_ntt', 'selector_ntt',
                                   'matrix', 'inverse_crt', 'conjugation'}
                if (set(stages) != expected_stages or row.get('root_stage_timings_available') is not True
                        or not all(isinstance(value, (float, int)) and math.isfinite(value) and value >= 0
                                   for value in stages.values())):
                    raise ValueError('invalid full factor-stage timings: '+function)
                root_seconds = row.get('root_generation_seconds', -1)
                reconstruction_seconds = row.get('basis_reconstruction_seconds', -1)
                refresh_seconds = row.get('refresh_seconds', -1)
                if (not all(math.isfinite(value) and value >= 0 for value in
                            (root_seconds, reconstruction_seconds, refresh_seconds))
                        or sum(stages.values()) > root_seconds + 1e-3
                        or root_seconds + reconstruction_seconds > refresh_seconds + 1e-3):
                    raise ValueError('full stage timing containment failed: '+function)
            if e['root_tile_stage_profiled']:
                validate_root_tile_timings(row)
            result = root/(function+'-replay.bin')
            e['ciphertext_identity'][function] = digest(result) == archived[function]
            if not e['ciphertext_identity'][function]:
                raise ValueError('full retained ciphertext differs: '+function)
            phase = root/(function+'-independent.json')
            run(function+'-independent', [sys.executable, str(frozen/'narrow_ckks_basis_checkpoint_validate.py'),
                '--checkpoint', str(checkpoint/'checkpoint.bin'), '--result', str(result),
                '--report', str(root/'replay.json'), '--function', function, '--output', str(phase)], 2*GIB)
            e['independent_phase_checks'][function] = json.loads(phase.read_text())
            if e['independent_phase_checks'][function]['independent_phase_accuracy_passed'] is not True:
                raise ValueError('independent basis phase check failed: '+function)
        if any(digest(args.reference_root/(f+'-replay.bin')) != archived[f] for f in functions):
            raise ValueError('reference ciphertext file changed')
        verify_frozen()
        e['all_checks_passed'] = True
    except Exception as error:
        e['error'] = str(error)
        e['all_checks_passed'] = False
    finally:
        # Output checks alone do not qualify the retained-key replay. Keep the
        # public success flag false while the final payload hashes are pending.
        e['evaluation_checks_passed'] = e['all_checks_passed']
        e['all_checks_passed'] = False
        if originals is not None and 'original_sha256' in e:
            try:
                e['stage'] = 'checking retained files after evaluation'
                save()
                after = {str(p): dict(bytes=p.stat().st_size, sha256=digest(p)) for p in originals}
                e['payload_rehashed_after_evaluation'] = True
                e['checkpoint_files_unchanged'] = after == e['original_sha256']
                e['checkpoint_files_checked'] = len(originals)
            except Exception as error:
                e['checkpoint_integrity_error'] = str(error)
                e['checkpoint_files_unchanged'] = False
        e['all_checks_passed'] = (
            e['evaluation_checks_passed'] and
            e.get('checkpoint_files_unchanged') is True)
        e['full_encrypted_root_reuse_validated'] = e['all_checks_passed']
        e['terminal'] = True
        e['stage'] = 'complete'
        save()
    return 0 if e['all_checks_passed'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
