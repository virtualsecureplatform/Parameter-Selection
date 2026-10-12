#!/usr/bin/env python3
"""Run and retain a full N32768 lean frequency-basis checkpoint once.

Requires current small checkpoint/replay and wide arithmetic evidence, checks
live-job memory reservations and disk space, and can wait for prerequisite
validations. Replays all functions with unchanged keys/input, compares exact
ciphertext bytes and independently decodes each refresh/square/normalization.
This does not qualify structured-key security or reproduce paper timings.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time

from narrow_ckks_corrected_selector_replay import GIB, digest, process
from narrow_ckks_basis_checkpoint_validate import validate


def save(path, data):
    pending = path.with_suffix(path.suffix+'.pending')
    pending.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    pending.replace(path)


def budget(reservations, cap, root):
    jobs = []
    for pid, start, limit, completion in reservations:
        job = process(pid, start, limit)
        if job is not None:
            jobs.append(job)
        elif completion is None or not completed(completion):
            # Retain the entire cap while an evaluator's parent validates or
            # launches a follow-up child; an absent PID alone is insufficient.
            jobs.append(dict(pid=pid,start_ticks=start,rss_bytes=0,cap_bytes=limit,
                             scope='reserved until controller completion'))
    available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                         if line.startswith('MemAvailable:')))*1024
    required = sum(max(0, job['cap_bytes']-job['rss_bytes']) for job in jobs)+cap+10*GIB+512*2**20
    free_disk = shutil.disk_usage(root).free
    # Selector payload179.7GB, setup/input archive and all public outputs,
    # plus filesystem margin. Recheck before each child.
    required_disk = 210*GIB
    return dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), jobs=jobs,
                available_memory_bytes=available, required_memory_bytes=required,
                free_disk_bytes=free_disk, required_disk_bytes=required_disk,
                child_cap_bytes=cap, parent_cap_bytes=512*2**20,
                eligible=available >= required and free_disk >= required_disk)


def completed(path):
    if not path.exists():
        return False
    data = json.loads(path.read_text())
    if data.get('terminal') is not True:
        return False
    passed = data.get('all_checks_passed') is True or data.get('all_requested_functions_passed') is True
    if passed:
        return True
    if (data.get('error') or 'failed' in data.get('stage','') or
            data.get('actual_exit_code') not in (None,0) or data.get('stage') == 'complete'):
        raise ValueError('prerequisite validation failed: '+str(path))
    # Some controllers mark the evaluator terminal before independent phase
    # checks finish. Wait for their positive gates rather than treating that
    # intermediate status as either success or failure.
    return False


def run(args):
    os.umask(0o077)
    resource.setrlimit(resource.RLIMIT_AS, (512*2**20, resource.getrlimit(resource.RLIMIT_AS)[1]))
    small = json.loads(args.small_validation.read_text())
    wide = json.loads(args.wide_validation.read_text())
    if (small.get('terminal') is not True or small.get('all_checks_passed') is not True or
            small.get('corrected_build_actual_exit_code') != 0 or
            small.get('corrected_replay_actual_exit_code') != 0 or
            small.get('replay_fresh_keys_generated') is not False or
            set(small.get('ciphertext_identity', {})) != {'sigmoid','gelu','tanh'} or
            not all(small['ciphertext_identity'].values()) or
            not all(small.get('independent_phase_checks', {}).get(f, {}).get('independent_phase_accuracy_passed') is True
                    for f in ['sigmoid','gelu','tanh'])):
        raise ValueError('small basis checkpoint provenance incomplete')
    if (wide.get('all_checks_passed') is not True or wide.get('actual_exit_code') != 0 or
            wide.get('coefficient_checks') != 1179648):
        raise ValueError('basis storage arithmetic evidence incomplete')
    source = args.repo/'test/ckks/ckks_narrow_benchmark.cpp'
    executable = args.repo/'build/test/ckks/ckks_narrow_benchmark'
    if digest(source) != small['source_sha256'] or digest(executable) != small['exe_sha256']:
        raise ValueError('basis source/executable changed after small validation')
    if args.destination.exists():
        raise ValueError('basis run destination must be new; existing checkpoints must be replayed directly')
    args.destination.mkdir(mode=0o700)
    frozen = args.destination/'ckks_narrow_benchmark'
    shutil.copy2(executable, frozen)
    result = dict(scope=__doc__, terminal=False, stage='checking prerequisites', all_checks_passed=False,
                  source_sha256=small['source_sha256'], executable_sha256=digest(frozen),
                  destination=str(args.destination), initial_fixture_generation_count=0,
                  replay_fresh_keys_generated=False, full_encrypted_basis_validation_passed=False,
                  security_validated=False, paper_performance_reproduced=False)
    save(args.output, result)
    reservations = []
    for item in args.reserved_job:
        parts = item.split(':',3)
        pid, start, cap = map(int,parts[:3])
        completion = Path(parts[3]) if len(parts)==4 else None
        reservations.append((pid,start,cap*GIB,completion))
    cap = 40*GIB
    while True:
        gates = {str(path):completed(path) for path in args.wait_validation}
        memory = budget(reservations, cap, args.destination)
        result.update(prerequisites=gates, memory_budget=memory)
        if all(gates.values()) and memory['eligible']:
            break
        if not args.wait:
            raise ValueError('basis launch prerequisites or budget not ready')
        result['stage'] = 'waiting for validated completion and memory headroom'
        save(args.output, result)
        print(result['stage'], json.dumps({'prerequisites':gates,'budget_eligible':memory['eligible']}), flush=True)
        time.sleep(30)
    checkpoint = args.destination/'checkpoint'
    baseline = None
    for phase in ['initial','replay']:
        # Protect against intervening builds and ensure the frozen executable
        # retains exactly the binary which passed small validation.
        if digest(frozen) != small['exe_sha256']:
            raise ValueError('frozen basis executable changed')
        memory = budget(reservations, cap, args.destination)
        # After creation, retained selectors consume disk space. The initial
        # reservation is no longer needed for a read-only replay.
        if phase == 'replay':
            memory['required_disk_bytes'] = 2*GIB
            memory['eligible'] = (memory['available_memory_bytes'] >= memory['required_memory_bytes']
                                  and memory['free_disk_bytes'] >= 2*GIB)
        if not memory['eligible']:
            raise ValueError('basis child budget changed before '+phase)
        report = args.destination/(phase+'.json')
        def limits():
            resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
            resource.setrlimit(resource.RLIMIT_STACK, (64*2**20,64*2**20))
            resource.setrlimit(resource.RLIMIT_CORE, (0,0))
        mode = '--lean-basis-fbs-'+('checkpoint' if phase == 'initial' else 'replay')
        with (args.destination/(phase+'.log')).open('x') as log:
            child = subprocess.Popen(['/usr/bin/time','-v','-o',str(args.destination/(phase+'.time')),
                str(frozen),mode,'all',str(checkpoint),str(report)],stdout=log,stderr=subprocess.STDOUT,
                preexec_fn=limits,env=dict(os.environ,OMP_NUM_THREADS='4',OMP_STACKSIZE='64M'))
            result.update(stage=phase+' full basis running', child_pid=child.pid, active_phase=phase,
                          memory_budget=memory)
            if phase == 'initial':result['initial_fixture_generation_count'] = 1
            save(args.output,result)
            print(result['stage'], child.pid, flush=True)
            code = child.wait()
        result[phase+'_actual_exit_code'] = code
        save(args.output,result)
        if code != 0:
            raise ValueError('full basis evaluator failed: '+phase)
        reported = json.loads(report.read_text())
        rows = reported.get('results', [])
        if (reported.get('accuracy_passed') is not True or len(rows) != 3 or
                {row.get('function') for row in rows} != {'sigmoid','gelu','tanh'} or
                reported.get('setup_generation_count') != (1 if phase == 'initial' else 0) or
                not all(row.get('n') == 32768 and row.get('accuracy_passed') is True and
                        row.get('fresh_keys_generated') is (phase == 'initial') and
                        row.get('retained_square_key_shared') is True for row in rows)):
            raise ValueError('full basis report gates failed: '+phase)
        result[phase+'_report'] = reported
        result[phase+'_independent'] = {}
        for function in ['sigmoid','gelu','tanh']:
            output = args.destination/(function+('-initial.bin' if phase == 'initial' else '-replay.bin'))
            checked = validate(checkpoint/'checkpoint.bin',output,report,function)
            if not checked['independent_phase_accuracy_passed']:
                raise ValueError('independent full basis phase failed: '+phase+' '+function)
            result[phase+'_independent'][function] = checked
            save(args.output,result)
        files = sorted(path for path in checkpoint.rglob('*') if path.is_file())
        hashes = {str(path):dict(bytes=path.stat().st_size,sha256=digest(path)) for path in files}
        if phase == 'initial':
            baseline = hashes
            save(args.destination/'retained-baseline.json',baseline)
        elif hashes != baseline:
            raise ValueError('basis checkpoint or selector payload changed during replay')
    identical = {f:digest(args.destination/(f+'-initial.bin')) == digest(args.destination/(f+'-replay.bin'))
                 for f in ['sigmoid','gelu','tanh']}
    if not all(identical.values()):
        raise ValueError('basis replay ciphertexts differ')
    result.update(terminal=True, stage='complete', all_checks_passed=True,
                  full_encrypted_basis_validation_passed=True, checkpoint_files_unchanged=True,
                  checkpoint_files_checked=len(baseline), ciphertext_identity=identical)
    save(args.output,result)
    print('full basis initial and retained replay independently passed',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['repo','destination','small-validation','wide-validation','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--reserved-job',action='append',default=[],help='PID:start_ticks:cap_GiB[:controller_validation.json]')
    parser.add_argument('--wait-validation',type=Path,action='append',default=[])
    parser.add_argument('--wait',action='store_true')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('basis controller output must be new; preserve existing run evidence')
    try:
        run(args)
    except Exception as error:
        data = json.loads(args.output.read_text()) if args.output.exists() else {}
        data.update(terminal=True,all_checks_passed=False,stage='basis controller failed',error=str(error))
        save(args.output,data)
        raise


if __name__ == '__main__':
    main()
