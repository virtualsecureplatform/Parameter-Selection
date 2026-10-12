#!/usr/bin/env python3
"""Run one retained packed setup (or reuse one), replay it, and check saved phases.

Private keys/selectors remain in the checkpoint directory. This diagnostic does
not qualify security or independently inspect all selector/key coefficients.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess

from narrow_ckks_packed_checkpoint_validate import validate
from narrow_ckks_packed_cpu_compare import digest


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def process_identity(pid):
    try:
        fields = Path('/proc/'+str(pid)+'/stat').read_text().rsplit(')',1)[1].split()
        return dict(pid=pid,start_ticks=int(fields[19]))
    except FileNotFoundError:
        return dict(pid=pid,start_ticks=None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--n',type=int,choices=[512,16384],required=True)
    parser.add_argument('--checkpoint',type=Path,help='reuse a completed compatible checkpoint')
    args = parser.parse_args()
    args.root.mkdir(mode=0o700)
    root = args.root.resolve()
    checkpoint = args.checkpoint.resolve() if args.checkpoint else root/'checkpoint'
    exe = root/'ckks_narrow_benchmark'
    shutil.copy2(args.exe,exe)
    if digest(exe) != digest(args.exe):
        raise ValueError('pinned benchmark copy mismatch')
    state = dict(scope=__doc__,utc=utc(),n=args.n,root=str(root),checkpoint=str(checkpoint),
                 fresh_setup_required=args.checkpoint is None,terminal=False,all_passed=False,
                 runner=process_identity(os.getpid()),exe_sha256=digest(exe),runs=[],
                 omp_threads=4,address_space_limit_bytes=(4 if args.n==512 else 64)*1024**3,
                 security_validated=False,goal_complete=False,
                 source_sha256={p.name:digest(p) for p in (Path(__file__),
                     Path(__file__).with_name("narrow_ckks_packed_checkpoint_validate.py"),
                     Path(__file__).with_name("narrow_ckks_prefix_arithmetic.py"))})
    def save():
        state['updated_utc']=utc()
        pending = args.output.with_suffix(args.output.suffix+'.pending')
        pending.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n')
        pending.replace(args.output)
    save()
    def limits():
        bound = state['address_space_limit_bytes']
        resource.setrlimit(resource.RLIMIT_AS,(bound,bound))
        resource.setrlimit(resource.RLIMIT_CORE,(0,0))
        resource.setrlimit(resource.RLIMIT_STACK,(64*1024**2,64*1024**2))
    def files():
        return {str(p.relative_to(checkpoint)):digest(p) for p in sorted(checkpoint.rglob('*')) if p.is_file()}
    def run(mode,label):
        report = root/(label+'.json');log = root/(label+'-time.log')
        command = ['/usr/bin/time','-v',str(exe),mode,str(checkpoint),str(report)]
        with log.open('wb') as out:
            process = subprocess.Popen(command,stdout=out,stderr=subprocess.STDOUT,
                env=dict(os.environ,OMP_NUM_THREADS='4',OMP_STACKSIZE='64M'),preexec_fn=limits)
            state['active_run']=dict(label=label,command=command,time_process=process_identity(process.pid))
            save()
            code=process.wait()
        state.pop('active_run',None)
        entry=dict(label=label,actual_exit_code=code,report=str(report),time_log=str(log))
        state['runs'].append(entry);save()
        if code!=0:
            raise ValueError('benchmark failed: '+label)
        evidence=validate(checkpoint/'checkpoint.bin',Path(str(report)+'.ciphertext'),
                          Path(str(report)+'.square.ciphertext'),report)
        (root/(label+'-phase-validation.json')).write_text(json.dumps(evidence,indent=2)+'\n')
        entry['independent_phase_validation']=evidence;save()
        return report
    try:
        stem = '--small-high-packed-bounded-' if args.n==512 else '--ll14-packed-bounded-'
        first = run(stem+('replay' if args.checkpoint else 'checkpoint'),'retained')
        state['stage']='hashing completed retained checkpoint';save()
        baseline=files()
        state['checkpoint_file_count']=len(baseline)
        state['checkpoint_hashes']=baseline
        state['stage']='replaying retained setup';save()
        second = run(stem+'replay','replay')
        report=json.loads(second.read_text())
        if report.get('refresh_keys_reused') is not True or report.get('post_refresh_square_key_reused') is not True:
            raise ValueError('replay regenerated evaluation keys')
        if report['packed_setup_seconds']!=0 or report['selector_reconstruction_seconds']!=0:
            raise ValueError('replay repeated packed setup')
        state['stage']='checking retained file and ciphertext identity';save()
        if baseline!=files():
            raise ValueError('retained checkpoint files changed across replay')
        for suffix in ('.ciphertext','.square.ciphertext'):
            if digest(Path(str(first)+suffix))!=digest(Path(str(second)+suffix)):
                raise ValueError('creation/replay ciphertext mismatch: '+suffix)
        state.update(terminal=True,all_passed=True,stage='completed',
                     checkpoint_files_unchanged=True,refresh_and_square_ciphertexts_identical=True,
                     refresh_and_square_keys_reused=True,full_ring_independent_phase_validated=args.n==16384)
    except BaseException as error:
        state.update(terminal=True,all_passed=False,stage='failed',error=repr(error));save();raise
    save();print('Retained packed setup/replay and independent phase checks passed at N'+str(args.n))


if __name__=='__main__':
    main()
