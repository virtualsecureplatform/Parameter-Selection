#!/usr/bin/env python3
"""Check completed packed identity reports against the unchanged C++ gates.

Checks reported metrics, profile, terminal /usr/bin/time output and cache cleanup.
Does not independently decrypt ciphertexts, qualify security, or compare speed
across profiles. Full-ring results cannot be substituted by small-ring controls.
"""
import argparse
import datetime
import json
from pathlib import Path
import re

from narrow_ckks_packed_cpu_compare import checked_report, digest


def validate(report_path, time_log, selector_directory, n=16384):
    report = json.loads(report_path.read_text())
    profile = checked_report(report, n)
    time = time_log.read_text()
    if ('Command terminated by signal' in time or
            not re.search(r'Exit status:\s*0\s*$', time)):
        raise ValueError('packed run did not complete with exit0')
    if selector_directory.exists() or selector_directory.is_symlink():
        raise ValueError('owned selector directory was not cleaned up')
    timings = {key: report[key] for key in
               ['packed_setup_seconds', 'selector_reconstruction_seconds',
                'key_setup_seconds', 'refresh_seconds']}
    if any(type(value) not in (float, int) or value <= 0 for value in timings.values()):
        raise ValueError('invalid packed phase timing')
    storage = {key: report[key] for key in
               ['packed_payload_bytes', 'packed_serialized_key_bytes',
                'packed_substitution_cache_bytes', 'dd_chain_storage_bytes', 'ntt_cache_bytes']}
    if any(type(value) is not int or value < 0 for value in storage.values()):
        raise ValueError('invalid packed storage measurement')
    if storage['packed_payload_bytes'] != 4 * 2 * n * 7 * 8:
        raise ValueError('packed indicator payload geometry mismatch')
    if storage['packed_serialized_key_bytes'] <= storage['packed_payload_bytes']:
        raise ValueError('packed key archive does not include substitution keys')
    return dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), scope=__doc__,
                terminal=True, actual_exit_code=0, reported_accuracy_validated=True,
                full_packed_report_validated=n == 16384, profile=profile,
                report=str(report_path), report_sha256=digest(report_path),
                time_log=str(time_log), time_log_sha256=digest(time_log),
                selector_directory=str(selector_directory), owned_cache_cleanup_passed=True,
                phase_seconds=timings, storage_bytes=storage,
                max_error=report['max_error'], rms_error=report['rms_error'],
                selector_reconstruction_max_error=report['selector_reconstruction_max_error'],
                post_refresh_product_max_error=report['post_refresh_product_max_error'],
                independent_phase_validated=False, security_validated=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['report', 'time-log', 'selector-directory', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--n', type=int, choices=[512, 16384], default=16384)
    a = parser.parse_args()
    result = validate(a.report, a.time_log, a.selector_directory, a.n)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Packed profile, fixed reported gates, terminal exit and cleanup validated at N'+str(a.n))


if __name__ == '__main__':
    main()
