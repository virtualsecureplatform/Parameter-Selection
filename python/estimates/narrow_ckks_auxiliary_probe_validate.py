#!/usr/bin/env python3
"""Check auxiliary8 encrypted probes against fixed 23-bit gates.

Small-ring diagnostics do not establish full-ring accuracy or key security.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def validate(path, function, n=512):
    r = json.loads(path.read_text())
    if n not in (512, 16384):
        raise ValueError('unsupported probe ring dimension')
    expected = dict(profile='auxiliary8-small-probe' if n == 512 else 'auxiliary8-full-probe', n=n, h=31,
                    function=function, auxiliary_bits=8, gaussian_sigma=3.2,
                    message_bound=4, target_precision_bits=23,
                    injected_failure=False)
    if n == 16384:
        expected.update(source_log_q=48, selector_log_q=385, structured_switch_key_log_q=56,
                        factor_key_log_q=356, s2c_input_key_log_q=93, output_conjugation_key_log_q=116,
                        input_log_q=85, input_log_delta=44, output_log_q=108, output_log_delta=48,
                        root_plain_log_delta=37, s2c_plain_log_delta=37, storage_bits=448,
                        series_tile_size=12, degree=40 if function == 'tanh' else 12,
                        compact_dd=True, compact_s2c=True, relinearize_before_rescale=True,
                        packed_selectors=False, file_backed_selectors=True,
                        input_pattern='random-with-endpoints')
    for key, value in expected.items():
        if type(r.get(key)) is not type(value) or r[key] != value:
            raise ValueError('probe profile mismatch: ' + key)
    if type(r.get('accuracy_passed')) is not bool:
        raise ValueError('accuracy flag is not boolean')
    fields = ('max_error', 'evaluation_max_error', 'fresh_input_max_error',
              'rms_error', 'complex_uniform_bound', 'post_refresh_product_max_error',
              'post_refresh_product_evaluation_max_error', 'refresh_seconds')
    for key in fields:
        if type(r.get(key)) not in (int, float) or not math.isfinite(r[key]) or r[key] < 0:
            raise ValueError('invalid probe metric: ' + key)
    tolerance = 2**-23
    endpoint = (1/(1+math.exp(-4)) if function == 'sigmoid' else
                2*(1+math.erf(4/math.sqrt(2))) if function == 'gelu' else math.tanh(4))
    amplitude = math.sqrt(2) * endpoint
    product_tolerance = 2*amplitude*tolerance+tolerance*tolerance+tolerance
    evaluation_product_tolerance = (2*max(0, amplitude-r['complex_uniform_bound'])*tolerance
                                    + tolerance*tolerance+tolerance)
    passed = (all(r[key] < tolerance for key in fields[:5]) and
              r['post_refresh_product_max_error'] < product_tolerance and
              r['post_refresh_product_evaluation_max_error'] < evaluation_product_tolerance)
    if r['accuracy_passed'] != passed:
        raise ValueError('accuracy flag disagrees with independent gates')
    return dict(report=str(path), report_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                accuracy_passed=passed, **{key: r[key] for key in fields},
                product_tolerance=product_tolerance,
                evaluation_product_tolerance_lower_bound=evaluation_product_tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--n', type=int, choices=(512, 16384), default=512)
    args = parser.parse_args()
    checked = {name: validate(args.root/(name+'.json'), name, args.n)
               for name in ('sigmoid', 'gelu', 'tanh')}
    passed = all(r['accuracy_passed'] for r in checked.values())
    result = dict(scope=__doc__, n=args.n, small_encrypted_validation_passed=passed if args.n == 512 else False,
                  full_encrypted_validation_passed=passed if args.n == 16384 else False, security_validated=False,
                  functions=checked)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('N', args.n, 'auxiliary8 accuracy', {name: r['accuracy_passed'] for name, r in checked.items()})
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
