#!/usr/bin/env python3
"""Validate retained aux8 function outputs with independent integer phases.

Also compares saved secrets, slots and input bytes to the original state.
Approximation/evaluation bounds remain C++ reported metrics; security and
the selector bit error's origin are outside this diagnostic.
"""
import argparse
import json
import math
from pathlib import Path
import struct

from narrow_ckks_auxiliary_probe_prefix import inspect as inspect_prefix
from narrow_ckks_auxiliary_square_validate import inspect as inspect_square


def profile_gate(report):
    function = report.get('function')
    if (function not in ('sigmoid', 'gelu', 'tanh') or report.get('n') not in (512, 16384)
            or report.get('auxiliary_bits') != 8 or report.get('gaussian_sigma') != 3.2
            or report.get('degree') != (40 if function == 'tanh' else 12)
            or report.get('series_tile_size') != 12 or report.get('fresh_keys_generated') is not False
            or report.get('input_reused') is not True or report.get('retained_keys_reused') is not True):
        raise ValueError('retained-function report profile mismatch')
    names = ['fresh_input_max_error', 'max_error', 'rms_error', 'evaluation_max_error',
             'complex_uniform_bound', 'target_amplitude', 'post_refresh_product_max_error',
             'square_rms_error', 'post_refresh_product_evaluation_max_error',
             'product_tolerance', 'evaluation_product_tolerance']
    if any(type(report.get(k)) not in (int, float) or not math.isfinite(report[k]) or report[k] < 0
           for k in names):
        raise ValueError('missing, nonfinite or negative retained-function metric')
    if report['product_tolerance'] <= 0 or report['evaluation_product_tolerance'] <= 0:
        raise ValueError('invalid retained-function product bound')
    return function


def same_retained_input(state, prefix, function, n):
    def read(file, count):
        value = file.read(count)
        if len(value) != count:
            raise ValueError('truncated retained-function identity input')
        return value
    with state.open('rb') as original, prefix.open('rb') as current:
        left, right = read(original, 81), read(current, 81)
        if left != right or left[0] not in (0, 1):
            raise ValueError('retained-function profile bytes differ')
        endian = '<' if left[0] else '>'
        for file, expected in [(original, 'sigmoid'), (current, function)]:
            length = struct.unpack(endian+'Q', read(file, 8))[0]
            if length != len(expected) or read(file, length) != expected.encode() or read(file, 1) != b'\x01':
                raise ValueError('retained-function tag or completion flag mismatch')
        # Two seven-limb secrets, complex target slots, two input polynomials.
        if read(original, n*232) != read(current, n*232):
            raise ValueError('retained secrets, slots or input ciphertext changed')
        read(current, n*112)  # The refreshed output is the allowed difference.
        if current.read(1):
            raise ValueError('unexpected retained-function prefix suffix')


def validate(state, directory):
    report = json.loads((directory/'report.json').read_text())
    function = profile_gate(report)
    same_retained_input(state, directory/'replayed-prefix.bin', function, report['n'])
    phase = inspect_prefix(directory/'replayed-prefix.bin')
    square = inspect_square(state, directory/'square.bin', function)
    if phase['n'] != report['n'] or phase['function'] != function or square['n'] != report['n']:
        raise ValueError('independent retained-function profile mismatch')
    for independent, mapping in [
        (phase, {k:k for k in ['fresh_input_max_error', 'max_error', 'rms_error']}),
        (square, {'post_refresh_product_max_error':'post_refresh_product_max_error',
                  'rms_error':'square_rms_error', 'target_amplitude':'target_amplitude',
                  'product_tolerance':'product_tolerance'})]:
        for key, reported in mapping.items():
            if not math.isclose(independent[key], report[reported], rel_tol=1e-11, abs_tol=1e-11):
                raise ValueError('independent retained-function metric mismatch '+reported)
    refresh_passed = max(phase[k] for k in ['fresh_input_max_error', 'max_error', 'rms_error']) < 2**-23
    cpp_bounds_passed = max(report[k] for k in ['evaluation_max_error', 'complex_uniform_bound']) < 2**-23
    cpp_bounds_passed &= report['post_refresh_product_evaluation_max_error'] < report['evaluation_product_tolerance']
    passed = refresh_passed and square['product_accuracy_passed'] and cpp_bounds_passed
    if report.get('accuracy_passed') is not passed:
        raise ValueError('retained-function accuracy gate disagrees')
    return dict(scope=__doc__, function=function, n=report['n'], independent_prefix=phase,
                independent_square=square, same_retained_secrets_slots_and_input=True,
                fixed_23bit_refresh_passed=refresh_passed, accuracy_passed=passed,
                evaluation_and_uniform_bounds_source='C++ report', fresh_keys_generated=False,
                security_validated=False, selector_bit_error_origin_established=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['state', 'directory', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    a = parser.parse_args()
    result = validate(a.state, a.directory)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Retained-function independent phases, input identity and accuracy gate agree')


if __name__ == '__main__':
    main()
