#!/usr/bin/env python3
"""Independently decode a retained aux8 post-refresh square ciphertext.

Checks the sigmoid target and the existing 23-bit propagated product bound.
Does not validate every retained evaluation-key coefficient or security.
"""
import argparse
import json
import math
from pathlib import Path
import struct

from narrow_ckks_prefix_arithmetic import PhaseDecoder, canonical_decode


def inspect(state, square, function='sigmoid'):
    if function not in ('sigmoid', 'gelu', 'tanh'):
        raise ValueError('unsupported retained square target function')
    def reader(file):
        def read(count):
            data = file.read(count)
            if len(data) != count:
                raise ValueError('truncated square diagnostic input')
            return data
        flag = read(1)[0]
        if flag not in (0, 1):
            raise ValueError('invalid portable archive byte order')
        return read, '<' if flag else '>'

    def polynomial(read, endian, n):
        limbs = struct.unpack(endian+str(n*7)+'Q', read(n*56))
        return [sum(limbs[7*i+j] << (64*j) for j in range(7)) for i in range(n)]

    with state.open('rb') as file:
        read, endian = reader(file)
        profile = struct.unpack(endian+'10Q', read(80))
        n = profile[1]
        if n not in (512, 16384) or profile != (0x4155583850523031, n, 85, 44, 385, 48, 37, 108, 8, 448):
            raise ValueError('unsupported retained square profile')
        length = struct.unpack(endian+'Q', read(8))[0]
        if length != 7 or read(length) != b'sigmoid' or read(1) != b'\x01':
            raise ValueError('square requires finished sigmoid state')
        secret = [x if x < 2**447 else x-2**448 for x in polynomial(read, endian, n)]
        if any(x not in (-1, 0, 1) for x in secret):
            raise ValueError('retained square secret is not ternary')
        read(n*56)  # Sparse secret is unused for this dense-key ciphertext.
        raw = struct.unpack(endian+str(n)+'d', read(n*8))
        slots = [complex(*raw[2*i:2*i+2]) for i in range(n//2)]
    if any(not math.isfinite(x.real) or not math.isfinite(x.imag) or
           abs(x.real) > 4 or abs(x.imag) > 4 for x in slots):
        raise ValueError('invalid square target domain')
    with square.open('rb') as file:
        read, endian = reader(file)
        a, b = polynomial(read, endian, n), polynomial(read, endian, n)
        if file.read(1):
            raise ValueError('unexpected square ciphertext suffix')
    got = canonical_decode(PhaseDecoder(secret).phase(a, b, 60), 48, direct=n == 512)
    def target_value(x):
        if function == 'sigmoid':
            return 1/(1+math.exp(-x))
        if function == 'gelu':
            return x*(1+math.erf(x/math.sqrt(2)))/2
        return math.tanh(x)
    target = [complex(target_value(x.real), target_value(x.imag)) for x in slots]
    errors = [abs(x-y*y) for x, y in zip(got, target)]
    amplitude = max(map(abs, target))
    tolerance = 2**-23
    bound = 2*amplitude*tolerance+tolerance*tolerance+tolerance
    return dict(n=n, function=function, log_q=60, log_delta=48,
                post_refresh_product_max_error=max(errors),
                rms_error=math.sqrt(math.fsum(x*x for x in errors)/len(errors)),
                target_amplitude=amplitude, product_tolerance=bound,
                product_accuracy_passed=max(errors) < bound,
                independent_phase_validated=True, security_validated=False,
                retained_key_coefficients_validated=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['state', 'square', 'report', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--function', choices=['sigmoid', 'gelu', 'tanh'], default='sigmoid')
    args = parser.parse_args()
    result = inspect(args.state, args.square, args.function)
    report = json.loads(args.report.read_text())
    for key in ['post_refresh_product_max_error', 'rms_error', 'target_amplitude', 'product_tolerance']:
        if not math.isfinite(report[key]) or not math.isclose(report[key], result[key], abs_tol=1e-11, rel_tol=1e-11):
            raise ValueError('independent square metric mismatch '+key)
    for key in ['n', 'function', 'log_q', 'log_delta', 'product_accuracy_passed']:
        if report[key] != result[key]:
            raise ValueError('independent square profile or gate mismatch '+key)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('independent retained square phase and metrics agree')


if __name__ == '__main__':
    main()
