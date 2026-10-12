#!/usr/bin/env python3
"""Independently verify retained auxiliary8 prepared coefficient ciphertexts.

Uses exact integer phase arithmetic and coefficient decoding, without TFHEpp
FFT routines. Reads private secrets; emits only error metrics and file hashes.
Does not validate evaluation keys, nonlinear refresh, or security.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct

from narrow_ckks_prefix_arithmetic import PhaseDecoder


def inspect(state, prepared):
    with state.open('rb') as file:
        def read(size):
            value = file.read(size)
            if len(value) != size:
                raise ValueError('truncated state prefix')
            return value
        flag = read(1)[0]
        if flag not in (0, 1):
            raise ValueError('invalid state endian flag')
        endian = '<' if flag else '>'
        metadata = struct.unpack(endian+'10Q', read(80))
        n = metadata[1]
        if n not in (512, 16384) or metadata != (0x4155583850523031, n, 85, 44, 385, 48, 37, 108, 8, 448):
            raise ValueError('wrong auxiliary8 state profile')
        length = struct.unpack(endian+'Q', read(8))[0]
        if length != 7 or read(length) != b'sigmoid' or read(1) != b'\x01':
            raise ValueError('requires finished sigmoid state')
        def polynomial():
            words = struct.unpack(endian+str(n*7)+'Q', read(n*56))
            return [sum(words[7*i+j] << (64*j) for j in range(7)) for i in range(n)]
        def secret():
            values = [v if v < 2**447 else v-2**448 for v in polynomial()]
            if any(v not in (-1, 0, 1) for v in values):
                raise ValueError('nonternary saved secret')
            return values
        secret()  # Validate the dense secret, but decrypt with the sparse one.
        sparse = secret()
        width = (n+30)//31
        if any(sum(v != 0 for v in sparse[x*width:min((x+1)*width, n)]) != 1 for x in range(31)):
            raise ValueError('invalid sparse block geometry')
        slots = struct.unpack(endian+str(n)+'d', read(n*8))
        prefix_bytes = file.tell()
    if any(not math.isfinite(v) or abs(v) > 4 for v in slots):
        raise ValueError('invalid saved message domain')
    raw = prepared.read_bytes()
    if len(raw) != 1+2*n*56 or raw[0] not in (0, 1):
        raise ValueError('invalid prepared archive length or endian flag')
    words = struct.unpack(('<' if raw[0] else '>')+str(2*n*7)+'Q', raw[1:])
    values = [sum(words[7*i+j] << (64*j) for j in range(7)) for i in range(2*n)]
    phase = PhaseDecoder(sparse).phase(values[:n], values[n:], 48)
    error = max(abs(math.ldexp(float(phase[i]), -44)-slots[2*i]) for i in range(n//2))
    error = max(error, max(abs(math.ldexp(float(phase[i+n//2]), -44)-slots[2*i+1]) for i in range(n//2)))
    return dict(n=n, function='sigmoid', prepared_max_error=error,
                phase_coefficients_checked=n, state_prefix_bytes_read=prefix_bytes,
                prepared_bytes=len(raw), prepared_sha256=hashlib.sha256(raw).hexdigest(),
                sparse_key_geometry_passed=True, log_q=48, log_delta=44,
                fixed_23bit_prepared_phase_passed=error < 2**-23,
                phase_arithmetic='exact integer negacyclic convolution',
                decoding='coefficient real/imaginary halves; no canonical FFT',
                evaluation_keys_validated=False, full_refresh_validated=False, security_validated=False)


def validate(state, prepared, report):
    result = inspect(state, prepared)
    reported = json.loads(report.read_text())
    if (reported.get('n') != result['n'] or reported.get('function') != 'sigmoid'
            or reported.get('auxiliary_bits') != 8 or reported.get('mode') != '--prepare'
            or reported.get('fresh_keys_generated') is not False):
        raise ValueError('wrong preparation report profile')
    value = reported.get('prepared_max_error')
    if (type(value) not in (int, float) or not math.isfinite(value) or value < 0
            or not math.isclose(value, result['prepared_max_error'], rel_tol=0, abs_tol=1e-13)):
        raise ValueError('independent prepared phase metric mismatch')
    result.update(report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
                  comparison_absolute_tolerance=1e-13)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.state, args.prepared, args.report)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('independent prepared coefficient phase agrees')


if __name__ == '__main__':
    main()
