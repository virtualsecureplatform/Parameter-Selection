#!/usr/bin/env python3
"""Inspect auxiliary8 N512/N16384 prefixes using independent phase arithmetic.

Reads secrets only from a private diagnostic file. Outputs metrics, never key
material. Does not inspect evaluation-key suffixes or retained-square keys.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
from narrow_ckks_prefix_arithmetic import PhaseDecoder, canonical_decode, canonical_component_decode

def inspect(path):
    with path.open('rb') as file:
        def read(size):
            data = file.read(size)
            if len(data) != size:
                raise ValueError('truncated diagnostic prefix')
            return data
        endian_flag = read(1)[0]
        if endian_flag not in (0, 1):
            raise ValueError('invalid archive endian flag')
        endian = '<' if endian_flag else '>'
        metadata = struct.unpack(endian+'10Q', read(80))
        n = metadata[1]
        if n not in (512, 16384) or metadata != (0x4155583850523031, n, 85, 44, 385, 48, 37, 108, 8, 448):
            raise ValueError('wrong auxiliary probe profile')
        length = struct.unpack(endian+'Q', read(8))[0]
        if length > 7:
            raise ValueError('invalid function length')
        function = read(length).decode('ascii')
        if function not in ('sigmoid', 'gelu', 'tanh') or read(1) != b'\x01':
            raise ValueError('invalid function or unfinished refresh')
        def polynomial():
            values = struct.unpack(endian+str(n*7)+'Q', read(n*56))
            return [sum(values[7*i+j] << (64*j) for j in range(7)) for i in range(n)]
        def secret():
            values = polynomial()
            signed = [x if x < 2**447 else x-2**448 for x in values]
            if any(x not in (-1, 0, 1) for x in signed):
                raise ValueError('secret is not ternary')
            return signed
        dense, sparse = secret(), secret()
        width = (n+30)//31
        if any(sum(x != 0 for x in sparse[i*width:min((i+1)*width, n)]) != 1 for i in range(31)):
            raise ValueError('sparse key violates one-per-block geometry')
        raw_slots = struct.unpack(endian+str(n)+'d', read(n*8))
        slots = [complex(*raw_slots[2*i:2*i+2]) for i in range(n//2)]
        ciphertexts = [(polynomial(), polynomial()) for _ in range(2)]
        prefix_bytes = file.tell()
    if any(not math.isfinite(x.real) or not math.isfinite(x.imag) for x in slots):
        raise ValueError('nonfinite saved slots')
    if any(abs(x.real) > 4 or abs(x.imag) > 4 for x in slots):
        raise ValueError('saved slots exceed the probe message domain')
    decoder = PhaseDecoder(dense)
    def decrypt(ciphertext, log_q, scale):
        a, b = ciphertext
        return canonical_decode(decoder.phase(a, b, log_q), scale, direct=n == 512)
    fresh = decrypt(ciphertexts[0], 85, 44)
    refreshed_phase = decoder.phase(*ciphertexts[1],108)
    refreshed = canonical_decode(refreshed_phase,48,direct=n == 512)
    isolated = canonical_component_decode(refreshed_phase,48,direct=n == 512)
    def target(x):
        if function == 'sigmoid':
            return 1/(1+math.exp(-x))
        if function == 'gelu':
            return x*(1+math.erf(x/math.sqrt(2)))/2
        return math.tanh(x)
    expected = [complex(target(x.real), target(x.imag)) for x in slots]
    errors = [abs(a-b) for a, b in zip(refreshed, expected)]
    real_errors = [abs(a.real-b.real) for a, b in zip(refreshed, expected)]
    imaginary_errors = [abs(a.imag-b.imag) for a, b in zip(refreshed, expected)]
    return dict(isolated_real_component_max_error=max(abs(a.real-b.real) for a,b in zip(isolated,expected)),
                isolated_imaginary_component_max_error=max(abs(a.imag-b.imag) for a,b in zip(isolated,expected)),
                real_component_max_error=max(real_errors),
                imaginary_component_max_error=max(imaginary_errors),
                real_component_rms_error=math.sqrt(math.fsum(x*x for x in real_errors)/len(real_errors)),
                imaginary_component_rms_error=math.sqrt(math.fsum(x*x for x in imaginary_errors)/len(imaginary_errors)),
                function=function, n=n, prefix_bytes=prefix_bytes,
                phase_arithmetic='integer scalar convolution' if n == 512 else 'exact 15-bit limbs with one bounded 31-bit NTT prime',
                fresh_input_max_error=max(abs(a-b) for a, b in zip(fresh, slots)),
                max_error=max(errors),
                rms_error=math.sqrt(math.fsum(x*x for x in errors)/len(errors)),
                sparse_key_geometry_passed=True, evaluation_key_suffix_validated=False,
                retained_square_validated=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = inspect(args.state)
    report = json.loads(args.report.read_text())
    if report['function'] != result['function']:
        raise ValueError('function mismatch')
    for name in ('fresh_input_max_error', 'max_error', 'rms_error'):
        reported = report.get(name)
        if (type(reported) not in (int, float) or not math.isfinite(reported) or reported < 0
                or not math.isfinite(result[name])
                or not math.isclose(result[name], reported, rel_tol=1e-11, abs_tol=1e-11)):
            raise ValueError('independent prefix metric mismatch: '+name)
    result.update(scope=__doc__, state_bytes=args.state.stat().st_size,
                  report_sha256=hashlib.sha256(args.report.read_bytes()).hexdigest(),
                  comparison_relative_tolerance=1e-11, comparison_absolute_tolerance=1e-11,
                  fixed_23bit_refresh_passed=all(result[k] < 2**-23 for k in
                      ('fresh_input_max_error', 'max_error', 'rms_error')))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('independent diagnostic prefix metrics agree')


if __name__ == '__main__':
    main()
