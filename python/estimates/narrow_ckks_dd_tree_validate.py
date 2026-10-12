#!/usr/bin/env python3
"""Independently decode public DD product-tree outputs with a retained key.

Reads only the original private ternary secret prefix. Public outputs contain
known plaintext targets and ciphertexts, never secrets or evaluation keys.
This checks tree precision, not complete refresh performance or security.
"""
import argparse
import json
import math
from pathlib import Path
import struct

from narrow_ckks_prefix_arithmetic import PhaseDecoder, canonical_decode


def validate(state, ciphertexts, report):
    def reader(file):
        def read(count):
            value = file.read(count)
            if len(value) != count:
                raise ValueError('truncated diagnostic archive')
            return value
        marker = read(1)
        if marker not in (b'\x00', b'\x01'):
            raise ValueError('invalid portable archive byte order')
        return read, '<' if marker == b'\x01' else '>'

    def polynomial(read, endian, n):
        words = struct.unpack(endian+str(n*7)+'Q', read(n*56))
        return [sum(words[i*7+j] << (64*j) for j in range(7)) for i in range(n)]

    with state.open('rb') as file:
        read, endian = reader(file)
        meta = struct.unpack(endian+'10Q', read(80))
        n = meta[1]
        if n != 16384 or meta != (0x4155583850523031, n, 85, 44, 385, 48, 37, 108, 8, 448):
            raise ValueError('incompatible retained state profile')
        length = struct.unpack(endian+'Q', read(8))[0]
        if length != 7 or read(length) != b'sigmoid' or read(1) != b'\x01':
            raise ValueError('invalid retained state tag')
        secret = polynomial(read, endian, n)
        secret = [x if x < 2**447 else x-2**448 for x in secret]
        if any(x not in (-1, 0, 1) for x in secret):
            raise ValueError('nonternary retained secret')
    decoder = PhaseDecoder(secret)
    cpp = json.loads(report.read_text())
    if not isinstance(cpp, list) or len(cpp) != 3:
        raise ValueError('expected three C++ tree reports')
    samples = []
    with ciphertexts.open('rb') as file:
        read, endian = reader(file)
        if struct.unpack(endian+'5Q', read(40)) != (0x5452454538303031, n, 108, 48, 3):
            raise ValueError('incompatible public tree archive')
        for kind in range(3):
            if struct.unpack(endian+'Q', read(8))[0] != kind:
                raise ValueError('unexpected public tree case order')
            values = struct.unpack(endian+str(n)+'d', read(n*8))
            target = [complex(values[2*i], values[2*i+1]) for i in range(n//2)]
            if any(not math.isfinite(z.real) or not math.isfinite(z.imag) for z in target):
                raise ValueError('nonfinite known tree target')
            a, b = polynomial(read, endian, n), polynomial(read, endian, n)
            decoded = canonical_decode(decoder.phase(a, b, 108), 48)
            errors = [abs(x-y) for x, y in zip(decoded, target)]
            maximum = max(errors)
            rms = math.sqrt(math.fsum(x*x for x in errors)/len(errors))
            row = cpp[kind]
            expected = dict(tree_kind=kind, projected_primary8_dd=True, fresh_keys_generated=False,
                            encrypted_leaves=31, transparent_beta_leaves=1, depth=5,
                            output_log_q=108, output_log_delta=48, slots_checked=n//2,
                            accuracy_passed=True)
            if any(type(row.get(k)) is not type(v) or row[k] != v for k, v in expected.items()):
                raise ValueError('C++ tree report profile mismatch')
            for name, value in [('maximum_error', maximum), ('rms_error', rms)]:
                if (type(row.get(name)) not in (int, float) or not math.isfinite(row[name])
                        or row[name] < 0 or not math.isfinite(value)
                        or not math.isclose(row[name], value, rel_tol=1e-11, abs_tol=1e-11)):
                    raise ValueError('independent tree metric mismatch: '+name)
            samples.append(dict(tree_kind=kind, maximum_error=maximum, rms_error=rms,
                                fixed_23bit_passed=maximum < 2**-23))
        if file.read(1):
            raise ValueError('unexpected public tree archive suffix')
    return dict(scope=__doc__, n=n, samples=samples, slots_checked=3*n//2,
                independent_integer_phases_validated=True, cpp_metrics_match=True,
                all_checks_passed=all(row['fixed_23bit_passed'] for row in samples),
                fresh_keys_generated=False, security_validated=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('state', 'ciphertexts', 'report', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.state, args.ciphertexts, args.report)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, allow_nan=False))
    return 0 if result['all_checks_passed'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
