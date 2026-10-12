#!/usr/bin/env python3
"""Decode basis-checkpoint outputs with independent integer phase arithmetic.

Reads diagnostic secrets without printing them. Checks the refresh, square and
normalized square; evaluation/approximation bounds remain C++ report metrics.
Does not qualify evaluation-key coefficients or security.
"""
import argparse
import json
import math
import struct
from pathlib import Path
from narrow_ckks_prefix_arithmetic import PRIME, ntt, PhaseDecoder, canonical_decode


class BasisPhaseDecoder:
    def __init__(self, secret):
        self.secret = secret
        self.n = len(secret)
        self.spectrum = None
        if self.n not in (512, 32768) or any(x not in (-1, 0, 1) for x in secret):
            raise ValueError('unsupported basis secret')
        # Fourteen-bit limbs bound the full N32768 signed convolution.
        if self.n * ((1 << 14) - 1) >= PRIME // 2:
            raise ValueError('basis phase scratch modulus is too small')

    def phase(self, a, b, q):
        if self.n == 512:
            return PhaseDecoder(self.secret).phase(a, b, q)
        if len(a) != self.n or len(b) != self.n or not 0 < q <= 640:
            raise ValueError('invalid basis ciphertext shape or level')
        modulus = 1 << q
        a = [x % modulus for x in a]
        if self.spectrum is None:
            self.spectrum = ntt(self.secret + [0] * self.n)
        product = [0] * self.n
        for shift in range(0, q, 14):
            transformed = ntt([(x >> shift) & 16383 for x in a] + [0] * self.n)
            convolution = ntt([x*y % PRIME for x, y in zip(transformed, self.spectrum)], True)
            for i in range(self.n):
                value = (convolution[i] - convolution[i+self.n]) % PRIME
                if value > PRIME // 2:
                    value -= PRIME
                product[i] += value << shift
        return [(y-x+modulus//2) % modulus-modulus//2 for x, y in zip(product, b)]


class Reader:
    def __init__(self, file):
        self.file = file
        flag = self.read(1)[0]
        if flag not in (0, 1):
            raise ValueError('invalid portable archive byte order')
        self.endian = '<' if flag else '>'

    def read(self, count):
        data = self.file.read(count)
        if len(data) != count:
            raise ValueError('truncated basis checkpoint data')
        return data

    def words(self, count):
        return struct.unpack(self.endian+str(count)+'Q', self.read(count*8))

    def polynomial(self, n, limbs):
        words = self.words(n*limbs)
        return [sum(words[i*limbs+j] << (64*j) for j in range(limbs)) for i in range(n)]


def validate(checkpoint, result, report_path, function):
    report = json.loads(report_path.read_text())
    rows = report.get('results', [report])
    matches = [r for r in rows if r.get('function') == function]
    if len(matches) != 1:
        raise ValueError('missing or duplicate basis function report')
    reported = matches[0]
    with checkpoint.open('rb') as file:
        reader = Reader(file)
        metadata = reader.words(12)
        magic, n, storage, selector, delta, bits, input_q, input_delta, bg, bbar, auxiliary, s2c = metadata
        if (magic != 0x4e41525242433031 or n not in (512, 32768) or
                (storage, selector, delta, bits, bg, bbar) not in ((576, 518, 47, 31, 18, 8), (640, 536, 49, 32, 4, 4)) or
                (input_q, input_delta, auxiliary, s2c) != (85, 44, 44, 37)):
            raise ValueError('unsupported basis checkpoint profile')
        limbs = storage // 64
        secret = [x if x < 1 << (storage-1) else x-(1 << storage)
                  for x in reader.polynomial(n, limbs)]
        reader.polynomial(n, limbs)  # Structured secret is outside this gate.
        values = struct.unpack(reader.endian+str(n)+'d', reader.read(n*8))
        slots = [complex(values[2*i], values[2*i+1]) for i in range(n//2)]
        input_ct = [reader.polynomial(n, limbs) for _ in range(2)]
        # Evaluation-key suffix is deliberately not read by this decoder.
    if any(not math.isfinite(x.real) or not math.isfinite(x.imag) or
           abs(x.real) > 4 or abs(x.imag) > 4 for x in slots):
        raise ValueError('invalid basis input domain')
    with result.open('rb') as file:
        reader = Reader(file)
        if reader.words(12) != metadata:
            raise ValueError('basis output profile mismatch')
        length = reader.words(1)[0]
        if length != len(function) or reader.read(length) != function.encode():
            raise ValueError('basis output function mismatch')
        ciphertexts = [[reader.polynomial(n, limbs) for _ in range(2)] for _ in range(3)]
        if file.read(1):
            raise ValueError('basis output has trailing data')
    q = selector - 39 - 6*delta
    scale = delta + bits
    if (reported['n'], reported['output_log_q'], reported['output_log_delta']) != (n, q, scale):
        raise ValueError('basis report profile mismatch')
    decoder = BasisPhaseDecoder(secret)
    def decode(ct, level, exponent):
        return canonical_decode(decoder.phase(*ct, level), exponent, direct=n == 512)
    fresh = decode(input_ct, input_q, input_delta)
    output = decode(ciphertexts[0], q, scale)
    square = decode(ciphertexts[1], q-scale, scale)
    normalized = decode(ciphertexts[2], input_q, input_delta)
    def target(x):
        if function == 'sigmoid': return 1/(1+math.exp(-x))
        if function == 'gelu': return x*(1+math.erf(x/math.sqrt(2)))/2
        if function == 'tanh': return math.tanh(x)
        raise ValueError('unsupported basis function')
    targets = [complex(target(x.real), target(x.imag)) for x in slots]
    tolerance = 2**-23
    bound = 2*max(map(abs, targets))*tolerance+tolerance*tolerance+tolerance
    metrics = dict(fresh_input_max_error=max(abs(a-b) for a,b in zip(fresh, slots)),
                   max_error=max(abs(a-b) for a,b in zip(output, targets)),
                   post_refresh_product_max_error=max(abs(a-b*b) for a,b in zip(square, targets)),
                   square_normalization_max_error=max(abs(a-b) for a,b in zip(normalized, square)))
    for name in ['max_error', 'post_refresh_product_max_error', 'square_normalization_max_error']:
        if not math.isclose(metrics[name], reported[name], rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError('independent basis metric mismatch '+name)
    if not math.isclose(bound, reported['product_tolerance'], rel_tol=1e-11, abs_tol=1e-11):
        raise ValueError('basis product bound mismatch')
    passed = (max(metrics[k] for k in ['fresh_input_max_error','max_error','square_normalization_max_error']) < tolerance
              and metrics['post_refresh_product_max_error'] < bound)
    return dict(scope=__doc__, n=n, function=function, metrics=metrics,
                product_tolerance=bound, independent_phase_accuracy_passed=passed,
                security_validated=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['checkpoint', 'result', 'report', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--function', choices=['sigmoid','gelu','tanh'], required=True)
    args = parser.parse_args()
    data = validate(args.checkpoint, args.result, args.report, args.function)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    raise SystemExit(0 if data['independent_phase_accuracy_passed'] else 1)
