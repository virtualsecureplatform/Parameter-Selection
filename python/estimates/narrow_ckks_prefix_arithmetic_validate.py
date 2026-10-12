#!/usr/bin/env python3
"""Validate full-ring prefix arithmetic with algebraic fixtures, not FBS runs.

All full-ring phase coefficients have independently constructed references.
These fixtures do not run encrypted functional bootstrapping or inspect keys.
"""
import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path
import random
import struct
import tempfile

from narrow_ckks_prefix_arithmetic import PRIME, PhaseDecoder, canonical_decode, ntt
from narrow_ckks_auxiliary_probe_prefix import inspect


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert PRIME > 2 and all(PRIME % d for d in range(3, math.isqrt(PRIME)+1, 2))
    randomizer = random.Random(20262119)
    toy_cases = 0
    for n in (8, 16, 32):
        for bits in (1, 15, 16, 85, 108):
            modulus = 1 << bits
            secret = [randomizer.choice((-1, 0, 1)) for _ in range(n)]
            a = [randomizer.randrange(modulus) + (1 << (bits+40)) for _ in range(n)]
            b = [randomizer.randrange(modulus) for _ in range(n)]
            reference = []
            for i in range(n):
                product = sum(a[(i-j) % n]*sign*(1 if i >= j else -1)
                              for j, sign in enumerate(secret))
                reference.append((b[i]-product+modulus//2) % modulus-modulus//2)
            decoder = PhaseDecoder(secret)
            assert decoder.phase(a, b, bits, True) == reference
            assert decoder.phase(a, b, bits) == reference
            toy_cases += 1
    # These controls reject unsupported shapes instead of applying an
    # incorrect root or silently weakening the signed-convolution bound.
    controls = []
    for name, operation in (
            ('non_power_of_two_transform', lambda: ntt([1, 2, 3])),
            ('nonternary_secret', lambda: PhaseDecoder([1, 2])),
            ('oversized_ring', lambda: PhaseDecoder([0]*32768)),
            ('wrong_ciphertext_shape', lambda: PhaseDecoder([1, -1]).phase([1], [0, 0], 85))):
        try:
            operation()
        except ValueError:
            controls.append(name)
        else:
            raise AssertionError('accepted '+name)
    n = 16384
    dense = [1]*n
    decoder = PhaseDecoder(dense)
    phases, ciphertexts, analytical_slots, decoded_slots = [], [], [], []
    for bits, scale, terms in (
            (85, 44, {0: 1/4, 3: 1/8, n-1: -1/32}),
            (108, 48, {0: 1/2, 2: 1/16, n-1: 1/64})):
        modulus = 1 << bits
        a_value = modulus-1
        a = [a_value+(1 << (bits+35))]*n
        phase = [0]*n
        for index, coefficient in terms.items():
            phase[index] = int(coefficient*(1 << scale))
        # (1+...+X^(N-1))^2 mod (X^N+1) has coefficient
        # 2*i-N+2. This independent closed form reaches the signed-limb
        # bound and checks carries without reusing either convolution code.
        b = [((2*i-n+2)*a_value+phase[i]) % modulus+(1 << (bits+34)) for i in range(n)]
        measured = decoder.phase(a, b, bits)
        assert measured == phase
        slots = canonical_decode(measured, scale)
        expected = []
        for index in range(n//2):
            z = cmath.exp(1j*math.pi*pow(5, index, 2*n)/n)
            value = sum(coefficient*(z**power if power != n-1 else -z.conjugate())
                        for power, coefficient in terms.items())
            expected.append(value)
        error = max(abs(a-b) for a, b in zip(slots, expected))
        assert error < 1e-11
        phases.append(phase)
        ciphertexts.append((a, b))
        analytical_slots.append(expected)
        decoded_slots.append(dict(log_q=bits, log_delta=scale, maximum_analytical_error=error))
    # A synthetic prefix exercises the full-size archive reader as well as
    # the math. It deliberately has no evaluation-key suffix and its output
    # is not a sigmoid FBS result.
    private = Path(tempfile.mkdtemp(prefix='aux8-prefix-algebraic-',
        dir='/home/work-bleaker/kmatsuoka/ckks-diagnostics'))
    state = private/'synthetic-prefix.bin'
    sparse = [0]*n
    for i in range(31):
        sparse[i*((n+30)//31)] = -1 if i % 2 else 1
    def polynomial(file, values):
        mask = (1 << 64)-1
        payload = bytearray()
        for value in values:
            value %= 1 << 448
            payload.extend(struct.pack('<7Q', *[(value >> (64*j)) & mask for j in range(7)]))
        file.write(payload)
    with state.open('wb') as file:
        file.write(b'\x01')
        file.write(struct.pack('<10Q', 0x4155583850523031, n, 85, 44, 385, 48, 37, 108, 8, 448))
        file.write(struct.pack('<Q', 7)+b'sigmoid'+b'\x01')
        polynomial(file, dense)
        polynomial(file, sparse)
        for slot in analytical_slots[0]:
            file.write(struct.pack('<2d', slot.real, slot.imag))
        for a, b in ciphertexts:
            polynomial(file, a)
            polynomial(file, b)
    result = inspect(state)
    expected_target = [complex(1/(1+math.exp(-z.real)), 1/(1+math.exp(-z.imag)))
                       for z in analytical_slots[0]]
    errors = [abs(a-b) for a, b in zip(analytical_slots[1], expected_target)]
    assert result['n'] == n and result['fresh_input_max_error'] < 1e-11
    assert abs(result['max_error']-max(errors)) < 1e-11
    assert abs(result['rms_error']-math.sqrt(math.fsum(x*x for x in errors)/len(errors))) < 1e-11
    assert result['max_error'] > 2**-23
    base = Path(__file__).parent
    report = dict(scope=__doc__, toy_cases=toy_cases, rejected_controls=controls,
                  full_phase_coefficients_checked=2*n, full_decoded_slots_checked=n,
                  full_analytical_slot_checks=decoded_slots,
                  full_archive_prefix_result=result,
                  synthetic_prefix_path=str(state), synthetic_prefix_only=True,
                  full_encrypted_fbs_accuracy_validated=False, security_validated=False,
                  source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                      (Path(__file__), base/'narrow_ckks_prefix_arithmetic.py',
                       base/'narrow_ckks_auxiliary_probe_prefix.py')})
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print('15 toy cases and 32768 full algebraic phase coefficients passed; synthetic full prefix checked')


if __name__ == '__main__':
    main()
