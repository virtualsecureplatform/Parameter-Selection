#!/usr/bin/env python3
"""Diagnose large phase errors in one retained auxiliary8 selector.

Independent integer convolution; reports aggregate metrics, never secrets or
phase coefficients. A matched impulse describes the error shape, not its cause.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct

from narrow_ckks_prefix_arithmetic import PhaseDecoder


def mask_impulse(high, dense):
    """Return a mask correction whose negacyclic secret product equals high."""
    divisor = math.gcd(*high)
    if not divisor:
        return None
    normalized = [x//divisor for x in high]
    if any(x not in (-1,0,1) for x in normalized):
        return None
    target = bytes(x+1 for x in normalized)
    sequence = bytes(x+1 for x in [-v for v in dense]+dense)
    for pattern, sign in ((target,1),(bytes(2-x for x in target),-1)):
        position = sequence.find(pattern)
        if 0 <= position <= len(dense):
            if position == 0:
                return 0, -sign*divisor
            return len(dense)-position, sign*divisor
    return None


def inspect(state, group, record, correction=None):
    n = 16384
    with state.open('rb') as f:
        flag = f.read(1)
        if flag not in (b'\x00', b'\x01'):
            raise ValueError('invalid state endian flag')
        endian = '<' if flag == b'\x01' else '>'
        metadata = struct.unpack(endian+'10Q', f.read(80))
        if metadata != (0x4155583850523031,n,85,44,385,48,37,108,8,448):
            raise ValueError('wrong state profile')
        length = struct.unpack(endian+'Q', f.read(8))[0]
        if length != 7 or f.read(7) != b'sigmoid' or f.read(1) != b'\x01':
            raise ValueError('wrong state function or unfinished')
        words = struct.unpack(endian+str(n*7)+'Q', f.read(n*56))
        dense = [sum(words[7*i+j] << (64*j) for j in range(7)) for i in range(n)]
        dense = [x if x < 2**447 else x-2**448 for x in dense]
    if any(x not in (-1,0,1) for x in dense):
        raise ValueError('nonternary dense secret')
    with group.open('rb') as f:
        flag = f.read(1)
        if flag not in (b'\x00',b'\x01'):
            raise ValueError('invalid selector endian flag')
        endian = '<' if flag == b'\x01' else '>'
        header = struct.unpack(endian+'9Q', f.read(72))
        if (header[:7] != (0x4e4152524f573031,n,385,48,37,529,448)
                or header[7] not in (0,1) or header[8] >= (5 if header[7] else 31)
                or not 0 <= record < 529 or group.stat().st_size != 73+529*2*n*56):
            raise ValueError('wrong selector metadata, length or record')
        f.seek(73+record*2*n*56)
        raw = f.read(2*n*56)
    words = struct.unpack(endian+str(2*n*7)+'Q', raw)
    values = [sum(words[7*i+j] << (64*j) for j in range(7)) for i in range(2*n)]
    phase = PhaseDecoder(dense).phase(values[:n],values[n:],385)
    # Every intended signed indicator encoding coefficient has magnitude <=
    # 2^48. Round away that bounded plaintext/noise to isolate errors >=2^64.
    high = [((x+2**63)//2**64)*2**64 for x in phase]
    residual = max(abs(x-y) for x,y in zip(phase,high))
    divisor = math.gcd(*high)
    impulse = mask_impulse(high, dense)
    result = dict(n=n, group=int(header[8]), tail=bool(header[7]), record=record,
                record_sha256=hashlib.sha256(raw).hexdigest(),
                phase_coefficients_checked=n,
                maximum_phase_bits=max(abs(x).bit_length() for x in phase),
                large_error_coefficients=sum(x != 0 for x in high),
                low_residual_maximum=residual,
                high_error_gcd_bits=divisor.bit_length(),
                high_error_gcd_is_power_of_two=bool(divisor and divisor&(divisor-1)==0),
                high_error_matches_single_mask_impulse=impulse is not None,
                intended_plaintext_and_noise_bound_compatible=residual < 2**49,
                phase_arithmetic='independent exact 15-bit limb NTT',
                failure_cause_established=False, security_validated=False)
    if impulse is not None:
        index, delta = impulse
        fixed = list(values[:n])
        fixed[index] = (fixed[index]+delta) % 2**385
        fixed_phase = PhaseDecoder(dense).phase(fixed,values[n:],385)
        result['in_memory_corrected_maximum_phase'] = max(abs(x) for x in fixed_phase)
        result['in_memory_corrected_phase_equals_low_residual'] = fixed_phase == [x-y for x,y in zip(phase,high)]
        result['mask_correction_changes_one_bit'] = (fixed[index]^values[index]).bit_count() == 1
        if correction is not None:
            # Coordinates belong to one public ciphertext, not a secret key.
            # Keep this experimental change private and leave the source intact.
            payload = dict(group=int(header[8]),tail=bool(header[7]),record=record,
                           record_sha256=result['record_sha256'],component=0,
                           coefficient=index,delta=str(delta),log_q=385)
            with correction.open('x') as output:
                correction.chmod(0o600)
                json.dump(payload,output,indent=2);output.write('\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--state',type=Path,required=True)
    p.add_argument('--group',type=Path,required=True)
    p.add_argument('--record',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--private-correction',type=Path)
    a = p.parse_args()
    result = inspect(a.state,a.group,a.record,a.private_correction)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('independent retained selector error shape inspected')


if __name__ == '__main__':
    main()
