#!/usr/bin/env python3
"""Independently validate packed identity checkpoint input, refresh and square.

Reads only the private checkpoint prefix, uses exact integer secret products,
and decodes canonical slots without TFHEpp. Evaluation keys and selector phases
are not inspected. N512 evidence cannot establish N16384 accuracy or security.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct

from narrow_ckks_packed_cpu_compare import EXPECTED, digest
from narrow_ckks_prefix_arithmetic import PhaseDecoder, canonical_decode


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
            raise ValueError('truncated packed checkpoint/ciphertext')
        return data

    def integers(self, count):
        return struct.unpack(self.endian+str(count)+'Q', self.read(8*count))

    def polynomial(self, n):
        limbs = self.integers(7*n)
        return [sum(limbs[7*i+j] << (64*j) for j in range(7)) for i in range(n)]

    def ciphertext(self, n):
        return self.polynomial(n), self.polynomial(n)


def inspect(checkpoint, refresh, square):
    with checkpoint.open('rb') as file:
        reader = Reader(file)
        profile = reader.integers(9)
        n = profile[1]
        if n not in (512, 16384) or profile != (
                0x4e415252504b3031, n, 85, 41, 379, 49, 34, 100, 448):
            raise ValueError('unsupported packed checkpoint profile')
        length = reader.integers(1)[0]
        if length != 4 or reader.read(length) != b'tanh':
            raise ValueError('unsupported packed identity checkpoint function tag')
        secret = [x if x < 2**447 else x-2**448 for x in reader.polynomial(n)]
        sparse = [x if x < 2**447 else x-2**448 for x in reader.polynomial(n)]
        if any(x not in (-1, 0, 1) for x in secret+sparse):
            raise ValueError('nonternary packed checkpoint secret')
        width = (n+30)//31
        if any(sum(x != 0 for x in sparse[i*width:min((i+1)*width,n)]) != 1
               for i in range(31)):
            raise ValueError('invalid packed structured secret geometry')
        raw = struct.unpack(reader.endian+str(n)+'d', reader.read(n*8))
        slots = [complex(*raw[2*i:2*i+2]) for i in range(n//2)]
        if any(not math.isfinite(x.real) or not math.isfinite(x.imag) or
               abs(x.real) > 0.125 or abs(x.imag) > 0.125 for x in slots):
            raise ValueError('invalid packed identity target domain')
        input_ct = reader.ciphertext(n)
        prefix_bytes = file.tell()
        file.seek(0)
        prefix_sha256 = hashlib.sha256(file.read(prefix_bytes)).hexdigest()
    decoder = PhaseDecoder(secret)

    def decode(ciphertext, q, scale):
        phase = decoder.phase(*ciphertext, q)
        got = canonical_decode(phase, scale, direct=n == 512)
        if any(not math.isfinite(x.real) or not math.isfinite(x.imag) for x in got):
            raise ValueError('nonfinite independent packed decode')
        return got

    def saved_ciphertext(path):
        with path.open('rb') as file:
            ciphertext = Reader(file).ciphertext(n)
            if file.read(1):
                raise ValueError('unexpected packed ciphertext suffix')
        return ciphertext

    fresh = decode(input_ct, 85, 41)
    refreshed = decode(saved_ciphertext(refresh), 100, 49)
    squared = decode(saved_ciphertext(square), 51, 49)
    errors = [abs(x-y) for x,y in zip(refreshed, slots)]
    metrics = dict(fresh_input_max_error=max(abs(x-y) for x,y in zip(fresh,slots)),
                   max_error=max(errors),
                   rms_error=math.sqrt(math.fsum(x*x for x in errors)/len(errors)),
                   post_refresh_product_max_error=max(abs(x-y*y) for x,y in zip(squared,slots)))
    tolerance = 2**-16
    if any(not math.isfinite(x) or not 0 <= x < tolerance for x in metrics.values()):
        raise ValueError('independent packed input/refresh/square accuracy gate failed')
    return dict(n=n, metrics=metrics, tolerance=tolerance,
                checkpoint_prefix_bytes_read=prefix_bytes,
                checkpoint_prefix_sha256=prefix_sha256, checkpoint=str(checkpoint),
                independent_phase_validated=True, evaluation_keys_loaded=False,
                selector_phases_validated=False, security_validated=False,
                full_ring_phase_validated=n == 16384)


def validate(checkpoint, refresh, square, report_path):
    report = json.loads(report_path.read_text())
    n = report.get('n')
    if type(n) is not int or n not in (512, 16384):
        raise ValueError('unsupported packed report ring')
    profile = dict(EXPECTED, n=n, packed_log_q=413+n.bit_length()-1,
                   selector_ntt_chunk_size=min(40,(n+30)//31))
    for name, expected in profile.items():
        value = report.get(name)
        if value != expected or type(value) is not type(expected):
            raise ValueError('incompatible packed report profile: '+name)
    if report.get('accuracy_passed') is not True:
        raise ValueError('packed report accuracy failed')
    if report.get('accuracy_failure_injected', False) is not False:
        raise ValueError('injected packed accuracy failure')
    if any(report.get(key) != 2**-16 for key in ('evaluation_tolerance','target_error_tolerance')):
        raise ValueError('changed packed report accuracy gate')
    result = inspect(checkpoint, refresh, square)
    if n != result['n']:
        raise ValueError('checkpoint/report ring mismatch')
    for name, actual in result['metrics'].items():
        value = report.get(name)
        if type(value) not in (int,float) or not math.isfinite(value) or not math.isclose(
                value, actual, abs_tol=1e-11, rel_tol=1e-11):
            raise ValueError('independent packed metric mismatch: '+name)
    value = report.get('post_refresh_product_evaluation_max_error')
    if type(value) not in (int,float) or not math.isfinite(value) or not math.isclose(
            value, result['metrics']['post_refresh_product_max_error'], abs_tol=1e-11, rel_tol=1e-11):
        raise ValueError('independent packed square evaluation metric mismatch')
    result['accuracy_passed'] = True
    result['artifacts'] = {name: dict(path=str(path),sha256=digest(path)) for name,path in
                           [('refresh',refresh),('square',square),('report',report_path)]}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('checkpoint','refresh','square','report','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    result = validate(args.checkpoint,args.refresh,args.square,args.report)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Independent packed input, refresh and square accuracy validated at N'+str(result['n']))


if __name__ == '__main__':
    main()
