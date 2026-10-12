#!/usr/bin/env python3
"""Exercise the independent packed phase validator with synthetic archive controls.

These files have no evaluation keys and are not encrypted refresh evaluations.
Retained production ciphertext validation is separate evidence.
"""
import argparse
import json
from pathlib import Path
import struct
import tempfile

from narrow_ckks_packed_checkpoint_validate import inspect, validate
from narrow_ckks_packed_cpu_compare import EXPECTED


def fixture(root, n, endian):
    root.mkdir()
    def write_poly(file, values):
        for value in values:
            value %= 1 << 448
            file.write(struct.pack(endian+'7Q',*((value >> (64*j)) & (2**64-1) for j in range(7))))
    secret = [1]+[0]*(n-1)
    sparse = [0]*n
    for i in range(31):
        sparse[i*((n+30)//31)] = (-1)**i
    def ciphertext(file, q, scale, real, imag):
        a = [((i+1)*0xdeadbeef123456789abcdef) % (1 << q) for i in range(n)]
        b = list(a)
        b[0] = (b[0]+int(real*(1 << scale))) % (1 << q)
        b[n//2] = (b[n//2]+int(imag*(1 << scale))) % (1 << q)
        write_poly(file,a);write_poly(file,b)
    checkpoint = root/'checkpoint.bin'
    with checkpoint.open('wb') as file:
        file.write(bytes([int(endian=='<')]))
        file.write(struct.pack(endian+'9Q',0x4e415252504b3031,n,85,41,379,49,34,100,448))
        file.write(struct.pack(endian+'Q',4)+b'tanh')
        write_poly(file,secret);write_poly(file,sparse)
        file.write(struct.pack(endian+str(n)+'d',*([1/32,1/64]*(n//2))))
        ciphertext(file,85,41,1/32,1/64)
    paths=[checkpoint]
    for label,q,real,imag in [('refresh',100,1/32,1/64),('square',51,3/4096,1/1024)]:
        p=root/(label+'.ciphertext');paths.append(p)
        with p.open('wb') as file:
            file.write(bytes([int(endian=='<')]))
            ciphertext(file,q,49,real,imag)
    return paths


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    cases=[];rejections=[]
    with tempfile.TemporaryDirectory(prefix='packed-phase-controls-') as tmp:
        root=Path(tmp)
        for n,endian in [(512,'<'),(512,'>'),(16384,'<')]:
            paths=fixture(root/('n'+str(n)+('le' if endian=='<' else 'be')),n,endian)
            result=inspect(*paths)
            assert max(result['metrics'].values()) < 1e-11
            cases.append(dict(n=n,endian=endian,maximum_error=max(result['metrics'].values())))
            if n==512 and endian=='<':baseline=paths
        checkpoint,refresh,square=baseline
        original=checkpoint.read_bytes()
        prefix=inspect(*baseline)
        report=dict(EXPECTED,accuracy_passed=True,evaluation_tolerance=2**-16,target_error_tolerance=2**-16,
                    **prefix['metrics'],post_refresh_product_evaluation_max_error=prefix['metrics']['post_refresh_product_max_error'])
        report_path=root/'report.json';report_path.write_text(json.dumps(report))
        assert validate(*baseline,report_path)['accuracy_passed']
        def reject(label,call):
            try:call()
            except ValueError:rejections.append(label)
            else:raise AssertionError('accepted '+label)
        for label,data in [('bad-byte-order',b'\x02'+original[1:]),
                           ('truncated-header',original[:10]),
                           ('wrong-magic',original[:1]+bytes(8)+original[9:]),
                           ('wrong-function',original[:81]+b'gelu'+original[85:]),
                           ('nonternary-secret',original[:85]+struct.pack('<Q',2)+original[93:]),
                           ('invalid-structured-secret',original[:85+512*56]+bytes(512*56)+original[85+1024*56:]),
                           ('truncated-input',original[:-16])]:
            checkpoint.write_bytes(data);reject(label,lambda:inspect(*baseline))
        checkpoint.write_bytes(original)
        for label,path,data in [('truncated-refresh',refresh,refresh.read_bytes()[:-8]),
                                ('square-suffix',square,square.read_bytes()+b'x'),
                                ('refresh-used-as-square',square,refresh.read_bytes()),
                                ('square-used-as-refresh',refresh,square.read_bytes())]:
            old=path.read_bytes();path.write_bytes(data);reject(label,lambda:inspect(*baseline));path.write_bytes(old)
        for label,field,value in [('reported-metric','max_error',0.001),('reported-ring','n',16384),
                                  ('reported-accuracy','accuracy_passed',1),('changed-gate','evaluation_tolerance',1),
                                  ('nonfinite-square','post_refresh_product_evaluation_max_error',float('nan'))]:
            bad=dict(report);bad[field]=value;report_path.write_text(json.dumps(bad));reject(label,lambda:validate(*baseline,report_path))
    evidence=dict(scope=__doc__,all_passed=True,synthetic_archive_cases=cases,rejected_controls=rejections,
                  rejection_count=len(rejections),production_full_ring_accuracy_validated=False,security_validated=False)
    args.output.write_text(json.dumps(evidence,indent=2)+'\n');print('Packed independent reader: synthetic archives and '+str(len(rejections))+' rejection controls passed')


if __name__=='__main__':main()
