#!/usr/bin/env python3
"""Inspect retained basis checkpoint geometry and optional payload integrity.

Reads public metadata and streams hashes without decoding or printing secrets.
Checks portable headers and fixed selector lengths. It does not validate the
checkpoint's evaluation-key suffix or ciphertext phases; successful encrypted
replay and independent accuracy gates remain required.
"""
import argparse
import hashlib
import json
import struct
from pathlib import Path


def header(path, words):
    with path.open('rb') as file:
        flag = file.read(1)
        if flag not in (b'\x00', b'\x01'):
            raise ValueError('invalid portable archive byte order: '+str(path))
        data = file.read(words*8)
        if len(data) != words*8:
            raise ValueError('truncated checkpoint metadata: '+str(path))
        return struct.unpack(('<' if flag == b'\x01' else '>')+str(words)+'Q', data)


def digest(path):
    with path.open('rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def inspect(directory, hash_payloads=False, reference=None):
    checkpoint = directory/'checkpoint.bin'
    metadata = header(checkpoint, 12)
    magic, n, storage, selector, delta, coefficient, q, scale, primary, secondary, auxiliary, plain = metadata
    if (magic != 0x4e41525242433031 or n not in (512,32768) or
            (storage,selector,delta,coefficient,primary,secondary) not in
            ((576,518,47,31,18,8),(640,536,49,32,4,4)) or
            (q,scale,auxiliary,plain) != (85,44,44,37)):
        raise ValueError('unsupported basis checkpoint profile')
    width = (n+30)//31
    tail_per_pack = (n//2)//(2*(width-1))
    tail_packs = (31+tail_per_pack-1)//tail_per_pack
    selector_bytes = width*2*n*(storage//8)
    minimum_prefix = 97+n*(4*(storage//8)+8)
    if checkpoint.stat().st_size <= minimum_prefix:
        raise ValueError('basis checkpoint has no evaluation-key suffix')
    files = [checkpoint]
    for tail, groups in [(False,31),(True,tail_packs)]:
        for group in range(groups):
            path = directory/'selectors'/(("tail-" if tail else "head-")+str(group)+'.bin')
            expected = (0x4e4152524f573031,n,selector,delta,39,width,storage,int(tail),group)
            if header(path,9) != expected:
                raise ValueError('basis selector header mismatch: '+str(path))
            if path.stat().st_size != 73+selector_bytes:
                raise ValueError('basis selector fixed length mismatch: '+str(path))
            files.append(path)
    records = {}
    for path in files:
        record = dict(bytes=path.stat().st_size)
        if hash_payloads or reference is not None:
            record['sha256'] = digest(path)
        records[str(path.relative_to(directory))] = record
    result = dict(scope=__doc__, n=n, storage_bits=storage, profile=list(metadata),
                  selector_groups_checked=len(files)-1, geometry_checked=True,
                  payloads_hashed=hash_payloads or reference is not None,
                  files=records, fresh_keys_generated=False, files_modified=False,
                  encrypted_replay_validated=False, evaluation_key_suffix_validated=False,
                  security_validated=False)
    if reference is not None:
        expected = json.loads(reference.read_text())
        if (expected.get('profile') != list(metadata) or expected.get('payloads_hashed') is not True
                or expected.get('files') != records):
            raise ValueError('basis checkpoint differs from retained integrity manifest')
        result['retained_integrity_matches'] = True
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--hash-payloads',action='store_true')
    parser.add_argument('--reference',type=Path)
    args = parser.parse_args()
    # Preserve a previous manifest; use a distinct output when verifying it.
    if args.output.exists():
        parser.error('inspection output must be new')
    data = inspect(args.directory,args.hash_payloads,args.reference)
    with args.output.open('x') as file:
        json.dump(data,file,indent=2,allow_nan=False)
        file.write('\n')
    print('Basis checkpoint headers and selector lengths checked; no keys generated')
