#!/usr/bin/env python3
"""Prepare a private one-bit selector correction for a controlled replay.

Leaves the original selectors intact. This is a diagnostic experiment, not
key generation or a production recovery mechanism.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import struct


def prepare(original, correction, destination):
    c = json.loads(correction.read_text())
    if (type(c.get('tail')) is not bool or type(c.get('group')) is not int
            or type(c.get('record')) is not int or type(c.get('coefficient')) is not int
            or c.get('component') != 0 or c.get('log_q') != 385):
        raise ValueError('unsupported correction metadata')
    delta = int(c['delta'])
    bit = abs(delta).bit_length()-1
    if (not 0 <= bit < 385 or abs(delta) != 1 << bit
            or not 0 <= c['record'] < 529 or not 0 <= c['coefficient'] < 16384
            or not 0 <= c['group'] < (5 if c['tail'] else 31)):
        raise ValueError('invalid one-bit correction')
    names = [f'head-{i}.bin' for i in range(31)]+[f'tail-{i}.bin' for i in range(5)]
    name = ('tail-' if c['tail'] else 'head-')+str(c['group'])+'.bin'
    source = original/name
    initial_stat = source.stat()
    with source.open('rb') as f:
        header_bytes = f.read(73)
        if header_bytes[0] not in (0,1):
            raise ValueError('invalid selector endian flag')
        endian = '<' if header_bytes[0] else '>'
        header = struct.unpack(endian+'9Q',header_bytes[1:])
        if header != (0x4e4152524f573031,16384,385,48,37,529,448,int(c['tail']),c['group']):
            raise ValueError('wrong selector group header')
        if initial_stat.st_size != 73+529*2*16384*56:
            raise ValueError('wrong selector file length')
        offset = 73+c['record']*2*16384*56
        f.seek(offset);record=f.read(2*16384*56)
    if hashlib.sha256(record).hexdigest() != c['record_sha256']:
        raise ValueError('original ciphertext hash mismatch')
    coefficient_offset = c['coefficient']*56
    words = struct.unpack(endian+'7Q',record[coefficient_offset:coefficient_offset+56])
    old = sum(x << (64*j) for j,x in enumerate(words))
    new = (old+delta) % 2**385
    if old ^ new != abs(delta):
        raise ValueError('requested correction is not exactly one bit')
    # Each limb follows the archive endian convention; limbs stay low to high.
    word, bit_in_word = divmod(bit,64)
    byte_in_word = bit_in_word//8 if endian == '<' else 7-bit_in_word//8
    byte_offset = offset+coefficient_offset+word*8+byte_in_word
    if any(not (original/n).is_file() for n in names):
        raise ValueError('incomplete original selector directory')
    destination.mkdir(mode=0o700)  # Refuse existing paths.
    for n in names:
        if n != name:
            os.link(original/n,destination/n)
    digest = hashlib.sha256()
    target = destination/name
    with source.open('rb') as f,target.open('xb') as out:
        target.chmod(0o600)
        while block := f.read(1024*1024):
            digest.update(block);out.write(block)
    assert source.stat().st_ino != target.stat().st_ino
    with target.open('r+b') as f:
        f.seek(byte_offset);old_byte=f.read(1)
        f.seek(byte_offset);f.write(bytes([old_byte[0] ^ (1 << (bit_in_word%8))]))
    with source.open('rb') as f:
        unchanged = hashlib.file_digest(f,'sha256').hexdigest() == digest.hexdigest()
    final_stat = source.stat()
    unchanged = unchanged and (initial_stat.st_size,initial_stat.st_mtime_ns) == (final_stat.st_size,final_stat.st_mtime_ns)
    if not unchanged:
        raise ValueError('original selector changed during diagnostic preparation')
    with target.open('rb') as f:
        corrected_hash = hashlib.file_digest(f,'sha256').hexdigest()
    return dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                scope=__doc__,directory=str(destination),selector_files=36,
                modified_group=name,modified_ciphertexts=1,modified_bits=1,
                original_group_sha256=digest.hexdigest(),corrected_group_sha256=corrected_hash,
                original_selectors_unchanged=True,other_groups_hardlinked=True,
                keys_generated=False,phase_validation_pending=True,
                functional_replay_validated=False,failure_cause_established=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,required=True)
    p.add_argument('--correction',type=Path,required=True)
    p.add_argument('--destination',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    result=prepare(a.original,a.correction,a.destination)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print('private one-bit correction prepared; original selectors unchanged')


if __name__ == '__main__':
    main()
