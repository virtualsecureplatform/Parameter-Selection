#!/usr/bin/env python3
"""Snapshot or verify every aux8 selector file; does not validate encryption.

A snapshot establishes file identity only from the moment it is recorded.
It cannot identify the origin of an earlier error, certify selector phases,
or provide authentication against an actor able to rewrite the manifest.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import stat
import struct

MAGIC = 0x4e4152524f573031


def geometry(n):
    if type(n) is not int or n not in (512, 16384):
        raise ValueError('unsupported auxiliary8 ring')
    width = (n + 30) // 31
    tails_per_pack = (n // 2) // (2 * (width - 1))
    return width, (31 + tails_per_pack - 1) // tails_per_pack


def identity(s):
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)


def inspect_group(path, n, tail, group):
    width, _ = geometry(n)
    flags = os.O_RDONLY | os.O_NOFOLLOW
    with os.fdopen(os.open(path, flags), 'rb') as f:
        before = os.fstat(f.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError('selector is not a regular file')
        size = 73 + width * 2 * n * 7 * 8
        if before.st_size != size:
            raise ValueError('selector file length mismatch')
        prefix = f.read(73)
        if len(prefix) != 73 or prefix[0] not in (0, 1):
            raise ValueError('invalid selector archive header')
        metadata = struct.unpack(('<' if prefix[0] else '>') + '9Q', prefix[1:])
        if metadata != (MAGIC, n, 385, 48, 37, width, 448, tail, group):
            raise ValueError('selector group profile mismatch')
        digest = hashlib.sha256(prefix)
        for chunk in iter(lambda: f.read(1 << 20), b''):
            digest.update(chunk)
        after = os.fstat(f.fileno())
        current = os.stat(path, follow_symlinks=False)
        if identity(before) != identity(after) or identity(after) != identity(current):
            raise ValueError('selector file changed while hashing')
        return dict(bytes=size, sha256=digest.hexdigest(),
                    device=before.st_dev, inode=before.st_ino,
                    mtime_ns=before.st_mtime_ns, ctime_ns=before.st_ctime_ns)


def snapshot(directory, n):
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError('selector directory must be a real directory')
    _, tails = geometry(n)
    groups = {}
    for tail, count in [(False, 31), (True, tails)]:
        for group in range(count):
            name = ('tail-' if tail else 'head-') + str(group) + '.bin'
            groups[name] = inspect_group(directory / name, n, tail, group)
    # Recheck earlier groups' identities before publishing a whole-directory
    # snapshot. This detects ordinary changes during the bounded-memory scan.
    for name, entry in groups.items():
        expected = tuple(entry[k] for k in ['device', 'inode', 'bytes', 'mtime_ns', 'ctime_ns'])
        if identity(os.stat(directory / name, follow_symlinks=False)) != expected:
            raise ValueError('selector directory changed during snapshot')
    return dict(format='aux8-selector-integrity-v1',
                utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                directory=str(directory.resolve()), n=n, groups=groups,
                payload_hash_algorithm='sha256', selector_phases_validated=False,
                error_origin_established=False, security_validated=False)


def verify(directory, manifest):
    if manifest.get('format') != 'aux8-selector-integrity-v1':
        raise ValueError('unsupported selector integrity manifest')
    current = snapshot(directory, manifest.get('n'))
    if set(current['groups']) != set(manifest.get('groups', {})):
        raise ValueError('selector manifest group set mismatch')
    for name, entry in current['groups'].items():
        expected = manifest['groups'][name]
        if entry['bytes'] != expected.get('bytes') or entry['sha256'] != expected.get('sha256'):
            raise ValueError('selector payload changed: ' + name)
    return dict(utc=current['utc'], n=current['n'], groups_checked=len(current['groups']),
                payloads_match_snapshot=True, selector_phases_validated=False,
                error_origin_established=False, security_validated=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['snapshot', 'verify'])
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--n', type=int, choices=[512, 16384])
    a = parser.parse_args()
    if a.mode == 'snapshot':
        if a.n is None:
            parser.error('snapshot requires --n')
        if a.manifest.resolve().parent == a.directory.resolve():
            parser.error('store the manifest outside the selector directory')
        result = snapshot(a.directory, a.n)
        # Never replace a baseline. A later observation needs a new filename.
        with a.manifest.open('x') as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write('\n')
        print('Recorded', len(result['groups']), 'selector group hashes')
    else:
        result = verify(a.directory, json.loads(a.manifest.read_text()))
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
