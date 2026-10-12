#!/usr/bin/env python3
"""Check the 24 weighted roots of an auxiliary8 degree12 complex sigmoid replay.

The C++ observer measures roots against the actual prepared sparse phase.
This checks trace consistency, not independent decryption or final accuracy.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def inspect(path, partial=False):
    raw=path.read_bytes()
    last=raw.rfind(b'\n')+1
    if last != len(raw) and not partial:
        raise ValueError('unfinished root trace line')
    complete=raw[:last]
    lines=complete.splitlines()
    if len(lines)>24 or (not partial and len(lines)!=24):
        raise ValueError('wrong root observation count')
    points=[]
    for i,line in enumerate(lines):
        row=json.loads(line)
        if set(row)!= {'frequency','weight_abs','max_error','max_abs'}:
            raise ValueError('wrong root trace schema')
        if any(type(row[k]) not in (int,float) or not math.isfinite(row[k]) for k in row):
            raise ValueError('nonfinite or nonnumeric root trace metric')
        if row['frequency'] != i%12+1:
            raise ValueError('unexpected root order or frequency')
        if row['weight_abs']<=0 or row['max_error']<0 or row['max_abs']<0:
            raise ValueError('invalid root metric sign')
        # Every expected root has modulus weight_abs. Reverse triangle
        # inequality bounds the maximum observed modulus by max_error.
        if abs(row['max_abs']-row['weight_abs']) > row['max_error']+1e-13:
            raise ValueError('root amplitude contradicts reported error')
        if i>=12 and row['weight_abs']!=points[i-12]['weight_abs']:
            raise ValueError('real/imaginary root weights differ')
        points.append(row)
    failures=[i for i,row in enumerate(points) if row['max_error']>=2**-23]
    return dict(scope=__doc__,root_count=len(points),expected_root_count=24,
                trace_complete=len(points)==24 and last==len(raw),
                ignored_unfinished_bytes=len(raw)-last,
                complete_prefix_sha256=hashlib.sha256(complete).hexdigest(),
                maximum_observed_root_error=max((p['max_error'] for p in points),default=None),
                maximum_observed_root_amplitude=max((p['max_abs'] for p in points),default=None),
                first_root_exceeding_23bit_gate=(dict(observation=failures[0],**points[failures[0]]) if failures else None),
                observed_roots_below_23bit=bool(points) and not failures,
                independent_phase_validated=False,final_refresh_accuracy_validated=False,
                failure_cause_established=False,security_validated=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--trace',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--partial',action='store_true')
    a=p.parse_args();r=inspect(a.trace,a.partial)
    a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('checked',r['root_count'],'root observations; complete',r['trace_complete'])


if __name__ == '__main__':
    main()
