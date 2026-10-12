#!/usr/bin/env python3
"""Construct signed projection dictionaries for each actual public rotation.

This establishes per-rotation compatibility and exact coordinates/signs only.
It does not construct a Babai collision search or reuse dictionaries isometrically.
"""
import argparse
from collections import Counter
import datetime
import hashlib
import itertools
import json
from math import fsum,log2,sqrt
from pathlib import Path

from narrow_ckks_block_mitm_dictionaries import Dictionary,split
from narrow_ckks_block_rotation_verify import polynomial


def public_mask(widths,guessed):
    if len(widths)!=len(guessed) or any(not 0<=g<=w or w<1 for w,g in zip(widths,guessed)):
        raise ValueError('invalid public mask')
    start=0;mask=0
    for w,g in zip(widths,guessed):
        mask|=((1<<g)-1)<<start
        start+=w
    return mask


def rotation_windows(widths,mask,rotation):
    n=sum(widths)
    if type(rotation) is not int or not 0<=rotation<n:
        raise ValueError('rotation out of range')
    selected=(mask>>rotation)|((mask&((1<<rotation)-1))<<(n-rotation))
    windows=[];start=0
    for w in widths:
        windows.append((selected>>start)&((1<<w)-1))
        start+=w
    return windows


def coordinates(widths,windows,rotation,entry):
    """Map selected-position ordinals to X^r*s coordinates and signed values."""
    n=sum(widths);result=[]
    for j,ordinal,sign in entry:
        if type(ordinal) is not int or ordinal<0 or sign not in (-1,1) or not 0<=j<len(widths):
            raise ValueError('invalid signed selected-position ordinal')
        window=windows[j]
        # Select the ordinal-th set bit independently of the suffix DP.
        for _ in range(ordinal):window&=window-1
        if not window:raise ValueError('selected-position ordinal out of range')
        position=sum(widths[:j])+(window&-window).bit_length()-1
        result.append(((position+rotation)%n,-sign if position+rotation>=n else sign))
    return tuple(sorted(result))


def self_test():
    widths=[2,2,1];n=sum(widths);starts=[sum(widths[:j]) for j in range(len(widths))]
    cases=0
    for guessed in itertools.product(*(range(w+1) for w in widths)):
        mask=public_mask(widths,guessed)
        for r in range(n):
            windows=rotation_windows(widths,mask,r);counts=[v.bit_count() for v in windows]
            actual=[set() for _ in range(len(widths)+1)]
            # Independent full signed secret enumeration, then X^r negacyclic
            # substitution and public-mask projection, without rotated bitsets.
            for positions in itertools.product(*(range(w) for w in widths)):
                for signs in itertools.product((-1,1),repeat=len(widths)):
                    projection=[]
                    for j,(p,sign) in enumerate(zip(positions,signs)):
                        original=starts[j]+p;rotated=(original+r)%n
                        if (mask>>rotated)&1:
                            projection.append((rotated,-sign if original+r>=n else sign))
                    actual[len(projection)].add(tuple(sorted(projection)))
            for weight,full in enumerate(actual):
                if not full:continue
                a,b,record=split(widths,counts,weight)
                partial_a=[a.unrank(i) for i in range(a.size)]
                partial_b=[b.unrank(i) for i in range(b.size)]
                for entry in partial_a+partial_b:
                    expected=[]
                    for j,ordinal,sign in entry:
                        selected=[p for p in range(widths[j]) if (mask>>((starts[j]+p+r)%n))&1]
                        original=starts[j]+selected[ordinal]
                        expected.append(((original+r)%n,-sign if original+r>=n else sign))
                    assert coordinates(widths,windows,r,entry)==tuple(sorted(expected))
                joined=Counter()
                for left,right in itertools.product(partial_a,partial_b):
                    if {j for j,_,_ in left}&{j for j,_,_ in right}:continue
                    joined[coordinates(widths,windows,r,left+right)]+=1
                assert joined.keys()==full
                assert set(joined.values())=={record['compatible_decompositions_per_projection']}
                assert len(full)==int(record['full_projection_count'])
                cases+=1
    return cases


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile-report',type=Path,required=True)
    parser.add_argument('--probability-report',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--summary-only',action='store_true',help='retain exact aggregates and coordinate samples, omitting per-rotation rows')
    args=parser.parse_args();toy=self_test()
    profile=json.loads(args.profile_report.read_text())
    posterior=json.loads(args.probability_report.read_text())
    anchor=profile['attacks']['bdd_mitm_hybrid']['best']
    widths=profile['block_widths'];guessed=anchor['guessed_per_block'];weight=posterior['guessed_weight']
    if (not posterior['complete'] or posterior['widths']!=widths or
        posterior['guessed_per_block']!=guessed or posterior['n']!=profile['n'] or
        posterior['h']!=profile['h'] or sum(widths)!=profile['n'] or
        len(widths)!=profile['h'] or sum(guessed)!=anchor['grid_zeta'] or
        posterior['profile_report_sha256']!=hashlib.sha256(args.profile_report.read_bytes()).hexdigest()):
        raise ValueError('posterior/public geometry mismatch')
    mask=public_mask(widths,guessed);cache={};points=[];samples=[]
    support_sum=0;entry_sum=0;roots=[];support_digest=hashlib.sha256()
    for r in range(profile['n']):
        windows=rotation_windows(widths,mask,r);counts=tuple(v.bit_count() for v in windows)
        if counts not in cache:
            groups=tuple(Counter(zip(widths,counts)).items())
            expected=polynomial(groups,weight,True)[weight]
            if expected:
                a,b,record=split(widths,counts,weight)
                assert int(record['full_projection_count'])==expected
                for dictionary in (a,b):
                    count=polynomial(tuple(Counter(zip(widths,dictionary.guessed)).items()),dictionary.weight,True)[dictionary.weight]
                    assert count==dictionary.size
                cache[counts]=record
            else:cache[counts]=None
        record=cache[counts]
        support_digest.update(f'{r}:{0 if record is None else record["full_projection_count"]}\n'.encode())
        if record is None:
            points.append(dict(rotation=r,supported=False));continue
        support=int(record['full_projection_count']);entries=int(record['enumerated_entries'])
        support_sum+=support;entry_sum+=entries;roots.append(sqrt(support))
        points.append(dict(rotation=r,supported=True,full_projection_count=str(support),
            enumerated_entries=str(entries),left_weight=record['left_weight'],
            left_forced_blocks=record['left_forced_blocks'],
            compatible_decompositions_per_projection=record['compatible_decompositions_per_projection']))
        if r in (0,1,profile['n']//2,profile['n']-1):
            a,b,_=split(widths,counts,weight)
            for name,dictionary in (('left',a),('right',b)):
                rank=dictionary.size//2;entry=dictionary.unrank(rank)
                assert dictionary.rank(entry)==rank
                mapped=coordinates(widths,windows,r,entry)
                assert all((mask>>p)&1 for p,_ in mapped)
                samples.append(dict(rotation=r,side=name,rank=str(rank),entry=entry,rotated_signed_coefficients=mapped))
    if support_sum!=int(posterior['rotation_labelled_candidate_count']):
        raise ValueError('all-rotation support count differs from saved exact report')
    report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=__doc__,
        sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),Path(__file__).with_name('narrow_ckks_block_mitm_dictionaries.py'),
             args.profile_report,args.probability_report)},
        n=profile['n'],h=profile['h'],weight=weight,zeta=anchor['grid_zeta'],toy_cases=toy,
        unique_count_geometries=len(cache),rotation_count=len(points),
        points=None if args.summary_only else points,points_omitted=args.summary_only,samples=samples,
        labelled_support_count=str(support_sum),rotation_support_sequence_sha256=support_digest.hexdigest(),
        separate_dictionary_entry_sum=str(entry_sum),
        separate_dictionary_entry_sum_log2=log2(entry_sum),
        square_root_labelled_support_log2=log2(support_sum)/2,
        separate_square_root_support_sum_log2=log2(fsum(roots)),
        complete=True,security_validated=False,full_attack_constructed=False,
        ring_dictionary_reuse_established=False,
        limitations=['Exact per-rotation projection decomposition only; no Babai collision search',
            'Separate enumeration is a concrete construction, not a lower bound on shared work',
            'Entry counts exclude lattice reduction, vector arithmetic, hashing and collision filtering',
            'Reported enumeration overhead is not a security margin; existing cost screen is unchanged'])
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('PASS:',toy,'signed negacyclic toy strata;',len(points),'actual rotations; labelled support exact; separate enumeration bits',report['separate_dictionary_entry_sum_log2'],flush=True)


if __name__=='__main__':main()
