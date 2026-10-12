#!/usr/bin/env python3
"""Exact compatible sparse dictionaries; no lattice attack or ring reuse claim.

Forced blocks have no zero projection. Assign each forced block to one side,
then split the remaining nonzero entries between two partial dictionaries.
Filter pairs sharing a block. Every full projection has exactly binomial(w-s,k)
compatible decompositions, including signs and positions. Rank/unrank constructs
individual entries without materializing the astronomically large dictionaries.
"""
import argparse
from collections import Counter
import datetime
import hashlib
import itertools
import json
from math import comb,log2,prod
from pathlib import Path

from narrow_ckks_block_rotation_probability import coefficients


class Dictionary:
    def __init__(self, guessed, forced, weight):
        self.guessed=tuple(guessed)
        self.forced=frozenset(forced)
        self.weight=weight
        if (type(weight) is not int or not 0<=weight<=len(guessed) or
            any(not 0<=j<len(guessed) or guessed[j]<1 for j in self.forced) or
            any(type(g) is not int or g<0 for g in guessed)):
            raise ValueError('invalid partial dictionary')
        self.suffix=[[0]*(weight+1) for _ in range(len(guessed)+1)]
        self.suffix[-1][0]=1
        for j in range(len(guessed)-1,-1,-1):
            for k in range(weight+1):
                self.suffix[j][k]=(0 if j in self.forced else self.suffix[j+1][k])
                if k:self.suffix[j][k]+=2*guessed[j]*self.suffix[j+1][k-1]
        self.size=self.suffix[0][weight]

    def unrank(self,rank):
        if type(rank) is not int or not 0<=rank<self.size:
            raise ValueError('dictionary rank out of range')
        left=self.weight
        result=[]
        for j,g in enumerate(self.guessed):
            zero=0 if j in self.forced else self.suffix[j+1][left]
            if rank<zero:continue
            rank-=zero
            if left==0:raise AssertionError('nonzero branch at zero weight')
            step=self.suffix[j+1][left-1]
            choice,rank=divmod(rank,step)
            assert choice<2*g
            result.append((j,choice//2,-1 if choice%2==0 else 1))
            left-=1
        assert left==0 and rank==0
        return tuple(result)

    def rank(self,entry):
        selected={j:(p,s) for j,p,s in entry}
        if len(selected)!=len(entry) or len(entry)!=self.weight or not self.forced<=selected.keys():
            raise ValueError('invalid dictionary entry')
        if any(not 0<=j<len(self.guessed) or not 0<=p<self.guessed[j] or s not in (-1,1)
               for j,(p,s) in selected.items()):
            raise ValueError('invalid signed block position')
        rank=0;left=self.weight
        for j in range(len(self.guessed)):
            if j not in selected:continue
            if j not in self.forced:rank+=self.suffix[j+1][left]
            p,s=selected[j]
            rank+=(2*p+int(s==1))*self.suffix[j+1][left-1]
            left-=1
        assert left==0 and 0<=rank<self.size
        return rank


def split(widths,guessed,weight):
    if len(widths)!=len(guessed) or any(not 0<=g<=w for w,g in zip(widths,guessed)):
        raise ValueError('invalid full block geometry')
    forced=[j for j,(w,g) in enumerate(zip(widths,guessed)) if g==w]
    _,support=coefficients(widths,guessed,weight)
    if support==0:raise ValueError('empty projection stratum')
    optional=[g for j,g in enumerate(guessed) if j not in forced]
    polynomial=[1]+[0]*weight
    for g in optional:
        for k in range(weight,0,-1):polynomial[k]+=2*g*polynomial[k-1]
    best=None
    for cut in range(len(forced)+1):
        for left_weight in range(cut,weight-(len(forced)-cut)+1):
            k=left_weight-cut
            left=prod(2*guessed[j] for j in forced[:cut])*polynomial[k]
            right=prod(2*guessed[j] for j in forced[cut:])*polynomial[weight-len(forced)-k]
            if not left or not right:continue
            candidate=(left+right,cut,left_weight,left,right)
            if best is None or candidate<best:best=candidate
    work,cut,left_weight,left,right=best
    left_guessed=[0 if j in forced[cut:] else g for j,g in enumerate(guessed)]
    right_guessed=[0 if j in forced[:cut] else g for j,g in enumerate(guessed)]
    a=Dictionary(left_guessed,forced[:cut],left_weight)
    b=Dictionary(right_guessed,forced[cut:],weight-left_weight)
    assert (a.size,b.size)==(left,right)
    multiplicity=comb(weight-len(forced),left_weight-cut)
    return a,b,dict(forced_blocks=forced,left_forced_blocks=forced[:cut],
        right_forced_blocks=forced[cut:],left_weight=left_weight,right_weight=weight-left_weight,
        left_entries=str(left),right_entries=str(right),enumerated_entries=str(work),
        full_projection_count=str(support),compatible_decompositions_per_projection=multiplicity,
        compatible_pair_count=str(multiplicity*support),all_pair_count=str(left*right),
        enumeration_log2=log2(work),square_root_support_log2=log2(support)/2,
        enumeration_over_square_root_support_bits=log2(work)-log2(support)/2,
        optimization_scope='All side weights and contiguous splits of the public forced-block order only')


def self_test():
    cases=0
    widths=[2,2,1]
    for guessed in itertools.product(*(range(w+1) for w in widths)):
        for weight in range(len(widths)+1):
            _,count=coefficients(widths,guessed,weight)
            if not count:continue
            a,b,record=split(widths,guessed,weight)
            options=[]
            for w,g in zip(widths,guessed):
                choices=[None] if g<w else []
                choices += [(p,s) for p in range(g) for s in (-1,1)]
                options.append(choices)
            full={tuple((j,*value) for j,value in enumerate(values) if value is not None)
                  for values in itertools.product(*options)
                  if sum(v is not None for v in values)==weight}
            assert len(full)==count
            partial_a=[a.unrank(i) for i in range(a.size)]
            partial_b=[b.unrank(i) for i in range(b.size)]
            for dictionary,entries in ((a,partial_a),(b,partial_b)):
                assert len(set(entries))==dictionary.size
                assert all(dictionary.rank(entry)==i for i,entry in enumerate(entries))
            joined=Counter()
            for left,right in itertools.product(partial_a,partial_b):
                if {j for j,_,_ in left}&{j for j,_,_ in right}:continue
                joined[tuple(sorted(left+right))]+=1
            assert joined.keys()==full
            assert set(joined.values())=={record['compatible_decompositions_per_projection']}
            assert sum(joined.values())==int(record['compatible_pair_count'])
            cases+=1
    return cases


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile-report',type=Path,required=True)
    parser.add_argument('--weight',type=int)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();cases=self_test()
    profile=json.loads(args.profile_report.read_text())
    anchor=profile['attacks']['bdd_mitm_hybrid']['best']
    guessed=anchor['guessed_per_block'];widths=profile['block_widths']
    weight=anchor['guessed_weight'] if args.weight is None else args.weight
    if (sum(widths)!=profile['n'] or len(widths)!=profile['h'] or
        sum(guessed)!=anchor['grid_zeta'] or any(type(w) is not int or w<1 for w in widths)):
        raise ValueError('profile block/mask mismatch')
    a,b,record=split(widths,guessed,weight)
    samples=[]
    for name,dictionary in (('left',a),('right',b)):
        for rank in sorted({0,dictionary.size//2,dictionary.size-1}):
            entry=dictionary.unrank(rank);assert dictionary.rank(entry)==rank
            samples.append(dict(side=name,rank=str(rank),entry=entry,
                                coefficient_positions=[sum(widths[:j])+p for j,p,_ in entry]))
    report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=__doc__,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        profile_report=str(args.profile_report),profile_sha256=hashlib.sha256(args.profile_report.read_bytes()).hexdigest(),
        n=profile['n'],h=profile['h'],zeta=anchor['grid_zeta'],rotation=0,weight=weight,
        widths=widths,guessed_per_block=guessed,split=record,samples=samples,
        exhaustive_toy_cases=cases,complete=True,security_validated=False,
        full_attack_constructed=False,ring_dictionary_reuse_established=False,
        limitations=['Constructs exact signed projection dictionaries at rotation zero only',
            'Pair compatibility filters overlapping blocks; this is not a Babai collision search',
            'Enumerated entries exclude vector arithmetic, hashing, collisions, decoding and lattice reduction',
            'No square-root speedup over all rotation-labelled candidates is established',
            'The entry-count overhead is not a security margin or a lower bound on other attack algorithms',
            'The dictionaries are represented implicitly; full sizes are not materialized'])
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('PASS:',cases,'exhaustive signed toy strata; actual full-size rank/unrank; enumeration overhead bits',record['enumeration_over_square_root_support_bits'])


if __name__=='__main__':main()
