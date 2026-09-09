"""Conditional, correlation-safe error-envelope composition for a 16-bit CLPX reverse.

The theorem takes marginal primitive tail envelopes as hypotheses. The adapter
uses Parameter-Selection standard deviations as *assumed sub-Gaussian scales*;
variance estimates alone do not establish those hypotheses. No independence
between calls, reused keys, or input coefficients is assumed in composition.
"""
from __future__ import annotations
import math
from . import clpx as noise
from .params import clpx as defaults
from .keyvariation import brroundnoise

Q = 1 << 64
MASK = Q-1

def distance_to_sign_boundary(x: int) -> int:
    x &= MASK
    return min(x, abs(x-Q//2), Q-x)

def oracle_margins(basebit=4):
    """Exact noiseless HomDecomp/bit-extraction trace over all 16-bit residues.

Only subtract-PBS decisions whose outputs feed requested bits are obligations.
All shifts and decisions are integer operations, not floating-point samples.
"""
    if basebit not in (2,4):
        raise ValueError('only the implemented basebit 2 and 4 circuits are supported')
    numdigit=16//basebit+1
    carry_digits=1 if basebit==4 else 2
    digits=numdigit+carry_digits
    source_bits=basebit*digits
    submins=[Q]*(numdigit-1)
    outmins=[Q]*16
    witnesses=[None]*(numdigit-1+16)
    for value in range(1<<16):
        acc=(value<<44)+(1<<42)  # value/2^20 plus the C++ 2^-22 bias
        sub=0
        for j in range(numdigit):
            cres=((acc << (source_bits-basebit*(j+1)))+
                  (sub-(1<<(63-basebit)) if j else 0)) & MASK
            centered=(cres+(1<<(63-basebit))) & MASK
            if j<numdigit-1:
                margin=distance_to_sign_boundary(centered)
                if margin<submins[j]:
                    submins[j]=margin; witnesses[j]=value
            if j:
                for k in range(basebit):
                    index=(j-1)*basebit+k
                    phase=(centered << (basebit-k-1)) & MASK
                    if int(phase>=Q//2)!=((value>>index)&1):
                        raise AssertionError((basebit,value,index,phase))
                    margin=distance_to_sign_boundary(phase)
                    if margin<outmins[index]:
                        outmins[index]=margin; witnesses[numdigit-1+index]=value
            sub=(1 if centered<Q//2 else -1)*(1<<(62-basebit))
    return dict(basebit=basebit,numdigit=numdigit,carry_digits=carry_digits,
                source_bits=source_bits,verified_residues=1<<16,
                hom_margins=[x/Q for x in submins],output_margins=[x/Q for x in outmins],
                hom_margin_numerators=submins,output_margin_numerators=outmins,
                denominator=Q,witnesses=witnesses)

def primitive_scales(half_n=760, half_alpha_bits=17, ks20_levels=7, ks20_basebit=2):
    class Half(defaults.lvlhalfparam):
        n=half_n
        alpha=2.0**(32-half_alpha_bits)
        α=alpha
        σ=alpha**2
    class Ring2(defaults.lvl2param):
        l=4
        lₐ=4
        ℬbit=10
        ℬₐbit=10
        ℬ=1<<10
        ℬₐ=1<<10
    class B2:
        domainP=Half
        targetP=Ring2
    class B1:
        domainP=Half
        targetP=defaults.lvl1param
    class K20(defaults.lvl2hparam):
        targetP=Half
        t=ks20_levels
        basebit=ks20_basebit
    class K1h(defaults.lvl1hparam):
        targetP=Half
    return {
        'ks20':math.sqrt(noise.identity_key_switch_variance(K20))/2**32,
        'ks21':math.sqrt(noise.identity_key_switch_variance(defaults.lvl21param))/2**32,
        'ks1h':math.sqrt(noise.identity_key_switch_variance(K1h))/2**32,
        'br2':math.sqrt(float(brroundnoise(B2)))/2**32,
        'br1':math.sqrt(float(brroundnoise(B1)))/2**32,
        'pbs2':math.sqrt(noise.pbs_variance(B2))/2**64,
        'pbs1':math.sqrt(noise.pbs_variance(B1))/2**32,
    }

def compose(margins, radii):
    """Deterministic sufficient conditions; all radii in normalized torus units.

Works for arbitrary dependent errors satisfying the supplied absolute envelopes.
Fid is handled modulo 1/16: its possible lift changes telescope in radix 2.
"""
    r=radii
    fid=(r['input']+2*r['ks20']+2*r['br2']+r['pbs2'])/8+r['pbs2']
    constraints=[dict(stage="digit_range",index=0,
                      error=3*(fid+r["pbs2"]),radius=1/32)]
    for i in range(16):
        constraints.append(dict(stage='digit_round',index=i,
            error=(2 if i==0 else 3)*fid+r['ks20']+r['br2'],radius=1/32))
    acc=16*r['pbs2']
    bb=margins['basebit']; nd=margins['numdigit']; source=margins['source_bits']
    for j in range(nd):
        cres=acc*2**(source-bb*(j+1))+r['ks21']+(r['pbs1'] if j else 0)
        mid=cres+r['ks1h']
        if j<nd-1:
            constraints.append(dict(stage='hom_subtract',index=j,
                error=mid+r['br1'],radius=margins['hom_margins'][j]))
        if j:
            for k in range(bb):
                i=(j-1)*bb+k
                constraints.append(dict(stage='output_bit_pbs',index=i,
                    error=mid*2**(bb-k-1)+r['br1'],radius=margins['output_margins'][i]))
    for i in range(16):
        constraints.append(dict(stage='output_decrypt',index=i,error=r['pbs1'],radius=1/8))
    return constraints

def estimate(*,input_sigma,basebit=4,half_n=760,half_alpha_bits=17,multiplications=256,
             ks20_levels=7,ks20_basebit=2):
    if input_sigma<0 or not math.isfinite(input_sigma) or multiplications<=0:
        raise ValueError('finite nonnegative input sigma and positive run length required')
    if not 1 <= ks20_basebit <= 16 or not 1 <= ks20_levels or ks20_basebit*ks20_levels > 32:
        raise ValueError("KS20 decomposition must retain 1..32 bits, with basebit 1..16")
    margins=oracle_margins(basebit)
    scales=primitive_scales(half_n,half_alpha_bits,ks20_levels,ks20_basebit)
    scales['input']=input_sigma
    # Extra deterministic scalar quantization allowances, separate from the
    # assumed centered stochastic rounding envelopes. FFT error is NOT covered.
    biases={k:0.0 for k in scales}
    for k in ('ks20','ks21','ks1h','br2','br1'): biases[k]=2**-32
    for k in ('pbs1','pbs2'): biases[k]=2**(-32 if k=='pbs1' else -64)
    unit=compose(margins,scales); bias_constraints=compose(margins,biases)
    limits=[]
    for c,b in zip(unit,bias_constraints):
        c['bias']=b['error']
        c['cutoff_limit']=(c['radius']-c['bias'])/c['error'] if c['error'] else math.inf
        limits.append(c['cutoff_limit'])
    maximum=min(limits)
    # Counts include unused executed HomDecomp work, conservatively. Every
    # scalar event need only obey a marginal bound; no product of probabilities.
    nd=margins['numdigit']; digits=nd+margins['carry_digits']
    counts={'input':16,'ks20':48,'br2':48,'pbs2':48,
            'ks21':digits,'ks1h':digits+nd-1,
            'br1':digits-1+16,'pbs1':digits-1+16}
    events=sum(counts.values())*multiplications
    target_cutoff=math.sqrt(2*math.log(2*events*2**40))
    # Reserve 1% error-radius slack; strict inequalities are required.
    cutoff=0.99*maximum
    log2_bound=min(0.0,math.log2(2*events)-cutoff**2/(2*math.log(2)))
    for c in unit:
        c['error_at_target_cutoff']=target_cutoff*c['error']+c['bias']
        c['passes_at_target_cutoff']=c['error_at_target_cutoff']<c['radius']
    return dict(scope='One 16-bit-output block; whole run of configured length. '
                      'Input noise and primitive sub-Gaussian tail hypotheses are explicit.',
        parameters=dict(basebit=basebit,numdigit=nd,source_bits=margins['source_bits'],
                        half_n=half_n,half_alpha_bits=half_alpha_bits,
                        ks20_levels=ks20_levels,ks20_basebit=ks20_basebit),
        multiplications=multiplications,primitive_scales=scales,primitive_bias_allowances=biases,
        primitive_events_per_product=counts,total_primitive_events=events,
        maximum_uniform_cutoff_sigma=maximum,used_cutoff_sigma=cutoff,
        cutoff_required_for_target=target_cutoff,
        conditional_log2_failure_bound=log2_bound,
        conditional_target_passes=log2_bound<=-40,
        bottleneck=min(unit,key=lambda c:c['cutoff_limit']),constraints=unit,
        oracle=margins,
        implementation_certified=False,
        unproved_hypotheses=[
          'Each primitive error has two-sided tail <= 2 exp(-t^2/2) at the supplied scale; Parameter-Selection variances alone do not prove this.',
          'The input CLPX coefficients have the specified marginal tails around a valid radix-2 encoding.',
          'Numerical FFT error is bounded within the primitive envelopes; the current adapter has no certified FFT envelope.',
        ])
