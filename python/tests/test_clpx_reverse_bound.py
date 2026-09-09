"""Deterministic oracle/envelope tests; these do not empirically prove tails."""
import math
import random
import unittest
from noiseestimation.clpx_reverse_bound import oracle_margins, compose, estimate


def fold(x):
    x=x%1
    return x if x<0.5 else 0.5-x

def rounded_digit(x):
    x=x%1
    sign=1 if x<0.5 else -1
    interval=math.floor(16*(x%0.5))
    return sign*min(interval,8-interval)

def scalar_reverse(value,basebit,radii,rng):
    """Plaintext scalar trace of the C++ circuit with bounded error injections."""
    def err(k): return rng.choice((-1,1))*radii[k]
    nd=16//basebit+1
    source=basebit*(nd+(1 if basebit==4 else 2))
    acc=0.; previous=0.
    for j in range(16):
        u=-value/2**(j+2)+err('input')
        first=fold(u+err('ks20')+err('br2'))+err('pbs2')
        v=fold(first+err('ks20')+err('br2'))/8+err('pbs2')
        d=rounded_digit(previous-2*v+err('ks20')+err('br2')+1/32)
        acc+=d/2**(20-j)+err('pbs2')
        previous=v
    acc+=2**-22
    sub=0.; result=0
    for j in range(nd):
        cres=acc*2**(source-basebit*(j+1))+err('ks21')
        if j: cres+=sub-2**(-basebit-1)
        inp=cres+err('ks1h')+2**(-basebit-1)
        if j<nd-1:
            sub=(1 if (inp+err('br1'))%1<0.5 else -1)*2**(-basebit-2)+err('pbs1')
        if j:
            # C++ uses a separate IKS for output extraction, shared over this digit's bits.
            out_in=cres+2**(-basebit-1)+err('ks1h')
            for k in range(basebit):
                phase=(out_in*2**(basebit-k-1)+err('br1'))%1
                output=(1/8 if phase>=0.5 else -1/8)+err('pbs1')
                actual=output>0
                result|=int(actual)<<((j-1)*basebit+k)
    return result

class ReverseBoundTests(unittest.TestCase):
    def test_full_domain_semantics_and_exact_margins(self):
        for bb,last in ((4,4097/262144),(2,16385/262144)):
            m=oracle_margins(bb)
            self.assertEqual(m['verified_residues'],65536)
            self.assertEqual(m['output_margins'][-1],last)
            self.assertTrue(all(x>0 for x in m['hom_margins']))

    def test_fold_composition_and_perturbation_at_wraps(self):
        for i in range(4096):
            x=i/4096
            actual=fold(fold(x))/8
            dist=lambda a,b:abs((a-b+1/32)%(1/16)-1/32)
            self.assertLess(dist(actual,x/8),1e-14)
            for error in (-1/4096,1/4096):
                self.assertLessEqual(dist(fold(fold(x)+error)/8,x/8),abs(error)/8+1e-14)

    def test_bounded_error_injections(self):
        for bb,n,a,levels,digitbits in ((4,760,17,7,2),(2,832,19,7,2),(2,760,17,2,8),(2,760,17,3,5)):
            report=estimate(input_sigma=2**(109.20114828158572/2-64),basebit=bb,
                            half_n=n,half_alpha_bits=a,ks20_levels=levels,ks20_basebit=digitbits)
            t=.9*report['maximum_uniform_cutoff_sigma']
            radii={k:t*v+report['primitive_bias_allowances'][k] for k,v in report['primitive_scales'].items()}
            self.assertTrue(all(c['error']<c['radius'] for c in compose(report['oracle'],radii)))
            rng=random.Random(20260906)
            for value in [0,1,2,3,32767,32768,65025,65535]+[rng.randrange(65536) for _ in range(512)]:
                self.assertEqual(scalar_reverse(value,bb,radii,rng),value,(bb,value))

    def test_fail_closed_and_run_union(self):
        sigma=2**(109.20114828158572/2-64)
        old=estimate(input_sigma=sigma)
        new=estimate(input_sigma=sigma,basebit=2,half_n=832,half_alpha_bits=19)
        one=estimate(input_sigma=sigma,basebit=2,half_n=832,half_alpha_bits=19,multiplications=1)
        fixed=estimate(input_sigma=sigma,basebit=2,ks20_levels=2,ks20_basebit=8)
        small=estimate(input_sigma=sigma,basebit=2,ks20_levels=3,ks20_basebit=5)
        self.assertLess(fixed['conditional_log2_failure_bound'],new['conditional_log2_failure_bound'])
        self.assertTrue(small['conditional_target_passes'])
        self.assertEqual(fixed['parameters']['half_n'],760)
        self.assertEqual(fixed['parameters']['half_alpha_bits'],17)
        self.assertFalse(old['conditional_target_passes'])
        self.assertTrue(new['conditional_target_passes'])
        self.assertFalse(new['implementation_certified'])
        self.assertAlmostEqual(new['conditional_log2_failure_bound']-one['conditional_log2_failure_bound'],8)
        with self.assertRaises(ValueError):oracle_margins(3)
        with self.assertRaises(ValueError):estimate(input_sigma=float('nan'))

if __name__=='__main__':unittest.main()
