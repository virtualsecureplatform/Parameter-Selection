"""Check exact real/imaginary phase splitting against analytical slot polynomials.

Synthetic algebra checks only; this does not validate encrypted refresh.
"""
import math,json,argparse
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output",type=Path,required=True)
args=parser.parse_args()
from narrow_ckks_prefix_arithmetic import canonical_component_decode,canonical_decode
cases=[]
for n in [16,512,16384]:
 for large_imaginary in [True,False]:
  phase=[0]*n;scale=48
  real_amplitude=.5 if large_imaginary else 2**57
  imag_amplitude=2**57 if large_imaginary else 1/64
  phase[0]=int(real_amplitude*2**scale);phase[n//2]=int(imag_amplitude*2**scale)
  mode_real=1/8 if large_imaginary else 2**56
  mode_imag=2**56 if large_imaginary else 1/32
  phase[1]+=int(mode_real*2**(scale-1));phase[n-1]-=int(mode_real*2**(scale-1))
  phase[n//2+1]+=int(mode_imag*2**(scale-1));phase[n//2-1]+=int(mode_imag*2**(scale-1))
  for direct in ([True,False] if n<=512 else [False]):
   got=canonical_component_decode(phase,scale,direct)
   target=[complex(real_amplitude+mode_real*math.cos(math.pi*pow(5,i,2*n)/n),imag_amplitude+mode_imag*math.cos(math.pi*pow(5,i,2*n)/n)) for i in range(n//2)]
   small_error=max(abs((a.real-b.real) if large_imaginary else (a.imag-b.imag)) for a,b in zip(got,target))
   large_relative_error=max(abs((a.imag-b.imag) if large_imaginary else (a.real-b.real)) for a,b in zip(got,target))/2**57
   assert small_error<1e-11 and large_relative_error<1e-11
   cases.append(dict(n=n,direct=direct,large_imaginary=large_imaginary,small_component_max_error=small_error,large_component_max_relative_error=large_relative_error))
for n in [0,2,3,6]:
 try:canonical_component_decode([0]*n,48)
 except ValueError:pass
 else:raise AssertionError('accepted dimension '+str(n))
args.output.write_text(json.dumps(dict(cases=cases,all_passed=True,invalid_dimensions_rejected=[0,2,3,6],scope='Analytical constant and cosine phase polynomials with extreme component scale separation; not encrypted refresh accuracy.',security_validated=False),indent=2)+'\n');print('Exact phase component split: analytical controls passed')
