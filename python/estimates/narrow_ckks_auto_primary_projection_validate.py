#!/usr/bin/env python3
"""Validate primary8 auto-control ciphertexts using independent integer phases.

Reads the original private diagnostic secret prefix without serializing secrets.
Exit 2 indicates a verified error exceeding the fixed 23-bit bound.
Usage: SCRIPT PRIVATE_STATE PUBLIC_CIPHERTEXTS CPP_REPORT
"""
import pathlib,sys,struct,json,math
from narrow_ckks_prefix_arithmetic import PhaseDecoder,canonical_decode
state,outputs,report=map(pathlib.Path,sys.argv[1:4])
with state.open('rb') as f:
 marker=f.read(1);assert marker in (b'\x00',b'\x01');endian='<' if marker==b'\x01' else '>'
 meta=struct.unpack(endian+'10Q',f.read(80));n=meta[1]
 assert meta==(0x4155583850523031,n,85,44,385,48,37,108,8,448) and n in (512,16384)
 length=struct.unpack(endian+'Q',f.read(8))[0];assert length==7 and f.read(length)==b'sigmoid' and f.read(1)==b'\x01'
 words=struct.unpack(endian+str(n*7)+'Q',f.read(n*56))
 dense=[sum(words[i*7+j]<<(64*j) for j in range(7)) for i in range(n)]
 dense=[x if x<2**447 else x-2**448 for x in dense];assert all(x in (-1,0,1) for x in dense)
decoder=PhaseDecoder(dense);samples=[];max_old=max_new=max_delta=max_input=0
with outputs.open('rb') as f:
 marker=f.read(1);assert marker in (b'\x00',b'\x01');endian='<' if marker==b'\x01' else '>'
 header=struct.unpack(endian+'5Q',f.read(40))
 auxiliary=header[0]==0x4155544f41383031
 q,delta=(348,48) if auxiliary else (356,49)
 assert header==((0x4155544f41383031 if auxiliary else 0x4155544f50383031),n,q,delta,3)
 if auxiliary:assert struct.unpack(endian+'Q',f.read(8))==(8,)
 def polynomial():
  data=f.read(n*56);assert len(data)==n*56
  words=struct.unpack(endian+str(n*7)+'Q',data)
  return [sum(words[i*7+j]<<(64*j) for j in range(7)) for i in range(n)]
 for sample in range(3):
  phases=[]
  if auxiliary:
   values=struct.unpack(endian+str(n)+'d',f.read(n*8))
   target=[complex(values[2*i],values[2*i+1]) for i in range(n//2)]
  else:target=[0j]*(n//2)
  for kind in ['input','old','projected']:
   a,b=polynomial(),polynomial();phase=decoder.phase(a,b,q);phases.append(phase)
   if kind=='input' and not auxiliary:assert all(x==0 for x in phase)
  input_error=max(abs(x-y) for x,y in zip(canonical_decode(phases[0],delta),target))
  old=max(abs(x-y.conjugate()) for x,y in zip(canonical_decode(phases[1],delta),target))
  new=max(abs(x-y.conjugate()) for x,y in zip(canonical_decode(phases[2],delta),target))
  difference=max(map(abs,canonical_decode([b-a for a,b in zip(phases[1],phases[2])],delta)))
  assert all(math.isfinite(x) for x in [input_error,old,new,difference])
  max_old=max(max_old,old);max_new=max(max_new,new);max_delta=max(max_delta,difference);max_input=max(max_input,input_error)
  samples.append({'sample':sample,'input_integer_phase_zero':not auxiliary,'old_max_slot_error':old,'new_max_slot_error':new,'max_slot_difference':difference,'input_max_slot_error':input_error})
  print(json.dumps(samples[-1]),flush=True)
 assert f.read(1)==b''
r=json.loads(report.read_text())
for name,value in [('old_max_slot_error',max_old),('new_max_slot_error',max_new),('max_slot_difference',max_delta)]+([('input_max_slot_error',max_input)] if auxiliary else []):assert math.isclose(r[name],value,rel_tol=1e-6,abs_tol=1e-12),(name,r[name],value)
e={'scope':'Independent exact integer phases from public diagnostic ciphertexts and original private retained ternary key; N512 scalar convolution or N16384 bounded 15-bit NTT. No key material serialized. Tests primary8 direct switching or the auxiliary8 path on nonzero encrypted inputs; not full packed reconstruction, whole refresh, or security.','n':n,'samples':samples,'old_max_slot_error':max_old,'new_max_slot_error':max_new,'max_slot_difference':max_delta,'fixed_23bit_passed':max_input<2**-23 and max_old<2**-23 and max_new<2**-23,'input_max_slot_error':max_input,'log_q':q,'log_delta':delta,'auxiliary_bits':8 if auxiliary else 0,'cpp_metrics_match':True,'fresh_keys_generated':False}
(outputs.parent/'independent-validation.json').write_text(json.dumps(e,indent=2)+'\n');print(json.dumps(e),flush=True)
if not e['fixed_23bit_passed']:sys.exit(2)
