# N=65536 scalar Binary-NTT BGV cycle certificate

The corrected scalar bootstrap is a low-to-high BGV circuit:

```text
23-limb ciphertext --BGV modulus drop--> one limb
one limb --centered p scaling + transition--> 23 limbs over plaintext p^2
--16-stage normalized Galois trace, with two drops--> 21 limbs
--bounded odd digit-removal polynomial--> 13 limbs
--exact division by p--> refreshed plaintext-p ciphertext
```

The selected parameters are:

```text
N                              65536
p                              65537
Q                              23 RNS primes, 1402 bits
q0                             one 61-bit RNS prime
operational secret             weight-32 ternary
evaluation-row error           p^2 * CBD(20)
low-to-high gadget             2 digits
Galois trace                   16 keys, 23 digits each
trace modulus drops            after stages 8 and 16
carry interval                 [-23,23]
digit-removal degree           93
bootstrap output               13 limbs
evaluation-key storage         approximately 8.18 GiB
```

The carry interval is deterministic under the accepted input bound. The body
rounding and 32 signed mask roundings contribute at most `16.5`; the scaled
old error contributes at most `6`, leaving strict room below `23`.

The checked worst-case recurrence gives:

```text
projected trace error          < 2^46.00
bootstrap output error         < 2^641.98
13-limb output capacity          2^775.61
two-limb preparation error     < 2^4.17
one-limb refreshed addition      36
multiply-and-drop error        < 2^4.17
one-limb input capacity          2^44.00
cycle contraction               27.41 bits
```

The unlimited-sample attack-cost proxy for the actual `p^2*CBD(20)` source at
the 1402-bit modulus is 133.44 bits. Reserving one bit for reduction accounting
leaves 132.44 bits. This remains a heuristic attack-cost proxy under the
Binary-NTT RLWE premise, not a conventional-RLWE reduction.

Run:

```bash
python3 python/proof/compact_cover_bgv_certificate.py
python3 -m unittest \
  python/proof/test_compact_cover_schedule.py \
  python/proof/test_compact_cover_bgv_certificate.py

sudo singularity exec --bind "$PWD:/work" python/sagemath.sif \
  sage -python /work/python/estimates/regular_cover_bgv_security.py \
  --degree 65536 --qbits 1402 --sigma 13582293620.514343 --mode rough
```

The scalar gate-manifest hash is

```text
9584b90e526fc67ca85c4ea1b6cea004ca64b30b70e4e6609d0961c7e6144843
```

and the certificate hash is

```text
69f97713b99002f8be8fc337b9899bd7e2969b5b27f2345577bd4e3a0cafb3f8
```
