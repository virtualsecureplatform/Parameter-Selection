# N=65536 scalar Binary-NTT BGV cycle certificate

The corrected scalar bootstrap is a low-to-high BGV circuit:

```text
20-limb ciphertext --BGV modulus drop--> one limb
one limb --centered p scaling + transition--> 20 limbs over plaintext p^2
--16-stage normalized Galois trace, with two drops--> 18 limbs
--bounded odd digit-removal polynomial--> 10 limbs
--exact division by p--> refreshed plaintext-p ciphertext
```

The selected parameters are:

```text
N                              65536
p                              65537
Q                              20 RNS primes, 1219 bits
q0                             one 61-bit RNS prime
operational secret             weight-32 ternary
evaluation-row error           p^2 * CBD(20)
low-to-high gadget             2 digits
Galois trace                   16 keys, 23 digits each
trace modulus drops            after stages 8 and 16
carry interval                 [-23,23]
digit-removal degree           93
bootstrap output               10 limbs
evaluation-key storage         approximately 7.09 GiB
```

The carry interval is deterministic under the accepted input bound. The body
rounding and 32 signed mask roundings contribute at most `16.5`; the scaled
old error contributes at most `6`, leaving strict room below `23`.

The checked worst-case recurrence gives:

```text
projected trace error          < 2^37.97
bootstrap output error         < 2^154.01
10-limb output capacity          2^592.78
two-limb preparation error     < 2^4.17
one-limb refreshed addition      36
multiply-and-drop error        < 2^4.17
one-limb input capacity          2^44.00
cycle contraction               27.41 bits
```

The unlimited-sample attack-cost proxy for the `p^2*CBD(20)` evaluation-row
source at the 1219-bit modulus is 162.94 bits. The public affine context is a
separate Binary-NTT row whose error is the weight-32 signed-ternary operational
secret. Flattening its single ring row to 65536 scalar samples and matching its
coefficient variance gives a 155.93-bit proxy. The full-source hybrid charges
one evaluation term and one context term, while the auxiliary-context triangle
charges the context term a second time. Their conservative sum, after unit-pivot
conditioning, is 154.92 bits. The exact product-ring unit failure is
below `2^-40.61`, so this conditioning changes the bit estimate by less than
`10^-12`. These remain heuristic attack-cost proxies
under the Binary-NTT RLWE premise, not a conventional-RLWE reduction; the
fixed-weight error is only variance-matched to the estimator's Gaussian model.

Run:

```bash
python3 python/proof/compact_cover_bgv_certificate.py
python3 -m unittest \
  python/proof/test_compact_cover_schedule.py \
  python/proof/test_compact_cover_bgv_certificate.py

sudo singularity exec --bind "$PWD:/work" python/sagemath.sif \
  sage -python /work/python/estimates/regular_cover_bgv_security.py \
  --degree 65536 --qbits 1219 --sigma 13582293620.514343 --mode rough

sudo singularity exec --bind "$PWD:/work" python/sagemath.sif \
  sage -python /work/python/estimates/regular_cover_bgv_security.py \
  --degree 65536 --qbits 1219 --sigma 0.02209708691207961 \
  --samples 65536 --mode rough
```

The scalar gate-manifest hash is

```text
9584b90e526fc67ca85c4ea1b6cea004ca64b30b70e4e6609d0961c7e6144843
```

and the certificate hash is

```text
155a50174b509e170a649674d38c55c450f475e43ae56d32a8865c3459eaab7e
```
