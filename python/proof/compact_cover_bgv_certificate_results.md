# N=65536 scalar compact-cover BGV certificate

The scalar-only specialization uses the identity

```text
p*(m + p*e) + p^2*e_ks = p*m + p^2*(e + e_ks).
```

One balanced full-modulus transition followed by exact public division by `p`
therefore refreshes a plaintext-`p` scalar ciphertext without coefficient to
slot conversion or probabilistic digit extraction.

The selected profile is:

```text
N                         65536
p                         65537
RNS limbs                 15
log2(Q)                   914.429
operational secret        weight-32 ternary
error                     CBD(20), sigma=sqrt(10), support [-20,20]
full-modulus gadget rows  5
gadget digit width        183 bits
bootstrap key             75 MiB plus header
```

Every RNS prime is one modulo both `2N` and `p^2`. The deterministic bounded
error calculation gives an output coefficient magnitude below approximately
`2^220.65`, versus `Q/2 > 2^913`. In normalized variance notation this is
`log2(V) <= -1387.57`, below the accepted input value `-900`; the refresh
contraction margin is about 487 bits.

The unlimited-sample LWE attack-cost proxy for the Binary-NTT source with the
actual `p^2*CBD(20)` evaluation-row phase error reports 243.53 bits at this
modulus. Reserving one bit for reduction accounting leaves 242.53 bits. This
is an attack-cost proxy under the Binary-NTT RLWE premise, not a
conventional-RLWE reduction.

Run:

```bash
python3 python/proof/compact_cover_bgv_certificate.py
python3 -m unittest \
  python/proof/test_compact_cover_schedule.py \
  python/proof/test_compact_cover_bgv_certificate.py

sudo singularity exec --bind "$PWD:/work" python/sagemath.sif \
  sage -python /work/python/estimates/regular_cover_bgv_security.py \
  --degree 65536 --qbits 915 --sigma 13582293620.514343 --mode rough
```

The canonical scalar gate manifest hash is

```text
209b8826908383c92bce2ea41f27eda9febbf250e69e2c5541f5c18b76e454f0
```

and the certificate hash is

```text
143be2a21f31cd5bb2bb3d6359b81f1103a05329639fd42152c3503669dd39b1
```
