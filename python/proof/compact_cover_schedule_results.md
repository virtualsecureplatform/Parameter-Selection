# N=65536 thin-BGV compact-cover schedule

The schedule extractor mirrors the non-cyclic power-of-two thin-bootstrap
branch in `Bootstrapping_BGV_BFV`.

For

```text
N = 65536
m = 131072
p = 65537
```

it obtains:

```text
Frobenius order d          = 2
plaintext slots            = 32768
hypercube dimensions       = (16384, 2)
hypercube generators       = (5, -1)
baby dimensions            = (2, 91)
giant dimensions           = (1, 181)
baby product               = 182
giant product              = 181
unique switch automorphisms= 362
generated subgroup size    = 65536
peak live ciphertexts      = 368
```

Thus the switch-key exponents generate the complete Galois group even though
the baby-step/giant-step evaluator keeps only a few hundred ciphertexts live at
once.

With one 64-bit RNS limb, a literal full-cover ciphertext is 64 GiB, whereas
the live scheduled frontier is approximately 0.36 GiB. With fifteen limbs,
the corresponding figures are approximately 960 GiB and 5.39 GiB.

This resolves the scheduling question: the native circuit has a narrow live
frontier despite generating the full group. The remaining problem is not
schedule extraction or memory accounting. It is the cryptographic theorem
that securely rebases and merges partial-cover frontiers while all transitions
share one ordinary Binary-NTT RLWE source and ultimately return to the same
advertised key.

Run:

```bash
python3 python/proof/compact_cover_schedule.py
python3 -m unittest python/proof/test_compact_cover_schedule.py
```
