# Algorithm 3 PBC/QH selector-noise search

The estimator follows direct `RLWE_lsb` PBC aggregation, the balanced selector
product tree, and every modulus switch. It establishes the noise of the LSB
selector product only; it does not establish the final MSB TFHE output.

## Paper-sized row

Section 7.2 gives `N=4096`, an initial modulus near 105 bits, a 50-bit level,
and a 12-bit modulus near the output. Footnote 11 additionally requires a
modulus `Q'''` strictly between the 50- and 12-bit levels. Search it with:

```sh
python3 python/estimates/shallowboot_noise.py \
  --lwe-h 29 --modulus-chain 105 50 24 12 --switch-after 3 4 5 \
  --search-middle 13:49
```

No integer bit-size for the middle modulus passes the average-case 128-bit
noise screen. This applies to the direct-RLWE PBC interpretation implemented
in TFHEpp. It does not prove every realization of the paper fails: the paper
does not give the concrete `Q'''` or a complete concrete noise calculation.

## Workable screened candidate

The local sparse-source proxy needs `h=37`, producing 40 factors and six
product layers. A joint noise search found:

```text
source:       n=1450, h=37, q=2^9, sigma=3.2
ring:         N=8192, RLWE sigma'=0.75, plaintext t=4
moduli:       2^151 -> 2^69 -> 2^52 -> 2^36
switch after: product layers 3, 4, and 5
PBC:          c=3, k=40
output key:   n=1450, h=60, Q'=2^15, sigma=3.2, basebit=1
```

Reproduce it with:

```sh
python3 python/estimates/shallowboot_noise.py \
  --ring-n 8192 --modulus-chain 151 69 52 36 \
  --switch-after 3 4 5
```

The analytic average-case screen reports 19.7 bits minimum pre-switch wrap
headroom and a post-conversion/LUT standard deviation below `2^-7` after
scaling to `q=512`. A secure-size TFHEpp regression additionally passes with a
dense ternary QH secret and aggregate PBC-equivalent Gaussian noise.

The 15-bit output key-switch value is not an `N=8192` NTT modulus because a splitting
prime must satisfy `2N | Q-1`. It is instead used by the
RLWE-to-LWE/output conversion. TFHEpp uses concrete NTT-prime products with bit
lengths 151, 69, 52, and 36. The RNS bases are nested prefixes so the exact BGV
modulus-switch correction is applicable.

The source sparse-LWE proxy was previously measured at 133.3 classical bits.
The exact `N=8192,Q=2^151,sigma'=0.75` QH-SS-RLWE-as-LWE proxy reports 188.7
classical bits (`dual_hybrid`), while the sparse source proxy reports 133.3
bits. The `arora_gb+guessing` and plain `dual` estimator routines could not
handle this secret distribution; successful uSVP, BDD, and dual-hybrid
estimates were all above 188 bits. QH-SS-RLWE itself remains the paper's
nonstandard assumption, so these are standard-LWE proxy results rather than a
reduction for that assumption. The candidate is also materially more expensive
than the paper's `N=4096` estimate.

The refreshed/output sparse key is intentionally different from the input key.
At the 15-bit key-switch modulus, even `n=1150,h=60,sigma=3.2` reports 138.6
bits; the executable conservatively retains dimension 1450.
Running the same selector-noise model with `--lwe-h 60` also passes the
six-layer chain, so the refreshed key can be bootstrapped again with a
63-bucket PBC schedule.
The paper-like 12-bit `h=37` output switch failed empirically from key-switch
noise, and the 13-bit `h=37` and `h=43` screens report 120.1 and 131.4 bits
respectively; the latter was still marginal in the two-branch runtime test.

## Encoding correction used by TFHEpp

The current ePrint PDF (SHA-256
`0ba34c2052e2fc717e17271ac04058096075bfd21e8922b7009c9328119eef12`)
has no published code or erratum located as of 2026-08-24. Its pseudocode
omits the use of `ACC`, but the intended encoding path can be completed as
follows:

1. compute and modulus-switch the `RLWE_lsb` selector product in plaintext
   ring `R_4`;
2. multiply it by `4^-1 mod Q`, mapping `m+4e` to an MSB `R_4` encoding
   without amplifying error;
3. multiply by the unscaled rotated sign LUT; and
4. extract and modulus-switch the resulting ordinary MSB LWE ciphertext.

The modular-inverse identity is the missing linear conversion. TFHEpp tests it
at `N=8192` with aggregate PBC-equivalent noise, the full four-stage tree, a
dense LUT, sample extraction, and output modulus switching.

## Full-key execution

`test/tfhe/shallowboot_secure_full.cpp` constructs either the complete
`n=1450,h=37,c=3,k=40` PBC key or the structured `n=1024,h=64` key, plus a
noisy 15-bit ring-to-LWE key-switch key.
The one-thread HEXL run on this host reports:

```text
PBC entries:             4290
estimated PBC key:       3.142 GiB
key generation:          4.84 s
bootstrap (HEXL, 1T):    0.793 s
bootstrap (HEXL, 8T):    0.230 s
bootstrap (HEXL):        0.179 s general (16T), 0.125 s structured (16T)
seeded structured:       0.122 s (8T)
peak resident memory:    4.44 GiB
output phases:           valid noisy ±128 branches in both profiles
```

The optimized one-thread scalar executable reports 1.475 s, so HEXL gives a
1.86x end-to-end speedup before thread-level parallelism. The nested-RNS
modulus drop, fused/cached PBC aggregation, bucket parallelism, coefficient
parallelism, and parallel output key switch reduce the earlier 1.676 s HEXL
runtime by 8.8x at 16 threads.

The existing one-thread block-binary benchmark reports 4.639 ms, making this
general implementation about 38.6 times slower and the seeded structured
profile about 26.4 times slower at their best measured thread
count. The 3+ GiB evaluation-key stream is already memory-bandwidth-bound by
16 threads; seeded key compression and a digit/tower-major contiguous key
layout are the main remaining optimization targets.

The final tower-drop implementation follows OpenFHE `ed361af`'s
`DCRTPolyImpl::ModReduce`: compute the plaintext correction in the last tower,
center-map it to surviving towers, multiply by the dropped-prime inverse, and
drop the tower. This removed CRT interpolation entirely. OpenFHE's bounded
thread controls and digit/tower-major key-switch loops also motivated the
bucket, coefficient, and output-key-switch parallel reductions used here.

Run both security proxies in the repository's Sage environment with:

```sh
sage -python python/estimates/shallowboot_lowdepth.py \
  --lwe-h 37 --ring-n 8192 --initial-qbits 151 \
  --post-switch-qbits 69 --windows 8 2 2 2 --ring-proxy-security
```

The refreshed key screen is reproduced with `--lwe-h 50 --lwe-qbits 14`
(without `--ring-proxy-security`).

## N=4096 search

A structured `n=1024,h=64` source makes the six-layer selector tree much less
noisy than the general PBC source. An abstract five-boundary schedule at
`Q=2^105` passes the average-case screen, and the ring proxy is 129.3 bits for
`sigma'=0.75`. Concrete `N=4096` NTT primes are at least about 16 bits, however,
so a nested six-level chain needs a 107-bit initial product. At 107 bits the
ring proxy is 126.7 bits for `sigma'=0.75`, 127.6 bits for `sigma'=1.0`, and
128.5 bits for `sigma'=1.5`.

The concrete `sigma'=1.5` chain with prefix-product bit lengths
`107 -> 90 -> 73 -> 56 -> 40 -> 25` was implemented and tested. Its full
aggregate-noise outputs were phases 212 and 310 instead of the target
neighborhoods around 384 and 128. Thus the average-case model understates the
correlated final-layer noise, while the narrow security margin leaves no room
to enlarge the modulus. This `N=4096` route is rejected under the current
QH-SS-RLWE proxy and correctness requirements.

## Additional tuning trials

- The structured `N=8192,n=1024,h=64` profile reduces PBC entries from 4290
  to 1088 and the estimated PBC key from 3.14 GiB to 0.797 GiB. Its best
  verified 16-thread runtime is 125 ms.
- BLAKE3-seeded `a` polynomials reduce that PBC key again to 0.398 GiB. The
  full noisy branches pass at phases 371 and 127. Regeneration is effectively
  hidden by reduced memory traffic: 122 ms at 8 threads versus 129 ms for the
  unseeded structured profile, while 16 seeded threads are slightly slower
  (126 ms) from bandwidth/parallel overhead.
- A balanced base-4 output key switch was implemented and passes an exact
  zero-noise algebra regression, but the secure noisy run landed near the
  decision boundary and was rejected. Binary decomposition remains active.
- Keeping the last two products at the same 52-bit level removed one modulus
  boundary analytically, but failed both general and structured aggregate-
  noise executions. The verified `151 -> 69 -> 52 -> 36` chain is retained.
- Shrinking the refreshed LWE dimension changed the sampled raw-key-switch
  error enough to produce marginal/failing branches. The executable retains
  the robust `n=1450,h=60,Q'=2^15` refreshed key.

This checkout built `python/sagemath.sif` using the supplied definition and
Singularity's rootless `--fakeroot` mode. The estimator used the vendored
lattice-estimator and the `BDGL16` reduction-cost model.
