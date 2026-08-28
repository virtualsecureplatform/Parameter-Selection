# Regular-cover BGV preliminary parameter screen

## Scope

This screen evaluates only quantities fixed by the formal regular-cover
reduction:

- cover cardinality and storage;
- dense-binary challenge-masking rows;
- pivot unit density and rejection failure;
- ordinary base-RLWE sample count; and
- the union-bound loss over cover components.

It does not yet reproduce the complete native BGV bootstrap noise recurrence.
Consequently, the noise contraction values passed to the examples below are
illustrative rather than a correctness certificate.

## Parameters from `Bootstrapping_BGV_BFV/Crypto/Params.m`

The checked implementation example has approximately

```text
base degree N = 512
log2(q)       = 900
error sigma   = sqrt(10)
```

Taking the worst-case full Galois cover `|Gamma|=512`, 128-bit statistical
challenge masking requires

```text
cover log2 cardinality = 235,929,600
public-key rows         = 471,859,454
ciphertext storage      = 471,859,200 bits
```

With 10,000 additional evaluation rows and 128 pivot trials, the reduction
uses approximately

```text
241,597,225,984 ordinary base-RLWE samples.
```

This is polynomial but impractical.

More importantly, flattening the ordinary Binary-NTT RLWE source to the same
LWE attack-cost proxy used elsewhere in this repository gives only

```text
best estimated attack cost: 2^11.68
```

for `N=512`, `log2(q)=900`, and `sigma=sqrt(10)`. Thus the example parameters
are not a secure instantiation of the theorem.

## Security-dimension comparison

Using the same modulus and error width with an unlimited-sample attack proxy:

```text
N = 32,768  -> approximately 89.94 bits: FAIL
N = 65,536  -> approximately 237.40 bits: PASS for a 128-bit target
```

The estimates are attack-cost proxies, not a reduction from Binary-NTT RLWE
to unstructured LWE.

For the conservative full cover `|Gamma|=65,536`, the proof-induced costs at
`N=65,536`, `log2(q)=900` are

```text
cover log2 cardinality = 3,865,470,566,400
public-key rows         = 7,730,941,133,054
ciphertext storage      = 7,730,941,132,800 bits
```

The source-sample count, including one million evaluation rows and pivot
trials, is approximately

```text
506,655,023,640,215,552.
```

Again, these values are polynomial and satisfy the proof's statistical terms,
but are far outside practical use.

## Commands

```bash
python3 python/proof/regular_cover_bgv_screen.py \
  --degree 512 --group-size 512 --qbits 900 --evaluation-rows 10000
```

```bash
singularity exec --bind "$PWD:/work" python/sagemath.sif \
  sage -python /work/python/estimates/regular_cover_bgv_security.py \
  --degree 65536 --qbits 900 --sigma 3.1622776601683795 --mode rough
```

## Remaining correctness input

Before any parameter is certified, the selected native bootstrap must provide:

- an exact per-component input-noise bound;
- an exact fresh-output-noise bound;
- a per-component failure probability;
- the exact number of switch, boot, and multiplication rows; and
- the modulus-switching/RNS schedule.

The regular cover preserves the per-component deterministic recurrence and
subtracts `log2(|Gamma|)` bits from the probabilistic failure margin. Once the
native values are supplied, the screen can evaluate the complete cover
candidate without another proof-side assumption.
