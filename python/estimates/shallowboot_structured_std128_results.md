# Structured shallow bootstrap security screen

This records the runnable TFHEpp candidate in
`TFHEpp/include/params/shallowboot.hpp`, plus an estimated shorter
post-key-switch blind-rotation candidate:

- source: structured binary, `n = 1024`, `q = 2^9`, one-hot `h = 64`,
  Gaussian `sigma = 3.2`;
- accumulator: TFHEpp's existing ternary `N = 1024`, 32-bit torus ring with
  normalized error `alpha = 2^-25` (an RLWE-as-LWE proxy).

The source-return key switch is outside this model. The current generic TFHEpp
key switch does not preserve enough margin for the `q = 512` representation,
so the shorter candidate is not enabled in TFHEpp yet.

## Reproduction

From `Parameter-Selection/python`:

```sh
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -e DOT_SAGE=/tmp/.sage \
  -v "$PWD:/work" -w /work sagemath/sagemath:latest \
  sage -python estimates/shallowboot_structured_std128.py
```

The recorded execution used the vendored lattice-estimator commit `3e48ef4`,
Sage container image
`sha256:8436018569d7738063789974afce157f59d60a5832710a5c2581218ec673dd34`,
and the `BDGL16` reduction-cost model.

## Results

| Instance | Distribution used by estimator | Minimum classical work factor |
| --- | --- | ---: |
| Source | Uniform weight-64 sparse binary proxy | 170.3 bits (`bdd_mitm_hybrid`) |
| Accumulator | Ternary RLWE-as-LWE proxy | 129.6 bits (`dual_hybrid`) |

The parameterized screen also evaluated shorter one-hot source candidates:

| Post-switch BR secret | Source proxy security | PBS log2 variance | Boolean log2 failure bound | Decision |
| --- | ---: | ---: | ---: | --- |
| `n=1024, h=64, q=512` | 170.3 | -17.36 | -477.73 | Baseline |
| `n=576, h=64, q=512` | 130.9 | -18.14 | -820.42 | Candidate |
| `n=512, h=64, q=512` | 121.6 | -18.30 | -914.53 | Reject: security |

The `n=576` row is the smallest tested whole-bucket candidate that passes
the 128-bit proxy screen. Its 64 buckets contain nine source coordinates each,
so it retains 64 external products while reducing the linear bootstrapping-key
work from 1,024 to 576 entries.

The accumulator therefore clears a 128-bit LWE-proxy screen by 1.6 bits.  The
source's 170.3-bit result is **not** a proof for the structured one-hot
distribution: it is the standard uniform sparse-binary proxy used in the
paper's comparison.  The actual one-hot secret support has `16^64 = 2^256`
elements, and the paper's Table 4 relies on an entropy-equivalence heuristic
for its concrete-security claim.  Treat this as a screened research parameter,
not a production security proof.
