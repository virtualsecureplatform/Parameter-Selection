# BatchBoot reassessment (TFHEpp implementation)

This is a reproducible screening report for
`TFHEpp/include/tfhe/batchboot.hpp`, run with the local
`lattice-estimator` and reduction-cost model `BDGL16`.

## Method

- The sparse input secret and accumulator secret are estimated separately.
  BatchBoot's EMP selector and automorphism keys are encryptions under the
  accumulator secret, whereas the input ciphertext is protected by the sparse
  source secret.
- The sparse binary source distribution is represented as
  `SparseTernary(h, 0, n)`.  This is the exact `0/1` distribution required by
  `BatchBootKeyGen`, with a public fixed Hamming weight `h`.
- Each instance uses `m=n` samples.  Security below is the minimum `log2(rop)`
  over the estimator's reported attacks, including `bdd_mitm_hybrid`.
- Noise is a normalized conservative variance screen.  It counts the concrete
  implementation's four selector external products per stored radix-4 digit.
  It excludes FFT floating-point error, so it is not a correctness proof.

## Results

| Source / target instance | Weakest LWE estimate | Result |
| --- | ---: | --- |
| Current level-1 source: `n=1024`, `q=2^32`, binary `h=34`, `alpha=2^-25` | 78.0 bits (`bdd_mitm_hybrid`) | Reject |
| Current level-1 target: `n=1024`, `q=2^32`, ternary, `alpha=2^-25` | 129.6 bits (`dual_hybrid`) | Passes alone |
| Level-2 source: `n=2048`, `q=2^64`, binary `h=34`, `alpha=2^-25` | 100.5 bits (`bdd_mitm_hybrid`) | Reject |
| Same source with `alpha=2^-21` | 106.0 bits (`bdd_mitm_hybrid`) | Reject |
| Same source with `alpha=2^-15` | 116.5 bits (`bdd_mitm_hybrid`) | Reject; ModSwitch sigma is 0.125 bins |
| Level-2 source: `n=2048`, `q=2^64`, binary `h=64`, `alpha=2^-25` | 138.8 bits (`bdd_mitm_hybrid`) | Pass |
| Level-2 target: `n=2048`, `q=2^64`, ternary, `alpha=2^-51` | 131.5 bits (`dual_hybrid`) | Pass |

For the recommended `h=64` candidate at 2,048 slots, `BatchBoot.py` reports
65 EMP calls, 24 external products per EMP, output standard deviation at most
`2^-19.28` for a 3-bit output, and source ModSwitch standard deviation
`1.22e-4` exponent bins.  The same bound leaves ample margin for an 8-bit
decoded output; rerun the script with `--output-bits 8` when setting the exact
application precision.

## Conclusion

Do not deploy the generic default level-1 BatchBoot combination: it fails both
the source-key security target and the full-packing 3-bit noise screen.

A candidate direction is a **separate** sparse binary source parameter with
`n=2048`, `q=2^64`, fixed public weight 64, and `alpha=2^-25`, combined with a
dense ternary level-2 accumulator target (`n=2048`, `q=2^64`, `l=4`,
`Bg=2^9`, `alpha=2^-51`).  This changes BatchBoot work from 35 to 65 EMP calls
relative to `h=34`, but the measured estimator margin is about 10.8 bits above
128.  It still needs an implementation-level parameter struct, a fresh
Monte-Carlo failure experiment, and a review of the desired multi-target and
evaluation-key sample model before being labelled a production parameter set.
