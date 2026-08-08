# BatchBoot versus block-binary throughput

## Scope

This is a comparison of two *different* ways to reduce the cost of TFHE
bootstrapping:

- **Block binary** accelerates one ordinary Boolean programmable bootstrap.
  It processes each `ell`-wide block of a structured binary LWE secret with one
  external product instead of processing the entries separately.
- **BatchBoot** amortizes one functional bootstrap over a packed batch of
  messages.  Its reported per-message time is total packed-job latency divided
  by the number of slots; it includes the batch accumulation and repacking
  machinery.

Consequently, their quoted `ms` numbers should not be read as a native
head-to-head TFHEpp result.

## 128-bit-screened native TFHEpp measurement (this checkout)

I added `benchmark/bench_batchboot_blockbinary`, then built it separately with
the default parameters and with `USE_BLOCK_BINARY=ON` plus
`USE_SUBSET_KEY=ON`.  Both builds used CMake `Release`, the repository's
`-O3 -march=native` settings, and AVX-512.  Seven warmed repetitions were run
with `taskset -c 0`; key generation and the warm-up operation were paused out
of the timed region.

The security screens for the configurations below are:

| Configuration | Conservative security screen |
| --- | ---: |
| Block binary: source `n=630`, `ell=2`, `alpha=0.0000925119974676756` | 131.9 bits for the block-binary TLWE source (local estimator with its block-aware and MitM checks); 135.2 bits for the dense ternary `n=1024`, `q=2^32` accumulator, conservatively screened at `alpha=2^-24` |
| BatchBoot: source `n=2048`, `q=2^64`, public binary `h=64`, `alpha=2^-25`; dense ternary target `n=2048`, `alpha=2^-51` | min(138.8, 131.5) = 131.5 bits |

The block-binary code uses accumulator fresh noise `alpha=2^-24.8`; screening
it at the larger `2^-24` noise is conservative. Both rows therefore clear a
common 128-bit screening threshold.

| Path timed | Timed latency | Normalized throughput | What it produces |
| --- | ---: | ---: | --- |
| Block-binary `HomNAND` | 5.332 ms | 187.6 Boolean gates/s | one refreshed Boolean TLWE |
| BatchBoot (secure candidate, 2,048 slots) | 112.85 s / batch | 18.15 packed outputs/s (55.10 ms/output) | 2,048 refreshed 1-bit packed outputs |

The block-binary row is the median of seven repetitions; the BatchBoot row is
one security-valid sampled-secret run.

With security matched, the block-binary Boolean path has **10.34x** higher
per-output throughput at the fully occupied 2,048-slot setting. The BatchBoot
job has higher latency because the secure source requires 65 EMP calls and the
target ring degree is 2,048. The block-binary coefficient of variation was
2.4%; the security-valid full BatchBoot run is one timed sample, so it has no
variance estimate. CPU frequency scaling remained enabled, so these are
indicative rather than publication-grade figures.

The 2,048-slot figure was measured directly with a randomly sampled,
uniform fixed-weight-64 binary source secret, then checked by decrypting every
output slot. It is not extrapolated. The full candidate uses six radix-4
digits and has a substantially larger evaluation key. The block-binary side is
the existing `ell=2`, `n=630` Boolean-NAND implementation. The two paths are
built from the same checkout and timed on the same pinned CPU core. The
BatchBoot harness now derives its label and repetition count from the CMake
configuration; further sampled-secret repetitions are needed for a median.

## Published throughput figures

| Scheme / source | Unit timed | Total latency | Normalized latency | Throughput | Key size |
| --- | --- | ---: | ---: | ---: | ---: |
| Block binary (Lee *et al.*, ePrint 2023/958) | one Boolean bootstrap | 6.49 ms | 6.49 ms / message | 154 messages/s | 60 MB |
| Baseline in that paper | one Boolean bootstrap | 10.53 ms | 10.53 ms / message | 95 messages/s | 109 MB |
| BatchBoot Boot2 (supplied Li *et al.* prepublication, Table 4) | 2,048 packed 2-bit messages | 7.32 s | 3.57 ms / message | 280 messages/s | 59.6 MB |

The direct arithmetic comparison is `6.49 / 3.57 = 1.82x`: BatchBoot Boot2
has the lower *amortized* latency in those publications.  This is only a
cross-paper indication, not a benchmark claim: the implementations, CPUs,
threading, parameter sets, ciphertext semantics, and message precision differ.
In particular, the supplied BatchBoot Table 4 was obtained with MOSFHE on a
Xeon Gold 6258R with AVX-512, while the block-binary number is from a TFHE
implementation described in its own paper.

For applications that need the answer from a **single** Boolean ciphertext as
soon as possible, block binary has the relevant latency (6.49 ms in its
publication).  BatchBoot has roughly 7.32 s job latency before its 2,048
answers are available; it wins only when a sufficiently large batch can be
formed and consumed as a batch.

## Work-count interpretation

For the current TFHEpp block-binary parameters, `ell=2` and `n=630`, so the
blind-rotation structure has 315 blocks.  The paper's mechanism therefore
reduces the blind-rotation external-product count from 630 to roughly 315
(plus its separate key-switching improvements).

The security-preserving BatchBoot candidate from `batchboot_results.md`
(`n=2048`, public binary Hamming weight 64, 2,048 slots) has a different work
shape: 65 EMP calls, each with 24 selector-product stages.  Each stage is
applied to all 2,048 slots, so this is 1,560 such stages **per output slot**
(or 3,194,880 slot-level product applications for the full job), before 11
automorphisms and the other batch operations.  The implementation keeps each
stage in Fourier form across the slots, which is why this count cannot be
equated to ordinary blind-rotation external products or converted to
milliseconds.  It is a work-shape description, not a speed estimate.

## Security and implementation compatibility

The existing generic BatchBoot defaults are not a fair secure competitor: the
reassessment estimates its `n=1024`, `h=34` binary source at about 78 bits.
The `n=2048`, `h=64` source candidate reaches about 138.8 bits, with the
dense level-2 target at about 131.5 bits.  It is a candidate only; it still
needs a Monte-Carlo failure experiment and a full evaluation-key/multi-target
analysis.

The supplied TFHEpp block-binary configuration has a block-binary LWE source
(`n=630`, `ell=2`) but its `lvl1param` has module rank `k=2`.  Generic
`BatchRingSwitchSecret` and `BatchRingSwitch` currently require an RLWE source
with `k=1` (`include/tfhe/batchboot.hpp`).  Thus the current block-binary
parameter family cannot simply be selected as the BatchBoot source; a combined
construction needs a new compatible parameter and key-switch design, followed
by a new security/noise analysis.

## A fair native TFHEpp benchmark

Do not time the correctness tests: `blockbinary_blindrotatefft.cpp`,
`blockbinary_nand.cpp`, and `batchboot.cpp` contain no throughput harness.
For a meaningful comparison, add a benchmark that:

1. generates keys outside the timed region, warms up FFT/NTT state, pins the
   same CPU cores, and fixes the thread count;
2. measures block-binary Boolean PBS end-to-end (blind rotation plus key
   switch) and reports median/p95 latency and messages/s;
3. measures BatchBoot's complete packed operation (ring switch, accumulation,
   automorphisms/repacking as applicable), reports both batch latency and
   `batch latency / slots` plus aggregate messages/s;
4. uses configurations with the same required precision and at least 128-bit
   security under the same estimator assumptions; and
5. reports bootstrapping/evaluation-key size and empirical decryption failure
   rate alongside timing.

The native harness now provides a 128-bit-screened full-occupancy measurement.
It still needs repeated sampled-secret runs, an occupancy sweep (8, 64, ...,
2,048 slots), and a failure-rate experiment before a deployment decision. On
this checkout, use block binary for latency-sensitive single-bit PBS and for
better measured per-output throughput at full BatchBoot occupancy. Do not use
the current level-1 BatchBoot source in either comparison because it misses the
stated security target.

## Sources

- C. Lee *et al.*, [Faster TFHE Bootstrapping with Block Binary Keys](https://eprint.iacr.org/2023/958), ePrint 2023/958.
- Z. Li *et al.*, supplied `sec26_prepub_li-zhihao.pdf`, Table 4.
