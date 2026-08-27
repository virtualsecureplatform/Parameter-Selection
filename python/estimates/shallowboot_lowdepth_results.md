# Algorithm 3 local parameter screen

The accompanying `shallowboot_lowdepth.py` records Algorithm 3's concrete
tree schedule and screens the sparse LWE input with the vendored
lattice-estimator (`BDGL16`).  The Binary-NTT RLWE component is the paper's
new assumption and is not covered by standard LWE estimation.

| Source parameters | Product windows | Sparse-LWE proxy | Result |
| --- | --- | ---: | --- |
| Paper row: `n=1450, h=29, q=512, sigma=3.2` | `8, 4` | 112.9 bits (`bdd_mitm_hybrid`) | Does not meet this checkout's 128-bit screen |
| Source-screened candidate: `n=1450, h=37, q=512, sigma=3.2` | `16, 4` | 133.3 bits (`bdd_mitm_hybrid`) | Source passes; validate the RLWE modulus/noise schedule |

The paper row has five pointwise product layers (`log2(8)+log2(4)`) and the
larger-source alternative has six.  In both cases the relinearization-free
Binary-NTT product tree performs no NTT/INTT after factors enter the tree.

The second row is not a production parameter set: increasing `h` changes the
noise growth and may require a modulus chain larger than the paper's
`105 -> 50`-bit schedule.  Its PBC choice is `c=3, k=h+3=40`; the `16, 4`
split yields three high-modulus groups and then a three-input low-modulus
tree, retaining only two low-modulus multiplication layers.  It is recorded
to make the discrepancy observable, not to claim a replacement security
proof.

For a **native** `Q≈2^105 -> Q''≈2^50` implementation, the script's
coefficient-domain noise screen gives a conservative low-stage error bound of
about `2^12.6`, against a Boolean plaintext scale of `2^47` (roughly
34.4 bits of margin).  This screen uses `d=4`, base about `2^27`, and the
paper's `sigma'=0.75`; it supports the candidate's noise budget but does not
substitute for a native-Q implementation or a security reduction for the new
Binary-NTT assumption.

An N=4096 experimental RNS run of this PBC shape currently fails the Boolean
sign check even with all sampled encryption errors set to zero.  The failure
occurs after the high-to-low boundary, so this is an implementation/correctness
blocker before a probabilistic noise estimate can be meaningful.  The current
wide-prime harness is therefore only a direct structured-key test harness, not
an executable PBC parameter set.

The direct structured harness also fails once its RLWE bootstrap and boundary
key-switch errors are sampled at the paper's `sigma'=0.75` (even before adding
the source-LWE `sigma=3.2` noise).  The high 105-bit modulus is currently
emulated with two NTT primes.  That representation is incompatible with the
paper's `INTT -> KeySW -> ModSW` boundary: the Binary-NTT secret has independent
NTT-slot representations per residue, whereas a coefficient-domain key switch
requires a single CRT coefficient secret.  Consequently no current low-depth
parameter record is safe to mark usable.  A native ~105-bit NTT implementation,
or an RNS secret/key-switch construction that preserves this compatibility, is
needed before parameter selection can certify one.

One prospective TFHEpp realization is a 128-bit torus with the existing
8-by-16-bit Double Decomposition: embed the native 105-bit level by
`2^23`, so `sigma'=0.75` becomes an encoded torus standard deviation of
`0.75*2^23`.  This retains a single wide secret across `INTT -> KeySW ->
ModSW`; it is the parameter mapping recorded as
`binary_ntt_source_screened_dd`.  The Algorithm-3-specific FullDD tensor
product now collapses to two components using the public quadratic hint and
passes a `lvl3simdparam` multiplication/decryption test without a
relinearization key.  A balanced four-factor tree and DD-native selected/dummy
PBC bucket aggregation also pass.  At the existing `alpha=2^-105`, an
eight-factor (three-layer) high tree passes, while a sixteen-factor/four-layer
tree does not; therefore a boundary remains necessary.  The experimental
128-to-prime boundary preserves each high-group plaintext, but the subsequent
low-prime BFV/QH tree is not yet correct on nontrivial ciphertexts.

A no-boundary experiment using `alpha=2^-120` (absolute sigma 256) is not a
secure substitute: the standard ternary RLWE-as-LWE proxy reports only 112.9
bits.  The secure DD realization must retain enough fresh error and fix the
boundary/low-stage arithmetic.

For the existing `alpha=2^-105` (`sigma=2^23` in the 128-bit torus), the same
ternary RLWE-as-LWE proxy reports **130.2 bits**.  Combined with the 133.3-bit
sparse-source screen, this is the secure parameter target for the DD route.

A DD-high / two-prime-RNS-low prototype was implemented with a CRT-wide
`2^128 -> p0*p1` boundary and CRT-wide BFV scale-and-round.  Each boundary
output decrypts correctly, and the same configuration builds with Intel HEXL
after fixing TFHEpp's orthogonal HEXL/FFT linkage.  However, the first low RNS
multiplication still fails in Debug, both with CRT-wide BFV rounding and with
inverse-scale BGV normalization, including when fresh high noise is disabled.
HEXL therefore supplies the intended per-prime acceleration but does not yet
resolve the post-boundary representation semantics.

The earlier reported RNS PBC bucket-aggregation mismatch was traced to a test
bug: the expected dummy polynomial incorrectly set every coefficient to one
instead of only coefficient zero.  Corrected zero-noise RNS/QH PBC tests pass.
Noisy randomized tests, however, show that the current simplified RLWE PBC key
accumulates too much selector-encryption error.  This persists with 105- and
124-bit high moduli, a two-prime ~105-bit low level, and structured one-hot
blocks as small as 17 entries.  The remaining implementation gap is therefore
the paper-compatible PBC bootstrapping-key representation/noise schedule, not
HEXL's per-prime NTT arithmetic.

Further HEXL/RNS parameter sweeps found no simultaneously secure and correct
setting for the simplified selector key.  Notable points:

- `N=4096, Q≈2^124, sigma'=0.75` has only 108.5 bits under the ternary
  RLWE-as-LWE proxy.
- `N=4096, Q≈2^105, sigma'=0.5` reaches 127.0 bits, narrowly below the
  128-bit target, and fails randomized correctness at trial 5.
- Structured one-hot sources with 32 blocks of size 17 have 130.8 bits of
  entropy and the desired five-layer tree, but noisy randomized bootstraps
  still fail.
- Adding a third high RNS prime and widening the low stage to two primes does
  not repair the selector-noise failure.

This reinforces that the current ordinary-RLWE selector entries are not a
faithful/noise-compatible substitute for the paper's PBC bootstrapping key.

A single-prime fallback (`Q_high ≈ 2^62`, `Q_low ≈ 2^50`) avoids the RNS
representation mismatch, but fails the N=4096 PBC Boolean check even with
sampled errors disabled.  It is therefore too small for the `h=37`, `16 -> 4`
product schedule and is not a usable replacement for the paper's ~105-bit
high modulus.

The paper's D2/QH-SS alternative was also exercised with a 53+52-bit RNS
high modulus, the same coefficient-small secret on both sides of the boundary,
and direct modulus switching.  The current experimental arithmetic fails the
Boolean check even with all sampled errors disabled.  Thus D2 removes the
RNS key-switch incompatibility but does not yet make this implementation a
validated realization of the native-Q construction.
