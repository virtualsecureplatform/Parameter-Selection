# TFHEpp quadratic-KDM proof screens

`lvl5boot_renyi.py` checks the first necessary numerical margin of the
fixed-weight, equal-covariance Rényi certificate formalized by
`FormalProof4FHE.RLWE.RankOneHNFLossinessRenyi`.

Run the default TFHEpp `lvl5bootparam` screen from the repository root:

```bash
python3 python/proof/lvl5boot_renyi.py --self-test
```

Use `--json` to obtain the exact rational certificate fields in a
machine-readable form.

## Checked quantity

For an even negacyclic ring degree and a secret uniform on

```text
T(n,w) = {s in {-1,0,1}^n : ||s||_0 = w},
```

the formally proved square covariance gives

```text
E[||S^2||_2^2]
  = 2 n^2 p2 + n (p - 3 p2),
p  = w/n,
p2 = w(w-1)/(n(n-1)).
```

For a spherical Gaussian terminal error with coefficient standard deviation
`sigma`, one scalar gadget row contributes

```text
(g/sigma)^2 E[||S^2||_2^2]
```

to the whitened quadratic energy. The hidden small-ratio denominator scales
both the quadratic message and terminal error and therefore cancels. The
linear masked-ratio energy, the other gadget rows, and every concentration or
bad-descriptor term are nonnegative in the uniform symmetric-prior
specialization.

The first-order Rényi margin can only be positive if

```text
E[F] < 2 log |T(n,w)|.
```

The program checks an exact, stronger obstruction. If `b` is the bit length
of `|T(n,w)|`, then `2 log |T(n,w)| < 2b`. It compares the top-row rational
energy directly with this integer upper bound, so floating-point logarithms
are used only for readable output.

## Result for `lvl5bootparam`

The implementation parameters are

```text
n = 32768, q = 2^640, w = 96, sigma = 2^33,
l = 35, Bg = 2^18.
```

The exact fixed-weight moment is

```text
E[||S^2||_2^2] = 600806592 / 32767
                 ~= 18335.7216712.
```

The first gadget row has `g = 2^622`. Its contribution alone has
`log2(E[F]) ~= 1192.1624`, whereas the complete twice-entropy budget is only
`1438.3107` as an ordinary real number (`log2 ~= 10.4902`). Thus the mean
condition misses by a factor of approximately `2^1181.67`; 33 of the 35
ordinary gadget rows individually exceed the exact entropy upper bound.

Even before the other rows and concentration penalties, the top row would
need `sigma > 2^623.836`, or normalized error greater than approximately
`2^-16.164`. That is incompatible with the narrow BFV correctness regime.

## Scope

This is a rigorous rejection of the direct uniform-fixed-weight,
equal-covariance Gaussian Rényi certificate applied to the ordinary rows from
which the current gadget key is encoded. It does **not** prove that quadratic
KDM is insecure, and it does not rule out a sharper analysis of the final
Double-Decomposition/FFT encoding, a different computational reduction, a
different reference channel, or a leakage-conditioned argument whose
posterior is analyzed by another method. Data processing only says the final
encoding leaks no more than the ordinary-row view; failure of an upstream
sufficient bound does not reverse that inequality. An exact reconstruction
or Rényi-invariance theorem for the encoding would transfer this obstruction.
The power-of-two `Theta` collision theorem is also a separate route.

TFHEpp first samples `l` ordinary RLWE rows and then applies the `lbar` Double
Decomposition and FFT representation. The screen counts those `l = 35`
independent rows; it does not incorrectly treat the `l*lbar = 1400` public
digits as independent error samples.

## Lean certificate

The matching module
`FormalProof4FHE.RLWE.TFHEppLvl5BootRenyiObstruction` kernel-checks the support
bound, specialized fixed-weight moment, top-row energy inequality, and the
resulting nonpositive first-order margin. The executable checker additionally
reports the sharp logarithms, all 35 row exponents, and threshold estimates.

## Symmetry-reduced Gaussian-cluster follow-up

`lvl5boot_gaussian_cluster.py` tests the local Gibbs/mixture idea without
enumerating the enormous secret support. It reads `lvl5bootparam` and the
multi-limb modular-Gaussian path directly from the sibling TFHEpp checkout:

```bash
python3 python/proof/lvl5boot_gaussian_cluster.py --self-test
python3 python/proof/lvl5boot_gaussian_cluster.py --json
```

Use `--tfhepp-root PATH` when TFHEpp is elsewhere. The JSON output includes
SHA-256 hashes of both source headers, so a saved result identifies the exact
implementation it screened.

For a fixed signed weight-`w` centre, candidates with exactly `k` replaced
support positions form an orbit of size

```text
2^w choose(w,k) choose(n-w,k).
```

The checker sums only the three orbit counts `k=0,1,2`, but that represents
all signs in a cloud of size

```text
192820470611922262353361155038354580963328,
log2(M) = 137.1463101125.
```

Radius one has only `117.5807` bits of capacity. Radius two therefore removes
the trivial cardinality obstruction: if all of its components overlapped
perfectly, Rényi slack `r >= 13.995` could in principle reach 128 bits.

It does not overlap. For every nonidentity `t`, a coefficient of `t-s` is odd
when the supports differ. When the support is unchanged, `(t-s)/2` has an odd
coefficient. Conditioning on all other mask coefficients makes one coefficient
of

```text
(t-s)*a + g*(t^2-s^2)
```

uniform on the whole 640-bit torus or on one parity coset. For the continuous
Gaussian model, the elementary theta estimate

```text
sum_{z in Z} exp(-z^2/(2 sigma^2)) < 4 sigma
```

therefore bounds every nonidentity top-row kernel by `2^-604` in expectation.
After summing all orbit members,

```text
log2 E[K_raw - 1] < -466.8537,
Pr[K_raw > 1 + 2^-128] < 2^-338.8537.
```

The barycentre factor can only decrease the certified effective size. Even as
`r` tends to infinity, the good-mask result supplies less than
`4.24e-39` security bits, rather than 128. This rejects the tested
support-radius-two Gaussian/Gibbs cluster; it does not prove insecurity or
exclude a nonlocal or computational reduction.

The source implementation is checked separately rather than silently replaced
by a continuous density. Multi-limb Gaussian noise is materialized through a
total saturating signed-`int64` lift, with exact support
`[-2^63,2^63-1]`. A selected mean coordinate at circular distance at least
`2^64` consequently has disjoint shifted noise support. Uniform-coset counting
and a union bound over the same cloud give

```text
Pr[some tested neighbor is within one noise diameter] < 2^-436.8537.
```

Thus the implemented finite-support channel has zero centre-to-neighbor
overlap on the good mask event; it does not repair the continuous-Gaussian
cluster bound. The companion Lean module
`FormalProof4FHE.RLWE.TFHEppLvl5BootGaussianClusterScreen` checks all orbit,
Markov, union-bound, Rényi-order, and signed-support arithmetic. The integer
theta-sum inequality and the source-to-uniform-coordinate identification stay
visible as the two analytic/refinement boundaries.

## Exact two-level Smith and complete-codebook screen

`lvl5boot_two_smith.py` applies the exact two-level principal-image size proved
in `FormalProof4FHE.RLWE.RankOneHNFLossinessTwoSmith`:

```text
ord_pi(Delta) = e*n + v,
|Delta R| = 2^(K*n - (e*n+v)).
```

Run it against the source-bound TFHEpp checkout with:

```bash
python3 python/proof/lvl5boot_two_smith.py --self-test
python3 python/proof/lvl5boot_two_smith.py --json
```

For distinct fixed-weight ternary secrets, the coefficientwise two-adic
exponent is zero when supports differ and one when they agree. A nonzero
primitive binary difference has Hasse depth at most `n-1`, so every secret
difference has total order at most `2*n-1 = 65535`. At `K=640`, its principal
image therefore has at least `20905985` bits.

The actual modular-Gaussian implementation first creates a signed-int64
coefficient. Two supported errors differ in fewer than `2^65` possible
integer values. Hence the complete degree-32768 error-difference support has
cardinality strictly below `2^(65*n)`, without assuming coefficient or row
independence. For an actual TFHEpp row, the public mask `A` is uniform and

```text
(t-s) * (A + g*(t+s))
```

is uniform on the principal image generated by `t-s`. Exact image counting
plus a union bound gives:

```text
radius-two cloud:       Pr[some shifted support overlaps] < 2^-18775927,
complete secret space: Pr[some shifted support overlaps] < 2^-18775027.
```

Thus changing the Hasse information set, evaluating a correlated Fourier
sum, or choosing a nonlocal subset of the actual fixed-weight codebook cannot
repair this direct statistical-cluster route for the uniform-mask row. This
is a proof-method obstruction, not an attack: ordinary RLWE is intended to be
computationally, rather than information-theoretically, hard.

The checker deliberately does not transfer that conclusion to the proof-only
lossy branch. TFHEpp has no concrete sampler for its DSPR/NTRU denominator,
numerators, and Hermite masks. A future instantiation must emit the
descriptor/leakage-weighted total-order histogram of

```text
(t-s) * (z_j + f*g_j*(t+s)).
```

With the coarse full-codebook support bound, the first total order not already
rejected at a 128-bit target is

```text
18840435 = 574*n + 31603.
```

This threshold is a diagnostic, not evidence that such descriptors occur.
The JSON report records it explicitly so a future sampler can be plugged in
without changing the criterion. It also contains exact binary Hasse dual-code
weight enumerators (depths 1, 2, 4, 8, and 16 by default) for the IID formula

```text
p_e^n * 2^-v * sum_w A_w * beta_e^w.
```

No numerical IID mass is assigned to `std::normal_distribution`: the C++
standard does not specify a portable exact PMF for that implementation path.

## Exact Hasse/descriptor certificate evaluator

`twosmith_exact_evaluator.py` implements the exact certificate interface from
`sketch/twosmith_exact_parameter_note.tex`. It uses `Fraction` throughout and
rejects JSON floating-point literals; probabilities must be integers or strings
such as `"17/32"`.

Run all deterministic cross-checks with:

```bash
python3 python/proof/twosmith_exact_evaluator.py --self-test
python3 python/proof/test_twosmith_exact_evaluator.py
```

Every regular command consumes one JSON object and emits exact JSON:

```bash
python3 python/proof/twosmith_exact_evaluator.py error-table pmf.json
python3 python/proof/twosmith_exact_evaluator.py local local.json
python3 python/proof/twosmith_exact_evaluator.py secret-hist secrets.json
python3 python/proof/twosmith_exact_evaluator.py ord polynomial.json
python3 python/proof/twosmith_exact_evaluator.py descriptor descriptor.json
python3 python/proof/twosmith_exact_evaluator.py aggregate histogram.json
```

The input roots are, respectively:

```json
{"K": 3, "pmf": {"-1": "1/4", "0": "1/2", "1": "1/4"}}
{"n": 8, "e": 0, "v": 3, "p": "1", "p_next": "1/2"}
{"n": 8, "w": 2, "valuations": [0, 1, 2, 3, 4, 5, 6, 7]}
{"K": 4, "polynomial": [1, -1, 0, 0]}
{"K": 4, "n": 4, "s": [1, 0, -1, 0], "t": [0, 1, 0, -1],
 "z": 2, "f": 1, "g": [1, 1, 0, 0]}
{"rows": [{"factors": {"0,0": "2"}}],
 "histogram": [{"tuple": [[0, 0]], "weight": "1"}]}
```

The local command emits the exact Lucas recursion trace and both the exact
Hasse factor and proved product upper bound. The secret command uses the dual
kernel recursion and emits the two signed fixed-weight branches. The descriptor
command performs negacyclic arithmetic modulo `2^K`, checks capped valuation
addition for the complete multiplier, and records `(e,v)`. The aggregation
command requires the full joint row tuple; it never manufactures a product of
marginal histograms.
