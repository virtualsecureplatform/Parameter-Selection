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
