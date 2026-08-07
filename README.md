# Parameter-Selection
To determine TFHE's parameter, run lwe-estimator.

## Running the Python scripts

Most legacy scripts under `python/` require SageMath and the lattice-estimator.
The newer BFV, CLPX, and GL noise estimators also run with ordinary Python.
An Apptainer (Singularity) container definition is provided for a reproducible
Sage environment.

### Prerequisites

- [Apptainer](https://apptainer.org/) (or Singularity) installed on your system
- Git submodules initialized:
  ```bash
  git submodule update --init --recursive
  ```

### Building the container

A pre-built `python/sagemath.sif` may already be present. To rebuild:

```bash
apptainer build python/sagemath.sif python/sagemath.def
```

### Running scripts

**Noise estimation scripts** (run from `python/`):

```bash
cd python
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/TFHEnoise.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/TFHEint.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/TFHElvl21.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/TFHElvl22.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/manyLUT.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/shortlwe.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/DirectPDF.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/CCbound.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/ConcreteCCbound.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/BFVnoise.py --preset tfhepp-lvl3simd-boot --B 15 --qbits-range 128:256:32
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/BFVvalidate.py
```

**Lattice security estimation scripts** (run from `python/`):

```bash
cd python
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/newTFHE.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/estimates/TFHE586.py
apptainer exec --bind "$(pwd):/work" sagemath.sif sage -python /work/estimates/verify128bit.py
apptainer exec --no-home --pwd /work --env DOT_SAGE=/tmp/.sage --bind "$(pwd):/work" sagemath.sif sage -python /work/estimates/TFHEprus_secure_128.py
```

The `--no-home --env DOT_SAGE=/tmp/.sage` form avoids Sage trying to create
`~/.sage` on read-only container paths in rootless Apptainer environments.

**Block-binary key estimates** from ePrint 2023/958:

```bash
cd python
apptainer exec --no-home --pwd /work --env DOT_SAGE=/tmp/.sage --bind "$(pwd):/work" sagemath.sif sage -python /work/estimates/block_binary.py paper-table --quiet
apptainer exec --no-home --pwd /work --env DOT_SAGE=/tmp/.sage --bind "$(pwd):/work" sagemath.sif sage -python /work/estimates/block_binary.py search --profile tfhe636 --lwe-only --target 128
```

`paper-table` prints the paper's Table 1 values next to the locally computed
wrapper estimates.  `--lwe-only` skips the TRLWE-as-LWE side and is the mode to
use when searching for the TLWE dimension at a fixed TLWE noise size.

**Parameter search** (run from the repository root):

```bash
apptainer exec --bind "$(pwd):/work" python/sagemath.sif sage -python /work/python/noiseestimation/search_lvl03param.py
```

**BFV average-case noise estimation** can also run with normal Python when
SciPy is installed:

```bash
python3 python/BFVnoise.py --preset tfhepp-lvl3simd-boot --B 15 --qbits-range 128:512:64 --error-std 3.19
python3 python/BFVnoise.py --preset tfhepp-lvl3simd-boot --B-range 1:15
python3 python/BFVnoise.py --preset tfhepp-lvl5-boot --B 15
python3 python/BFVvalidate.py
```

**GL-SHIP Double-Decomposition noise estimation** uses only the Python
standard library:

```bash
python3 python/GLnoise.py
python3 python/GLnoise.py --preset n512p17 --stages
python3 python/GLnoise.py --preset n1024p17 --model correlated
python3 python/GLnoise.py --optimize-tree-scale
python3 python/GLnoise.py --preset n512p17 --arithmetic legacy --json
python3 python/GLnoise.py --preset n512p17 --target-precision 12 --strict
python3 python/GLvalidate.py
```

`GLnoise.py` expands the GL paper's abstract `epsilon_HE` term for TFHEpp's
coefficient-domain implementation. It tracks fresh and grouped-StC phase
noise, dense-to-sparse switching, encrypted masked columns, X-only HMux,
the five-level balanced product tree, output conjugation, two Gaussian
channels, and the final Y reconstruction. Double Decomposition is modeled
with active primary digits, exact Bbar limbs, evaluation-key error, and
modulus-down rounding. The default fused-DD mode follows the same operation
boundaries as hybrid RNS: masked candidates remain under `P*Q` until the one
ModDown in Equation (15), each HMux stage accumulates its body/mask and radix
branches before one ModDown, and product-tree nodes relinearize at their input
modulus before rescaling the resulting two-component ciphertext. DD replaces
the decomposition and recomposition representation, not this algorithm.

`--arithmetic legacy` reproduces the former TFHEpp path, which performed a
ModDown for every candidate and HMux switch and rescaled all four product
tensor components before relinearization. Its rescale floor is
`r00 + (r01+r10)s + r11*s^2`; the fused path has the ordinary `r0+r1*s`
floor. `--masked-moddown` remains available as a single-boundary diagnostic.
The `independent` model is an average-case variance sum; `correlated` is a
worst-aligned sensitivity screen.

The paper's precision experiment samples each real and imaginary input slot
uniformly from `[-1,1]`. The estimator propagates the corresponding GL
coefficient variance through grouped StC instead of assigning magnitude one
independently to every polynomial coefficient. Full-bootstrap sine and
initial-encoding errors use the same distribution; the direct half-bootstrap
continues to report a deterministic bounded-message sine error.

The three presets copy `Q`, `P`, StC size, gap, sparse-secret, and window data
from ePrint 2026/811 and copy TFHEpp's DD base/storage widths. The paper does
not publish the complete per-prime schedule, so q0 and product-tree scales are
reconstructed from the matched ePrint 2025/784 SHIP profiles and remain CLI
overrides. The output includes torus-storage, 128-bit security-ceiling,
phase-wrap, outside-depth, and precision margins. `--optimize-tree-scale`
spends only the already available output-depth headroom and selects the first
uniform tree scale that reaches the requested precision. With fused DD, the
reconstructed `n512p17` and `n1024p17` profiles reach the paper's 14.94- and
15.88-bit measurements at 47- and 50-bit tree scales, respectively, without
increasing `Q` or `P`. The reconstructed `n256p17` schedule remains below
target even at its maximum feasible 26-bit tree scale and is deliberately
reported as unresolved rather than production-ready.

This is a conservative parameter screen, not an RLWE security estimate or a
correctness proof. The default precision target is the paper's measured
hybrid-RNS result; a `BELOW TARGET` result means that the selected DD schedule
does not reproduce that precision, not that DD intrinsically needs a larger
ring or that the paper's RNS measurement is wrong.

**CLPX scheme-switch noise estimation** can be run with normal Python when
SciPy is installed:

```bash
python3 python/CLPXnoise.py
python3 python/CLPXnoise.py --direction tlwes-to-clpx --paper-ss2clpx --validbit 8
python3 python/CLPXnoise.py --direction clpx-to-tlwes --validbit 16 --numdigit 9 --basebit 2
python3 python/CLPXnoise.py --compare-reverse-options --validbit 128
python3 python/CLPXnoise.py --direction switched-multiplication --paper-ss2clpx --validbit 8 --max-mults 8 --mult-chain square
python3 python/CLPXnoise.py --direction all --validbit 8 --max-mults 16
```

`CLPXnoise.py` follows the TFHEpp operation sequence in
`../TFHEpp/include/clpx/bfv-clpx.hpp` for `TLWES2CLPXIKS` and
`CLPX2TLWESIKSanybit`.  It composes the existing TFHE bootstrapping,
identity-key-switching, and annihilate-packing formulas from
`python/noiseestimation/keyvariation.py`.  The default CLPX preset in
`python/noiseestimation/params/clpx.py` mirrors the local TFHEpp default
`include/params/128bit.hpp` CLPX test path.  Programmable bootstrapping is
modeled as refreshing the output encryption noise; the script also reports the
largest internal PBS-input variance and a semantic interval margin, because
CLPX digit extraction can fail semantically even when the final refreshed TLWE
noise is small.

For TFHE-to-CLPX, the estimator also includes the paper's
`q^2 / 2^(2w)` variance contribution from the finite-width approximation of
`Delta_b`.  Failure probabilities are evaluated in the log domain so results
below floating-point `erfc` range remain finite.  The reported PBS input margin
uses the source TLWE message interval; it is not a requirement that input noise
fit within one modulus-switch rounding bin.

The paper model also keeps the two meanings of `w` separate: `Delta_b`'s
truncation width is a scheme-switch argument, whereas relinearization uses the
TFHEpp gadget base (`P.B`) and level count (`P.l`).  Margins for the reverse
path use their operation-specific decision intervals: `1/32` for rounded
radix-2 digits, half a decomposition bin for HomDecomp and final bit
extraction, and the corresponding carry-LUT interval.

`--compare-reverse-options` compares the implemented 16-bit-block path with
`basebit=4` candidates.  Packing the four final bit LUTs into one PBS is not a
valid option: negacyclic PBS requires `f(x+1/2)=-f(x)`, while the lower-bit
functions repeat after a half-torus shift.  Their distinct input shifts are
therefore necessary.  The default lvl2 PBS gadget (`l=4`, `Bgbit=9`) does not
leave a sufficient HomDecomp margin for `basebit=4`.  The reverse-only
candidate keeps four rows and selects
`Bgbit=10`, reducing lvl2 PBS variance by about 5.2 bits without increasing
the bootstrapping-key row count.  This raises the estimated `basebit=4`
HomDecomp margin from about 2.5 to 7.5 variance bits.  The same preset also
raises the current path's conservative aggregate semantic-failure estimate
from about `2^-56` to well beyond `2^-128`; it is therefore used for both
optimized candidates in the comparison.

The `--direction switched-multiplication` mode is an approximate depth screen
for `CLPXMult` (`TRLWEMultWithoutRelinerizationCLPX + Relinearization`) using
the estimated `TLWES2CLPXIKS` output noise as the initial CLPX noise, with the
same `--validbit`, `--num-multi`, `--shift`/`--shiftnum`, and `--w` arguments,
unless `--input-log2-var` is provided.  The older `--direction multiplication`
spelling is kept as an alias.  Use `--paper-ss2clpx` for the Nagai et al.
setting implemented by TFHEpp's `SS2CLPX.hpp`: CLPX base `b=2`, `Lutnum=4`,
`shiftnum=5` (TFHEpp template `shift=4`), and `w=20`.  In this mode the
post-switch multiplication estimate uses Equations (44)-(48) from the paper and
reports one supported multiplication for the 8-bit default before the next
CLPX-to-TFHE switch.  Without `--paper-ss2clpx`, the estimator keeps the older
TFHEpp default CLPX path (`plain_modulus=8`) and treats multiplication as a
bounded-digit screening model rather than the paper's direct post-switch path.

`BFVnoise.py` implements the invariant-noise variance formulas from `600.pdf`
("Improving and Automating BFV Parameters Selection: An Average-Case Approach").
The default TFHEpp bootstrap preset estimates the final digit-removal
`PolyEval` over `PrimePower2Param`, so its plaintext modulus is
`114689^2`.  By default, `BFVnoise.py` builds the same bounded
digit-removal polynomial as `GetLowestDigitRemovalPolynomialOverRange(p, B)`
and evaluates its actual degree and scalar coefficients.  Use
`--poly-source degree` to run a degree-only sweep.

For `PolyEval`, the default `--circuit-model dependent` applies the Section 7
dependent-ciphertext bounds from `600.pdf`.  As in the paper's identical-input
examples, this omits the unknown `Var((nu*nu')|i)` term.  TFHEpp's
double-decomposition relinearization is still approximated by the paper's
key-switching variants, so treat the output as screening data.

The TFHEpp presets use TFHE-style normalized fresh noise by default, so changing
`--qbits` alone keeps the normalized error fixed.  Use `--error-std` when you
want the BFV paper's fixed absolute error model for ciphertext-modulus sweeps.

`BFVvalidate.py` reproduces the OpenFHE-based validation parameters from
`600.pdf` Tables 7 and 10: `t=65537`, `sigma=3.19`, `chi_s=chi_u=U3`, Hybrid
key switching, HPSPOVERQ multiplication, and `log2(q) ~= 60` for
encryption/addition or `log2(q) ~= 120` for one multiplication.  It checks the
paper's average-case "our" column, not the experimental OpenFHE samples.

**Geometric-LWE-Estimator scripts** (run from `python/`; note the `cwd` must be set inside the submodule for sage's relative `load()` paths to resolve):

```bash
cd python
apptainer exec --bind "$(pwd):/work" sagemath.sif bash -c "cd /work/Geometric-LWE-Estimator/section_5_1 && sage /work/leakylwr.sage.py"
```

## Notation correspondence (TFHEpp ↔ Python ↔ papers)

This repo keeps two “views” of parameters:

- **TFHEpp**: `../TFHEpp/include/params/*.hpp` (preferred names)
- **Python noise estimator**: `python/noiseestimation/params/*.py` and `python/noiseestimation/keyvariation.py`

The table below summarizes the intended correspondence and meaning.

| Concept | TFHEpp name (C++) | Python name | Typical paper notation | Meaning / notes |
|---|---|---|---|---|
| TLWE/TRLWE polynomial degree | `n`, `nbit` | `n`, `nbit` | `N` | `n = 2^nbit` for ring variants |
| GLWE dimension | `k` | `k` | `k` | Number of polynomials in secret key (TRLWE has `k+1` components) |
| Torus modulus | implicit via `using T = ...` | `q` | `q` or `2^w` | Python explicitly sets `q = 2^w`; TFHEpp’s `q` is `2^{digits(T)}` |
| Fresh noise (stdev) | `α` | `α` | `α` or `σ` | TFHEpp `α` is normalized (torus); Python stores `α` in integer-torus units (`α = α_norm * q`) and often uses `σ = α^2` |
| Fresh noise (variance) | (derived) | `σ` | `σ^2` | Python convention: `σ = α^2` (variance in integer-torus units) |
| TRGSW main decomposition levels | `l` | `l` | `ℓ` | Number of gadget digits for the “body” part |
| TRGSW nonce decomposition levels | `lₐ` | `lₐ` | `ℓ` | Levels for the “mask/nonce” part (TFHEpp can use distinct params for each half) |
| TRGSW main base (bits) | `Bgbit` | `ℬbit` | `log2(B)` | `Bg = 2^{Bgbit}`; Python uses `ℬ = 2^{ℬbit}` |
| TRGSW nonce base (bits) | `Bgₐbit` | `ℬₐbit` | `log2(B)` | `Bgₐ = 2^{Bgₐbit}`; Python uses `ℬₐ = 2^{ℬₐbit}` |
| TRGSW main base value | `Bg` | `ℬ` | `B` | Power-of-two base |
| TRGSW nonce base value | `Bgₐ` | `ℬₐ` | `B` | Power-of-two base |
| **Double Decomposition** auxiliary levels | `l̅`, `l̅ₐ` | `l̅`, `l̅ₐ` | `ℓ̅` / “#limbs” | Enables DD external product / blind rotation in TFHEpp (e.g. `lvl3param` in `128bit.hpp`) |
| **Double Decomposition** auxiliary base (bits) | `B̅gbit`, `B̅gₐbit` | `B̅gbit`, `B̅gₐbit` | `K` (limb size) | Auxiliary base is `2^{B̅gbit}` (paper `K` bits) |
| Key switching digits | `t` | `t` | `t` or `ℓ_ks` | Number of decomposition digits in KS key |
| Key switching base (bits) | `basebit` | `basebit` | `log2(β_ks)` | KS base is `2^{basebit}` |
| Secret key distribution range | `key_value_min/max` | (via coefficients below) | (depends) | TFHEpp samples secrets uniformly in `[min,max]` |
| Secret key mean/variance | (derived) | `expectation_key_coefficient`, `variance_key_coefficient` | `μ_s`, `σ_s^2` | Used by the estimator when modeling key-dependent noise terms |
| BFV plaintext modulus | `plain_modulus` | `t` | `t` | For BFV bootstrap finalization, `PrimePower2Param` uses `t = p^2 = 114689^2` |
| BFV ciphertext modulus bits | `std::numeric_limits<T>::digits` | `q_bits` | `log2(q)` | TFHEpp `lvl3simdparam` uses 128-bit torus coefficients |
| BFV digit-error bound | `bfv_bootstrap_digit_error_bound` | `B` | `B` | Defines the bounded low digit removed by the final BFV bootstrap polynomial |

## References

The noise estimator (`python/noiseestimation/keyvariation.py`) is based on the following papers. PDFs are stored in the `references/` directory.

- Ilaria Chillotti, Damien Ligier, Jean-Baptiste Orfila, and Samuel Tap, "Improved Programmable Bootstrapping with Larger Precision and Efficient Arithmetic Circuits for TFHE," IACR ePrint 2021/729. https://eprint.iacr.org/2021/729
- Thomas de Ruijter, Jan-Pieter D'Anvers, and Ingrid Verbauwhede, "Don't be mean: Reducing Approximation Noise in TFHE through Mean Compensation," IACR ePrint 2025/809. https://eprint.iacr.org/2025/809
- Ruida Wang, Jincheol Ha, Xuan Shen, Xianhui Lu, Chunling Chen, Kunpeng Wang, and Jooyoung Lee, "Refined TFHE Leveled Homomorphic Evaluation and Its Application," IACR ePrint 2024/1318. https://eprint.iacr.org/2024/1318
- Mariya Georgieva Belorgey, Sergiu Carpov, Nicolas Gama, Sandra Guasch, and Dimitar Jetchev, "Revisiting Key Decomposition Techniques for FHE: Simpler, Faster and More Generic," IACR ePrint 2023/771. https://eprint.iacr.org/2023/771
- Craig Gentry and Yongwoo Lee, "Fully Homomorphic Encryption for Matrix Arithmetic," IACR ePrint 2025/1935. https://eprint.iacr.org/2025/1935
- Jung Hee Cheon, Guillaume Hanrot, Jongmin Kim, and Damien Stehlé, "SHIP: A Shallow and Highly Parallelizable CKKS Bootstrapping Algorithm," IACR ePrint 2025/784. https://eprint.iacr.org/2025/784
- Rostin Shokri and Nektarios Georgios Tsoutsos, "Low-Depth Bootstrapping for Matrix-Native FHE," IACR ePrint 2026/811. https://eprint.iacr.org/2026/811
