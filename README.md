# Parameter-Selection
To determine TFHE's parameter, run lwe-estimator.

## Running the Python scripts

The scripts under `python/` require SageMath and the lattice-estimator. An Apptainer (Singularity) container definition is provided for a reproducible environment.

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
```

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

**CLPX scheme-switch noise estimation** can be run with normal Python when
SciPy is installed:

```bash
python3 python/CLPXnoise.py
python3 python/CLPXnoise.py --direction tlwes-to-clpx --paper-ss2clpx --validbit 8
python3 python/CLPXnoise.py --direction clpx-to-tlwes --validbit 8 --numdigit 4 --basebit 2
python3 python/CLPXnoise.py --direction switched-multiplication --paper-ss2clpx --validbit 8 --max-mults 8 --mult-chain square
python3 python/CLPXnoise.py --direction all --validbit 8 --max-mults 16
```

`CLPXnoise.py` follows the TFHEpp operation sequence in
`../TFHEpp/include/bfv-clpx.hpp` for `TLWES2CLPXIKS` and
`CLPX2TLWESIKSanybit`.  It composes the existing TFHE bootstrapping,
identity-key-switching, and annihilate-packing formulas from
`python/noiseestimation/keyvariation.py`.  The default CLPX preset in
`python/noiseestimation/params/clpx.py` mirrors the local TFHEpp default
`include/params/128bit.hpp` CLPX test path.  Programmable bootstrapping is
modeled as refreshing the output encryption noise; the script also reports the
largest internal PBS-input variance and a bin margin, because CLPX digit
extraction can fail semantically even when the final refreshed TLWE noise is
small.

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
