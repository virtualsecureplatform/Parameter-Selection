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
python3 python/BFVnoise.py --preset tfhepp-lvl3simd-boot --B 15 --qbits-range 128:256:32
python3 python/BFVnoise.py --preset tfhepp-lvl3simd-boot --B-range 1:15 --scalar-mode unsigned-average
python3 python/BFVvalidate.py
```

`BFVnoise.py` implements the invariant-noise variance formulas from `600.pdf`
("Improving and Automating BFV Parameters Selection: An Average-Case Approach").
The default TFHEpp bootstrap preset estimates the final digit-removal
`PolyEval` over `PrimePower2Param`, so its plaintext modulus is
`114689^2`.  The default digit-removal degree is `4*B+1`, matching
`GetLowestDigitRemovalPolynomialOverRange(p, B)`.

Current limitations: the estimator uses the paper's independent-ciphertext
average-case multiplication model.  Evaluating a high-degree polynomial in one
ciphertext reuses dependent powers, and TFHEpp's double-decomposition
relinearization is only approximated by the paper's key-switching variants.
Treat the output as screening data until the dependent-circuit model is added.

`BFVvalidate.py` reproduces the OpenFHE-based validation parameters from
`600.pdf` Table 7: `t=65537`, `sigma=3.19`, `chi_s=chi_u=U3`, Hybrid key
switching, HPSPOVERQ multiplication, and `log2(q) ~= 60` for encryption/addition
or `log2(q) ~= 120` for one multiplication.  It checks the paper's average-case
"our" column, not the experimental OpenFHE samples.

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
