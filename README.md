# Parameter-Selection
To determine TFHE's parameter, run lwe-estimator.

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

## References

The noise estimator (`python/noiseestimation/keyvariation.py`) is based on the following papers. PDFs are stored in the `references/` directory.

- Ilaria Chillotti, Damien Ligier, Jean-Baptiste Orfila, and Samuel Tap, "Improved Programmable Bootstrapping with Larger Precision and Efficient Arithmetic Circuits for TFHE," IACR ePrint 2021/729. https://eprint.iacr.org/2021/729
- Thomas de Ruijter, Jan-Pieter D'Anvers, and Ingrid Verbauwhede, "Don't be mean: Reducing Approximation Noise in TFHE through Mean Compensation," IACR ePrint 2025/809. https://eprint.iacr.org/2025/809
- Ruida Wang, Jincheol Ha, Xuan Shen, Xianhui Lu, Chunling Chen, Kunpeng Wang, and Jooyoung Lee, "Refined TFHE Leveled Homomorphic Evaluation and Its Application," IACR ePrint 2024/1318. https://eprint.iacr.org/2024/1318
- Mariya Georgieva Belorgey, Sergiu Carpov, Nicolas Gama, Sandra Guasch, and Dimitar Jetchev, "Revisiting Key Decomposition Techniques for FHE: Simpler, Faster and More Generic," IACR ePrint 2023/771. https://eprint.iacr.org/2023/771
