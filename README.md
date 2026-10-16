# Parameter-Selection
To determine TFHE's parameter, run lwe-estimator.

## Running the Python scripts

Most legacy scripts under `python/` require SageMath and the lattice-estimator.
The newer BFV, CLPX, and GL noise estimators also run with ordinary Python.
An Apptainer (Singularity) container definition is provided for a reproducible
Sage environment.

The scalar TFHE → BFV multiplication → TFHE circuit has a standalone
conditional whole-run error analysis:

```sh
python3 python/BFVswitch.py --products 256 --failure-bits 40
python3 python/BFVswitch.py --operand-bits 9 --products 256 --failure-bits 40
python3 python/BFVswitch.py --operand-bits 16 --tune-16bit --products 256 --failure-bits 40
python3 -m unittest discover -s python -p test_bfv_switch.py -v
sage -python python/estimates/bfv_switch_security.py --output /tmp/bfv-switch-security.json
```

The error report explicitly assumes primitive sub-Gaussian tails and FFT
rounding envelopes; it is not an unconditional proof. See the companion
`TFHEpp/docs/BFVSchemeSwitching.md` for the circuit, measurements, and assumptions.
`--operand-bits` models wider scalar profiles with a full double-width product.
The original fixed profile meets the conditional target through 9-bit operands,
not at 10 bits or above. Width 16 now selects a separate N=8192 profile with a
paired two-bit forward conversion. The conditional whole-256-product bound
is approximately 2^-83.9. `--bitwise-forward` models the previous tuned circuit
(approximately 2^-54.3); `--tune-16bit` reproduces its original bounded,
bitwise-forward parameter search. Neither option changes compiled C++ code.
`--legacy-profile` reproduces the old
untuned width-16 failure. An experimental correctness pass is not a failure bound.

### Narrow CKKS parameter tuning

The ePrint 2026/2119 narrow CKKS DD port is tuned here. Its average-case
refresh screen uses the existing GL estimator's variance arithmetic:

```sh
python3 python/narrowCKKSnoise.py --n 512
python3 python/narrowCKKSnoise.py --search --max-key-modulus 218
python3 -m unittest discover -s python -p test_narrow_ckks_noise.py -v
sage -python python/estimates/narrow_ckks_security.py --output /tmp/narrow-security.json
```

The retuned N16384 CPU profile is qualified in
[`narrow_ckks_ll14_retuned_dense_security.json`](python/estimates/narrow_ckks_ll14_retuned_dense_security.json).
It uses continuous sigma 32768 for q440 selectors, sigma 4096 for q433 DD
and q435 native automorphisms, and native primary digits of 32 bits.
The weakest dense classical estimates are 129.30/130.08/130.65 bits under
MATZOV (selectors/native/DD), and 131.62/132.37/133.07 under BDGL16.
The MATZOV sparse proxy is 143.81 bits. Exact rotation event counts across
several guessing dimensions and layouts, with an optimistic union-success
bound and dictionary size assumption, reached 129.09 bits in the tested grids.
These are model-dependent screens; the exact error PMF, scalar GSA,
structured-RLWE and circular/KDM assumptions remain explicit in the report.
No result is a security certificate or a global attack lower bound.

All three real functional evaluations (GeLU, Tanh and Sigmoid) and their
following multiplications pass the original accuracy gates, confirmed by
independent integer decoding. Sampled residual checks verify the actual
saved key noise. The TFHEpp `ckks_narrow_ll14_retuned_*` setup/replay targets
use distinct archive profiles, reject historical weak keys, and reuse saved
keys and input between evaluations. Historical weak public key material is
not upgraded by this profile.

```sh
sage -python python/estimates/narrow_ckks_security.py --case dd_ll14_native_automorphism_keys --case dd_ll14_direct_dd_master --noise-stddev 4095.6010789816987 --cost-model MATZOV --max-beta 4096 --output /tmp/retuned-dense.json
sage -python python/estimates/narrow_ckks_security.py --case dd_ll14_basis_direct_selectors --noise-stddev 32767.601060377354 --cost-model MATZOV --max-beta 4096 --output /tmp/retuned-selectors.json
python3 python/estimates/narrow_ckks_block_ring_optimistic.py --self-test
sage -python python/estimates/narrow_ckks_block_ring_optimistic.py --estimator python/lattice-estimator --output /tmp/retuned-structured.json
```

`narrow_ckks_block_ring_optimistic.py` computes exact event masses and signed
supports for every ring rotation, grouping equivalent block-count geometries.
Its exhaustive small-ring test checks the counts against signed secret
projections. It compares rotation-labelled support with the largest individual
rotation support, an optimistic lower bound on dictionary union size. The
latter is a sensitivity calculation and does not demonstrate a feasible attack.
Full commands, grids, hashes and source revisions are retained in the report.

The historical low-noise retained-cache CPU profile is assessed in
[`narrow_ckks_current_cpu_security.json`](python/estimates/narrow_ckks_current_cpu_security.json).
At N=16384, the weakest classical dense-family estimates are:

| Family | log2(q) | BDGL16 bits | MATZOV bits |
|---|---:|---:|---:|
| Native automorphism keys | 435 | 128.96 | 126.76 |
| DD master | 433 | 129.56 | 127.32 |
| Direct selectors | 440 | 128.97 | 126.79 |

This profile fails the 128-bit screen under MATZOV, the pinned estimator's
default model. The script preserves its historical BDGL16 default; select
`--cost-model MATZOV` explicitly. TFHEpp truncates continuous normal samples
toward zero: configured sigma3.2 gives integer RMS2.828301687834186, and
sigma102.4 gives RMS102.00190833992261. The estimator uses variance-matched
discrete Gaussian proxies, not the exact error distribution.

```sh
sage -python python/estimates/narrow_ckks_security.py --case dd_ll14_native_automorphism_keys --case dd_ll14_direct_dd_master --noise-stddev 2.828301687834186 --cost-model MATZOV --max-beta 4096 --output /tmp/current-dense.json
sage -python python/estimates/narrow_ckks_security.py --case dd_ll14_basis_direct_selectors --noise-stddev 102.00190833992261 --cost-model MATZOV --max-beta 4096 --output /tmp/current-selectors.json
```

The report embeds all generic attack results, bounded structured sensitivities,
commands, pinned revision, source/evidence hashes and regression checks.
The q56 uniform-sparse proxy gives143.89 bits under BDGL16; bounded block/rotation
sensitivities reach137.05 bits under MATZOV. Neither is a security certificate.
The latter reuses geometry-only posterior histograms and recomputes decoding
with `narrow_ckks_block_ring_joint_decode.py --noise-stddev`; omitting that option
preserves the saved profile's noise. Unsearched geometry, decoder correlations,
ring/KDM assumptions and the exact error distribution remain unqualified.
Retune with margin and requalify numerical accuracy and performance before
claiming128-bit security. This assessment generated no keys or inputs and
changed no CKKS parameters.

The noise search uses an average-case BSGS model for S2C and sparse
encapsulation; FFT error and security remain excluded. It follows
TFHEpp's default tensor-rescale-before-relinearization order and
unsuppressed gadget-key noise. Pass `--relinearize-before-rescale` to screen
the explicit DD API that switches the unscaled tensor at the input modulus
then rounds its two-component ciphertext. Both orders exclude FFT error
and have separate DD key levels; it must not inherit SEAL's hybrid-RNS
noise behavior. An empty candidate list means that this bounded search found
no candidate, not an impossibility result. The eight-RMS-plus-sine criterion
is an engineering screen, not a Gaussian tail or whole-run failure proof.
`--auxiliary-switch-bits` screens scale-raised DD automorphisms and a
raised S2C/encapsulation pipeline. It preserves fresh input
noise, includes final division rounding, and reports the largest evaluation
key modulus and structured switch-key modulus. The narrow refresh pipeline
supports this variant through the schedule's optional auxiliary scale. The
screen still excludes FFT error and structured-secret security. Search bounds can be set with
`--min-scale-bits`, `--max-scale-bits`, `--min-plain-bits`, and
`--max-plain-bits`; defaults preserve the earlier search ranges.
The factor screen includes the longest tail-alignment binary rotation chain,
with one division per chain when the auxiliary scale is enabled. Earlier
saved reports predate this term and remain historical measurements/screens.
`python3 python/estimates/narrow_ckks_throughput_candidates.py` screens two
prospective HT schedules with nine retained late-DD multiplications. Their
operands are fresh encryptions of one, preserving the signal amplitude;
the screen includes scale reduction, input/key error, final rounding, and
complex-slot capacity. N=512 encrypted diagnostics in TFHEpp cover nine
retained levels; validation at each configured full ring size is outstanding.
The screen reports selector residency, raw and compact DD-chain storage, and
raw and compact S2C coefficient storage separately. TFHEpp's optional
file-backed selector source loads one group at a time. Its compact DD format
stores balanced secondary digits at their signed integer width; compact S2C
stores checked int64 coefficients and expands one BSGS group for evaluation.
Both preserve the original arithmetic exactly. For HT15, DD-chain storage
drops from about 630 GiB to 5.63 GiB, and S2C coefficient storage from 112 GiB
to 8 GiB. Full-size encrypted validation remains outstanding. Evaluation
keys, derived FFT/NTT caches, and temporary buffers are additional memory;
the selector-only peak is not the whole process's memory budget.
Derived cache payloads are also reported explicitly. The optional bounded
tail-rotation evaluator reuses one binary FFT key across matching batch
columns, applies bits in the original order, and divides the auxiliary lift
once at the end. The C++ storage audit and Python model agree: LL14's six
tail keys occupy 13.5 GiB together or 2.25 GiB one at a time; HT15's nine
occupy about 195.65 GiB together or 21.74 GiB one at a time. Additional
staging is needed for full HT15; full-ring accuracy requires encrypted
validation independently of these storage estimates.
Selector NTT caches can also be generated in chunks. The evaluator retains
modular sums across chunks, then performs inverse NTT, CRT and rounding once
for the complete group. `narrowCKKSnoise.py --selector-chunk-size 40
--batch-columns 2` reports both transient NTT storage and those modular sums.
For HT15, the 40-selector NTT payload is about 1.64 GiB versus 43.4 GiB for
a whole group. Raw selector residency and other evaluation workspaces still
count toward the process budget.
The staged FFT evaluator keeps factor keys during factor extraction, releases
them before preparing the product-tree keys, and releases those before output
conjugation. Its caches belong to the evaluation call; the public setup stays
unchanged. For full HT15, the largest derived-key phase is about 45.01 GiB.
The prospective whole-process envelope is about 109.62 GiB with a 10 GiB
reserve. Small-ring encrypted tests pass; full-size validation is outstanding.

`noiseestimation.narrow_frequency_basis.minimum_positive_frequency_basis`
plans the optional Section 7 construction for degrees 1 through 40. It finds
a smallest positive basis by exhaustive search and returns each frequency's
singleton or two-summand reconstruction. Degree 12 uses four blind rotations;
degree 40 uses nine. The saved `python/estimates/narrow_ckks_frequency_basis.json`
reports the plans. This construction needs one extra multiplication level and
applies complex coefficients after reconstructing the unweighted roots.
TFHEpp implements public ciphertext reconstruction with compact DD keys,
fixed-point weights, and a conjugation key at the resulting output level.
The output scale increases by the coefficient bit count. Gaussian h=31
small-ring series and retained-square tests pass the unchanged 23-bit gate;
full-ring basis validation and parameter qualification remain outstanding.
The N=512, h=31 Gaussian basis path also passes sigmoid, GeLU and windowed
tanh on [-4,4], including complex endpoints and retained DD squares, with
the unchanged 23-bit evaluation/target gates. Its report is in TFHEpp at
`docs/narrow-dd-small-gaussian-basis-fbs.json`.
Run `python3 python/test_narrow_frequency_basis.py` to check coverage,
small-degree minima against enumeration, and complex series reconstruction.

`frequency_basis_schedule_budget` accounts for the extra product level,
increased logical coefficient scale, fixed-point quantization at exact roots,
and retained-square amplitude. It separately checks whether an output or its
square can normalize to the existing refresh input. Saved prospective budgets
are in `python/estimates/narrow_ckks_frequency_basis_budgets.json`; they exclude
noise, FFT accuracy and security. The full-ring 505-bit selector candidate
supports a guarded square but cannot refresh that square. The 536-bit candidate
has the required normalization headroom and needs 640-bit native storage when
including its optional packed key. Neither full-ring candidate is validated.

TFHEpp's `ckks_narrow_benchmark --basis-storage [report.json]` reports the
compiled full-ring 536-bit candidate's payloads. `--small-basis-fbs-file
sigmoid|gelu|tanh|all new-cache-directory [report.json]` checks that candidate at
N=512; `--basis-fbs-file` runs N=32768. These use staged FFT caches, file-backed
selectors and all four or nine positive basis frequencies in each component's
tile. The full candidate's conservative process estimate is 68.58 GiB with a
10 GiB reserve; the proposed virtual-memory cap is 80 GiB. Its selector files
occupy 199.69 GB, so use the workspace's NFS filesystem for full-size runs.
The memory estimate is not measured full-size RSS. All three small-ring checks
pass the 23-bit gates, retained squares and dense-input normalization; another
nonlinear refresh also requires the squared values to fit its message domain.
Four pinned generic dense attacks pass the 128-bit screen (weakest 199.79 bits
at the optional packed modulus 590); sparse-switch proxy estimates are
nonfinite. Structured-secret and circular/KDM security remain unvalidated.

The `all` option generates one public setup and reuses one encrypted input for
all three functions. The separate `--small-lean-basis-fbs-file all` profile
uses 576-bit storage, primary digits of 18 bits and secondary digits of 8 bits;
its selector/root/coefficient scales are 518/47/31 bits. All three N=512
checks pass the same 23-bit gates, squares and dense normalization.
`--lean-basis-storage` reports the full compiled layout; the conservative
process envelope is 41.63 GiB with a proposed 48 GiB virtual-memory cap and
179.72 GB of selector files. The full N=32768 `--lean-basis-fbs-file all` run passed the fixed 23-bit
accuracy gates for all three functions, retained squares, and dense normalization;
it exited successfully and cleaned up its selectors. Security remains unqualified.
Independent layout and level checks are in `test_narrow_frequency_basis.py`;
budget, noise and dense-screen results are saved in
`python/estimates/narrow_ckks_lean_basis_{budget,noise,security}.json`.
The unchanged sparse-switch modulus of 92 bits retains the inconclusive
uniform-sparse proxy result. Generic dense screens at moduli 572 and 523 bits
do not qualify the structured secret or ring and circular/KDM assumptions.

`python/estimates/narrow_ckks_throughput_candidates.py --lean-dd` keeps the
HT level schedules unchanged while screening primary18/secondary8 digits.
Both N=512 HT diagnostics pass all nine retained products after disabling
GNU IPA-ICF for TFHEpp's benchmark target; full-ring checks remain pending.
The independent compiled-layout envelopes are 17.43 GiB for HT14 and
61.06 GiB for HT15, including a 10 GiB runtime reserve.

`python/estimates/narrow_ckks_ll14_fbs_plaintext_candidate.py` screens an LL14
FBS candidate with selector385/rootDelta48/plain37, outputQ108 and optional
packedQ436. It preserves the 437-bit modulus budget and a guarded square.
The unweighted factor/tree RMS excludes trigonometric coefficient
amplification and FFT error. Encrypted validation is tracked separately;
the full-ring GeLU baseline's 22.63-bit result remains a failed 23-bit control.

The plain37 script accepts repeated `--validation-report REPORT` arguments
and `--output PATH`. It independently checks the full N=16384 schedule,
fixed 23-bit target/evaluation gates, complex uniform approximation bound,
and propagated square gates. Missing functions keep aggregate encrypted
validation false; a failed report also returns exit status 1. Use a separate
output path when checking a repeat so the original failure evidence remains
available. Full GeLU and the untraced sigmoid repeat pass these checks; the
original plain37 sigmoid report fails and its cause remains unresolved.
Full tanh now passes the fixed 23-bit refresh/evaluation, complex uniform-bound,
and propagated-square checks. The complete repeat aggregate is recorded in
`python/estimates/narrow_ckks_ll14_plain37_repeat_validation.json`; the
original-report aggregate in `narrow_ckks_ll14_plain37_original_validation.json`
retains the sigmoid failure. Both contain all three functions with none
missing. The repeat aggregate is an accuracy result and does not establish
the original failure's cause. These checks do not qualify security.

`python/estimates/narrow_ckks_ll13_fbs_candidate.py` screens the low precision
FBS schedule selector217/rootDelta27/plain24/outputQ58 on [-4,4]. It uses
expanded selectors: optional packedQ254 exceeds the 218-bit LL13 budget.
It normalizes three bits before sigmoid/tanh squares and six before GeLU,
leaving productQ31 with scales24/21 and three guard bits. Fixed paper
precision gates are9.04/7.89/7.17 bits, respectively. The earlier key218
candidate passed small accuracy tests but its generic dense screen fell
below128 bits; those controls are retained. Dense217 passes the generic
screen at128.30 bits. Sparse-switch49 has an incomplete uniform-sparse
proxy result and does not qualify the public block constraint. Full encrypted
sigmoid validation passes at N=8192 with9.41 bits versus9.04 required,
including its retained square. Full GeLU also passes at8.18 bits versus7.17
required. Full tanh passes at 8.45 bits versus 7.89 required, including its
retained square. All three full-ring accuracy checks pass; structured-secret
security remains incomplete.

The candidate script also accepts repeated `--validation-report REPORT`
arguments and an optional `--output PATH`. It checks full N=8192 profile
metadata and the fixed target/evaluation and propagated square gates.
Partial results list missing functions; `encrypted_validation_passed` becomes
true only when sigmoid, GeLU, and tanh all pass. This flag covers accuracy;
`security_validated` remains false. It rejects small-ring or altered-gate
reports instead of promoting them to full-profile evidence.

To reproduce the LL13 process envelope, generate TFHEpp's `--storage REPORT`
and pass it to the candidate script with `--storage REPORT --memory-report
OUTPUT`. The script rejects a different schedule or a mismatch between the
compiled factor/DD/group payloads and the independent model. The current
estimate is12.91 GiB including10 GiB reserve, with a proposed16 GiB limit;
it is not measured full-ring RSS.

The security script defaults to
this repository's pinned lattice-estimator and keeps incomplete sparse-secret
proxy results separate from the dense-key screens. Full encrypted measurements
remain in TFHEpp's opt-in `ckks_narrow_benchmark` executable.

The streamed LL14 candidate search is reproducible with:

```sh
python3 python/narrowCKKSnoise.py --n 16384 --source-bits 48 --input-scale-bits 41 \
  --s2c-bits 37 --search --target-precision 16 --max-key-modulus 437 --streamed
```

The first candidate has selector modulus 379 bits, tree scale 49 bits,
plaintext scale 34 bits, and a 427-bit packed-key modulus. Its full NTT
cache would occupy about 223 GiB; preparing one public block at a time
needs about 6.2 GiB of cache. This estimate excludes ciphertext/key storage
and temporary arithmetic scratch. Smaller-ring Gaussian validation passed
in TFHEpp; security-sized encrypted validation remains required.

Functional outputs need an amplitude budget for the retained multiplication
level. A fixed-scale DD square leaves centered coefficient capacity
`2^(Qout - 2*scale_bits - 1)`. For real and imaginary slot components bounded
by `R`, a conservative square coefficient bound is `2*R^2`. This check is
separate from the standard sine refresh noise model:

```sh
python3 python/narrowCKKSnoise.py --n 16384 --source-bits 48 --input-scale-bits 44 \
  --s2c-bits 37 --selector-bits 379 --scale-bits 49 --plain-bits 34 \
  --message-bound 4 --range-only --retained-component-bound 4.01
```

The existing 100-bit output modulus fails this range check. Three guard bits
select output modulus 108, selector modulus 387, and packed-key modulus 435
at N=16384. These values fit the 437-bit modulus budget but require new
noise/security validation. `--retained-component-bound` also constrains the
ordinary parameter search and adds range evidence to single-profile reports.
Include approximation and encrypted evaluation error in the component bound.
The `dd_ll14_fbs_dense_selectors` and `dd_ll14_fbs_dense_packed_setup` security
cases screen the revised exposed key moduli; structured-secret security is
still a separate assessment.

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
apptainer exec --no-home --pwd /work --env DOT_SAGE=/tmp/.sage --bind "$(pwd):/work" sagemath.sif sage -python /work/estimates/shallowboot_structured_std128.py
apptainer exec --no-home --pwd /work --env DOT_SAGE=/tmp/.sage --bind "$(pwd):/work" sagemath.sif sage -python /work/estimates/shallowboot_lowdepth.py
python3 python/estimates/shallowboot_noise.py
```

The shallow-bootstrap screen covers TFHEpp's executable `n=1024`, `q=512`,
one-hot `h=64` source and its current `N=1024`, 32-bit-torus accumulator
proxy.  It reports the source's uniform-sparse-binary proxy separately from
the accumulator proxy: the one-hot structured-secret assumption still needs
the additional concrete analysis discussed in ePrint 2026/1730.

`shallowboot_lowdepth.py` records Algorithm 3's `n=1450`, `h=29`, `N=4096`,
`105 -> 50` bit schedule and screens its sparse-LWE source. Its Binary-NTT
RLWE security remains the paper's new assumption and is intentionally not
treated as a standard lattice-estimator result.
The checked local estimates and the source-security discrepancy are recorded
in [`python/estimates/shallowboot_lowdepth_results.md`](python/estimates/shallowboot_lowdepth_results.md).

`shallowboot_noise.py` follows the implemented Algorithm-3 PBC/QH-RLWE phase
variance, checks modulus wrap headroom, and applies every modulus switch. Use
`--lwe-h 29 --search-middle 13:49` to search the paper-sized chain's omitted
intermediate modulus. The recorded noise-feasible candidate uses
`--ring-n 8192 --modulus-chain 151 69 52 36 --switch-after 3 4 5`.
This estimator requires ordinary Python, not Sage.
The older boundary diagnostic printed by `shallowboot_lowdepth.py` excludes
ciphertext-product noise and must not be interpreted as a correctness margin.
Recorded output and the resulting diagnosis are in
[`python/estimates/shallowboot_noise_results.md`](python/estimates/shallowboot_noise_results.md).
The recorded `N=8192,Q=2^151,sigma'=0.75` ring proxy is 188.7 classical bits.
The general input `n=1450,h=37,q=512` proxy is 133.3 bits. The faster
structured input profile uses `n=1024,h=64` with 16-position one-hot blocks;
its uniform-sparse proxy is 170.3 bits and its actual distribution follows the
paper's structured-entropy assumption. The refreshed output key uses
`n=1450,h=60,Q'=2^15,sigma=3.2` (the smaller `n=1150` proxy is already 138.6
bits). The structured executable optionally stores every PBC `a` polynomial
as a 256-bit BLAKE3 seed, reducing its estimated PBC key from 0.797 to 0.398
GiB without increasing the measured online runtime.

**BatchBoot reassessment** (run from `python/`):

```bash
# Conservative functional-BatchBoot noise screen; no Sage dependency.
python3 BatchBoot.py --profile tfhepp-lvl1 --slots 1024 --hamming-weight 34

# Estimate both the sparse source and the accumulator target secrets.
apptainer exec --no-home --pwd /work --env DOT_SAGE=/tmp/.sage \
  --bind "$(pwd):/work" sagemath.sif \
  sage -python /work/BatchBoot.py --profile tfhepp-lvl1 \
  --slots 1024 --hamming-weight 34 --security
```

`BatchBoot.py` follows the concrete radix-4 selector count in
`TFHEpp/include/tfhe/batchboot.hpp`: each EMP uses four external products per
stored radix-4 digit, including the final partially used digit when
`log2(slots)` is odd.  It treats the sparse source key and accumulator target
key as different LWE instances, since evaluation keys are encrypted under the
target key.  The reported failure value is a conservative Gaussian screening
bound; it does not cover FFT numerical error or serve as a correctness proof.
The built-in profiles are TFHEpp structs, not proposed production settings:
for a separate sparse input parameter, set `--source-alpha-bits` and rerun
`--security` for both source and target instances.
The latest local assessment and its exact estimator assumptions are recorded in
[`python/estimates/batchboot_results.md`](python/estimates/batchboot_results.md).
The same-checkout throughput comparison with TFHEpp's existing block-binary
implementation is recorded in
[`python/estimates/batchboot_vs_blockbinary.md`](python/estimates/batchboot_vs_blockbinary.md).

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

**TFHEpp fixed-weight quadratic-KDM proof screen** uses exact rational
arithmetic and ordinary Python:

```bash
python3 python/proof/lvl5boot_renyi.py --self-test
python3 python/proof/lvl5boot_renyi.py --json
python3 python/proof/lvl5boot_gaussian_cluster.py --self-test
python3 python/proof/lvl5boot_gaussian_cluster.py --json
```

This checks the first necessary mean-energy margin of the equal-covariance
Rényi certificate against the actual `lvl5bootparam` secret, noise, and 35
pre-Double-Decomposition gadget rows. The current parameter fails that
specific sufficient condition already at its top gadget row. The matching
Lean module checks the exact obstruction arithmetic; see
`python/proof/README.md` for the result and its encoding/posterior scope.

**TFHEpp subset-key delayed-projection screen** source-binds the current word
sizes, dimensions, and error widths, rejects the lifted invertible-minor solver
by its exact noise budget, and evaluates the canonical-ternary second-moment
failure bound for genuinely short high-modulus preimages:

```bash
python3 python/proof/tfhe_subset_joint_screen.py --self-test
python3 python/proof/tfhe_short_preimage_screen.py --self-test
python3 python/proof/tfhe_short_preimage_screen.py --json
```

The screen now proves information-theoretic random-matrix existence for the
complete target family and checks the disjoint-block Gram/noise budget. It is
not a complete security certificate: efficient public recovery has an SIS
barrier, and the exact finite C++ error-law comparison remains open.

The Gaussian-cluster follow-up source-binds the actual 640-bit TFHEpp
parameters and evaluates a radius-two fixed-weight orbit cloud without
enumeration. That cloud has `137.1463` bits of perfect-overlap capacity, but
the uniform mask makes its expected nonidentity overlap smaller than
`2^-466.85` in total. Markov gives `K <= 1+2^-128` except with probability
below `2^-338`, so optimizing the Rényi order cannot produce a meaningful
lossiness bound. For the implemented signed-`int64` noise channel, the same
neighbors are support-disjoint except with probability below `2^-436`. The
matching Lean module checks the finite arithmetic and support-separation
theorem; details and precise scope are in `python/proof/README.md`.

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
