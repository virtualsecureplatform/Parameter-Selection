#!/usr/bin/env python3
"""Source-bound arithmetic screen for the parity-aligned lvl02 candidate.

The candidate is intentionally a modified-format parameter profile.  It uses
the exact correlated BRK/aligned-KSK error law proved in FormalProof4FHE; it is
not executable through TFHEpp's native independent subset KSK.  Its suffix
source is no longer a suffix-subspace assumption: the exact parity reduction
identifies it with ordinary degree-1024 ternary RLWE using twice the BRK ring
sample count.

The screen binds the proposal to ``params/source_aligned_lvl02.hpp`` and
computes every finite dimension and correctness exponent with integer or
rational arithmetic.  Heuristic hardness estimates live in the companion
Sage script ``python/estimates/tfhe_lvl02_correlated_candidate.py``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Candidate:
    prefix_dimension: int
    ring_degree: int
    ring_rank: int
    gadget_levels: int
    gadget_basebit: int
    torus_bits: int
    fresh_sigma: int
    finite_error_coarse_count: int
    finite_error_coarse_scale_log2: int
    finite_error_dither_bits: int

    @property
    def suffix_dimension(self) -> int:
        return self.ring_degree - self.prefix_dimension

    @property
    def rows_per_control(self) -> int:
        return (self.ring_rank + 1) * self.gadget_levels

    @property
    def brk_rows(self) -> int:
        return self.prefix_dimension * self.rows_per_control

    @property
    def aligned_width(self) -> int:
        return self.brk_rows * self.ring_degree

    @property
    def suffix_rlwe_samples(self) -> int:
        return 2 * self.brk_rows

    @property
    def digit_radius(self) -> int:
        return 1 << (self.gadget_basebit - 1)

    @property
    def factor_energy_bound(self) -> int:
        return self.aligned_width * self.digit_radius**2

    @property
    def threshold(self) -> int:
        return 1 << (self.torus_bits - 4)

    @property
    def fresh_variance_bound(self) -> int:
        return self.fresh_sigma**2 * self.factor_energy_bound

    @property
    def chernoff_exponent(self) -> Fraction:
        return Fraction(self.threshold**2, 2 * self.fresh_variance_bound)

    @property
    def decomposition_bits(self) -> int:
        return self.gadget_levels * self.gadget_basebit

    @property
    def finite_error_variance(self) -> Fraction:
        coarse = self.finite_error_coarse_count * (
            1 << (2 * self.finite_error_coarse_scale_log2)
        )
        dither = Fraction(
            (1 << (2 * self.finite_error_dither_bits)) - 1, 6
        )
        return Fraction(coarse) + dither

    @property
    def finite_error_stddev(self) -> float:
        return math.sqrt(float(self.finite_error_variance))

    @property
    def finite_error_mgf_proxy(self) -> int:
        coarse = self.finite_error_coarse_count * (
            1 << (2 * self.finite_error_coarse_scale_log2)
        )
        dither = ((1 << (2 * self.finite_error_dither_bits)) - 1) // 3
        return coarse + dither

    @property
    def finite_error_chernoff_exponent(self) -> Fraction:
        return Fraction(
            self.threshold**2,
            2 * self.finite_error_mgf_proxy * self.factor_energy_bound,
        )


SELECTED = Candidate(
    prefix_dimension=1024,
    ring_degree=2048,
    ring_rank=1,
    gadget_levels=18,
    gadget_basebit=2,
    torus_bits=64,
    fresh_sigma=1 << 42,
    finite_error_coarse_count=255,
    finite_error_coarse_scale_log2=38,
    finite_error_dither_bits=32,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def _log2_two_sided_tail(exponent: Fraction) -> float:
    return 1.0 - float(exponent) / math.log(2.0)


def _fresh_sigma_ceiling(candidate: Candidate, failure_bits: int) -> int:
    denominator = (
        2.0
        * candidate.factor_energy_bound
        * (failure_bits + 1)
        * math.log(2.0)
    )
    return math.floor(candidate.threshold / math.sqrt(denominator))


def _require(source: str, pattern: str, description: str) -> None:
    if re.search(pattern, source, flags=re.MULTILINE) is None:
        raise ValueError(f"candidate header does not bind {description}")


def analyse(
    root: Path,
    failure_bits: int = 128,
    formal_proof_root: Path | None = None,
) -> dict[str, Any]:
    if failure_bits <= 0:
        raise ValueError("failure bits must be positive")

    header = root / "include/params/source_aligned_lvl02.hpp"
    key_header = root / "include/tfhe/key.hpp"
    parity_header = root / "include/tfhe/sourcealignedlvl02.hpp"
    trlwe_header = root / "include/tfhe/trlwe.hpp"
    if formal_proof_root is None:
        formal_proof_root = root.resolve().parent / "FormalProof4FHE"
    parity_proof = (
        formal_proof_root
        / "FormalProof4FHE/TFHE/TFHEppCandidateLvl02ParitySecurity.lean"
    )
    finite_error_proof = (
        formal_proof_root
        / "FormalProof4FHE/TFHE/SourceAlignedProofErrorSampler.lean"
    )
    parameter_proof = (
        formal_proof_root
        / "FormalProof4FHE/TFHE/TFHEppSourceAlignedParameterScreen.lean"
    )
    header_source = header.read_text(encoding="utf-8")
    key_source = key_header.read_text(encoding="utf-8")
    parity_source = parity_header.read_text(encoding="utf-8")
    trlwe_source = trlwe_header.read_text(encoding="utf-8")
    parity_proof_source = parity_proof.read_text(encoding="utf-8")
    finite_error_proof_source = finite_error_proof.read_text(encoding="utf-8")
    parameter_proof_source = parameter_proof.read_text(encoding="utf-8")
    candidate = SELECTED

    bindings = {
        "prefix_dimension": r"SourceAlignedLvl0Param[\s\S]*?\bn\s*=\s*1024\s*;",
        "ring_nbit": r"SourceAlignedLvl2Param[\s\S]*?\bnbit\s*=\s*11\s*;",
        "gadget_levels": r"SourceAlignedLvl2Param[\s\S]*?\bl\s*=\s*18\s*;",
        "gadget_basebit": r"SourceAlignedLvl2Param[\s\S]*?\bBgbit\s*=\s*2\s*;",
        "error_sigma": r"SourceAlignedLvl2Param[\s\S]*?error_sigma_log2\s*=\s*42\s*;",
        "binary_domain": r"SourceAlignedLvl0Param[\s\S]*?key_value_min\s*=\s*0\s*;",
        "ternary_target": r"SourceAlignedLvl2Param[\s\S]*?key_value_min\s*=\s*-1\s*;",
        "paired_error_marker": (
            r"SourceAlignedLvl2Param[\s\S]*?"
            r"source_aligned_parity_error\s*=\s*true\s*;"
        ),
        "finite_error_marker": (
            r"SourceAlignedLvl2Param[\s\S]*?"
            r"source_aligned_finite_error\s*=\s*true\s*;"
        ),
        "finite_error_coarse_count": (
            r"source_aligned_error_coarse_count\s*=\s*255\s*;"
        ),
        "finite_error_coarse_scale": (
            r"source_aligned_error_coarse_scale_log2\s*=\s*38\s*;"
        ),
        "finite_error_dither_bits": (
            r"source_aligned_error_dither_bits\s*=\s*32\s*;"
        ),
    }
    for description, pattern in bindings.items():
        _require(header_source, pattern, description)

    subset_chain = all(
        fragment in key_source
        for fragment in (
            "std::get<Key<lvl1param>>(keys)[i]",
            "std::get<Key<lvl0param>>(keys)[i]",
            "std::get<Key<lvl2param>>(keys)[i]",
            "std::get<Key<lvl1param>>(keys)[i]",
        )
    )
    parity_key_law = all(
        token in parity_source
        for token in (
            "sourceAlignedPrefixIndex",
            "sourceAlignedSuffixIndex",
            "sourceAlignedParityKeyGen",
            "isSourceAlignedParityKeyPair",
            "sourceAlignedBkfftGen",
        )
    )
    parity_error_law = all(
        token in trlwe_source
        for token in (
            "usesSourceAlignedParityError",
            "sampleSourceAlignedParityErrorPair",
            "joinSourceAlignedParityError",
        )
    )
    finite_error_law = all(
        token in trlwe_source
        for token in (
            "SourceAlignedProofErrorCoins",
            "decodeSourceAlignedProofError",
            "sampleSourceAlignedProofErrorCoins",
            "sampleSourceAlignedProofError",
            "std::popcount",
            "dither_positive",
            "dither_negative",
        )
    )
    exact_parity_reduction = all(
        token in parity_proof_source
        for token in (
            "smallTernaryRLWEProblem",
            "RLWE.problem (2 ^ 64) 1024 (2 * brkRowCount)",
            "paritySuffixAdvantage_eq_smallTernaryRLWE",
            "endpointAdvantage_le_smallRLWE_and_prefix",
            "2 * epsilonRLWE + 4 * epsilonPrefix",
        )
    )
    exact_finite_mgf = all(
        token in finite_error_proof_source
        for token in (
            "def coarseCount : ℕ := 255",
            "def coarseScaleLog : ℕ := 38",
            "def ditherBits : ℕ := 32",
            "expectation_exp_decodeReal_le_two_pow",
            "expectation_exp_dotProduct_noiseVector_le",
            "sphericalCertificate",
        )
    ) and all(
        token in parameter_proof_source
        for token in (
            "finiteSamplerCovariance_eq",
            "finiteFreshNoiseCertificate",
            "finiteSampler_adaptiveAbsTail",
        )
    )
    all_headers = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in (root / "include").rglob("*.hpp")
        if path != header
    )
    native_aligned_format = any(
        token in all_headers
        for token in ("SourceAlignedKSK", "source_aligned_ksk", "AlignedKSK")
    )

    exponent = candidate.chernoff_exponent
    sigma_ceiling = _fresh_sigma_ceiling(candidate, failure_bits)
    input_alpha = 0.000_092_511_997_467_675_6
    input_sigma = input_alpha * (1 << 16)

    return {
        "scope": {
            "statement": (
                "Modified correlated source-aligned lvl02 candidate; not a "
                "claim about TFHEpp's native independent subset KSK."
            ),
            "tfhepp_root": str(root.resolve()),
            "formal_proof_root": str(formal_proof_root.resolve()),
            "failure_target_bits": failure_bits,
            "source_binding_ok": (
                parity_key_law
                and parity_error_law
                and finite_error_law
                and exact_parity_reduction
                and exact_finite_mgf
            ),
            "lvlh2_exclusion": (
                "The current lvlhalf key is sampled independently of lvl2, "
                "so lvlh2 is not an instance of the subset-prefix theorem."
            ),
        },
        "source_evidence": {
            "candidate_header": str(header.relative_to(root)),
            "candidate_header_sha256": _sha256(header),
            "key_header_sha256": _sha256(key_header),
            "parity_header": str(parity_header.relative_to(root)),
            "parity_header_sha256": _sha256(parity_header),
            "trlwe_header_sha256": _sha256(trlwe_header),
            "parity_proof": str(parity_proof.relative_to(formal_proof_root)),
            "parity_proof_sha256": _sha256(parity_proof),
            "finite_error_proof": str(
                finite_error_proof.relative_to(formal_proof_root)
            ),
            "finite_error_proof_sha256": _sha256(finite_error_proof),
            "parameter_proof": str(
                parameter_proof.relative_to(formal_proof_root)
            ),
            "parameter_proof_sha256": _sha256(parameter_proof),
            "default_subset_chain_lvl0_through_lvl1_to_lvl2": subset_chain,
            "proof_aligned_parity_key_law_detected": parity_key_law,
            "proof_aligned_paired_error_law_detected": parity_error_law,
            "exact_finite_error_decoder_detected": finite_error_law,
            "exact_finite_error_mgf_certificate_detected": exact_finite_mgf,
            "exact_ordinary_ternary_rlwe_reduction_detected": (
                exact_parity_reduction
            ),
            "native_correlated_aligned_ksk_detected": native_aligned_format,
            "candidate_registered_in_secret_key": False,
        },
        "candidate": {
            **asdict(candidate),
            "suffix_dimension": candidate.suffix_dimension,
            "rows_per_control": candidate.rows_per_control,
            "brk_rows": candidate.brk_rows,
            "suffix_rlwe_samples": candidate.suffix_rlwe_samples,
            "aligned_width": candidate.aligned_width,
            "digit_radius": candidate.digit_radius,
            "factor_energy_bound": candidate.factor_energy_bound,
            "decomposition_bits": candidate.decomposition_bits,
            "fresh_variance_bound": candidate.fresh_variance_bound,
            "conservative_threshold": candidate.threshold,
            "fresh_sigma_ceiling_for_requested_tail": sigma_ceiling,
            "chernoff_exponent": _fraction(exponent),
            "log2_two_sided_tail_upper_bound": _log2_two_sided_tail(exponent),
            "finite_error_law": {
                "formula": "2^38 * sum_{j=1}^{255} R_j + (U - V)",
                "coarse_count": candidate.finite_error_coarse_count,
                "coarse_scale_log2": (
                    candidate.finite_error_coarse_scale_log2
                ),
                "dither_bits": candidate.finite_error_dither_bits,
                "variance": _fraction(candidate.finite_error_variance),
                "standard_deviation": candidate.finite_error_stddev,
                "proved_mgf_proxy": candidate.finite_error_mgf_proxy,
                "configured_mgf_proxy": candidate.fresh_sigma**2,
                "proved_proxy_below_configured": (
                    candidate.finite_error_mgf_proxy
                    <= candidate.fresh_sigma**2
                ),
                "sharper_chernoff_exponent": _fraction(
                    candidate.finite_error_chernoff_exponent
                ),
                "log2_sharper_two_sided_tail_upper_bound": (
                    _log2_two_sided_tail(
                        candidate.finite_error_chernoff_exponent
                    )
                ),
            },
        },
        "hardness_instances": {
            "prefix_first": {
                "kind": "binary-secret LWE",
                "n": candidate.prefix_dimension,
                "q": 1 << candidate.torus_bits,
                "m": candidate.aligned_width,
                "sigma": candidate.finite_error_stddev,
                "error_model": (
                    "variance-matched discrete-Gaussian attack proxy for the "
                    "exact finite implementation law"
                ),
                "multiplicity": 2,
            },
            "prefix_second": {
                "kind": "binary-secret LWE",
                "n": candidate.prefix_dimension,
                "q": 1 << candidate.torus_bits,
                "m": candidate.aligned_width,
                "sigma": candidate.finite_error_stddev,
                "error_model": (
                    "variance-matched discrete-Gaussian attack proxy for the "
                    "exact finite implementation law"
                ),
                "multiplicity": 2,
            },
            "suffix_ordinary_ternary_rlwe": {
                "kind": "ordinary degree-1024 ternary RLWE",
                "n": candidate.suffix_dimension,
                "ring_degree": candidate.suffix_dimension,
                "ring_rank": 1,
                "q": 1 << candidate.torus_bits,
                "m": candidate.suffix_rlwe_samples,
                "sigma": candidate.finite_error_stddev,
                "error_model": (
                    "variance-matched discrete-Gaussian attack proxy for the "
                    "exact finite implementation law"
                ),
                "multiplicity": 2,
                "formal_reduction": (
                    "Exact parity scalarization and odd-secret reduction to "
                    "ordinary ternary RLWE."
                ),
                "warning": (
                    "The formal source reduction is exact.  Treating this "
                    "RLWE instance as coefficient LWE with a variance-matched "
                    "discrete Gaussian in lattice-estimator is still "
                    "heuristic attack evidence, not a theorem."
                ),
            },
            "input_tlwe": {
                "kind": "binary-secret LWE",
                "n": candidate.prefix_dimension,
                "q": 1 << 16,
                "m": candidate.aligned_width,
                "sigma": input_sigma,
                "warning": "The sample count is deliberately conservative.",
            },
        },
        "checks": {
            "balanced_subset_split": candidate.prefix_dimension
            == candidate.suffix_dimension,
            "same_decomposition_coverage_as_default_lvl2": (
                candidate.decomposition_bits == 36
            ),
            "fresh_sigma_below_exact_component_ceiling": (
                candidate.fresh_sigma <= sigma_ceiling
            ),
            "fresh_component_meets_requested_tail": (
                _log2_two_sided_tail(exponent) <= -failure_bits
            ),
            "proof_aligned_parity_key_law_implemented": parity_key_law,
            "proof_aligned_paired_error_law_implemented": parity_error_law,
            "exact_finite_error_decoder_implemented": finite_error_law,
            "exact_finite_error_mgf_certificate_formalized": exact_finite_mgf,
            "finite_error_proxy_matches_parameter": (
                candidate.finite_error_mgf_proxy <= candidate.fresh_sigma**2
            ),
            "exact_ordinary_ternary_rlwe_reduction_formalized": (
                exact_parity_reduction
            ),
            "modified_format_implemented": native_aligned_format,
        },
        "remaining_obligations": [
            "complete ordinary bootstrap rounding/noise composition",
            (
                "ordinary ternary-RLWE/LWE hardness for the exact finite "
                "error law; the variance-matched estimator output is "
                "heuristic evidence"
            ),
            "widened correlated aligned-KSK generator and evaluator",
            "standard ideal-uniform-bit model for configured CSPRNG outputs",
        ],
        "conclusion": {
            "arithmetic_correctness_candidate": True,
            "native_tfhepp_parameter_certified": False,
            "reason": (
                "The suffix is now exactly ordinary ternary RLWE, and the "
                "TFHEpp key/error samplers match the parity law.  The exact "
                "finite sampler now has the required MGF certificate.  The "
                "isolated correction arithmetic passes, but the correlated "
                "aligned KSK and complete bootstrap noise composition remain "
                "absent."
            ),
        },
    }


def self_test(root: Path, formal_proof_root: Path | None = None) -> None:
    report = analyse(root, formal_proof_root=formal_proof_root)
    candidate = report["candidate"]
    assert candidate["suffix_dimension"] == 1024
    assert candidate["rows_per_control"] == 36
    assert candidate["brk_rows"] == 36_864
    assert candidate["suffix_rlwe_samples"] == 73_728
    assert candidate["aligned_width"] == 75_497_472
    assert candidate["digit_radius"] == 2
    assert candidate["factor_energy_bound"] == 301_989_888
    assert candidate["chernoff_exponent"]["numerator"] == 1024
    assert candidate["chernoff_exponent"]["denominator"] == 9
    finite_error = candidate["finite_error_law"]
    assert finite_error["variance"]["denominator"] == 2
    assert finite_error["proved_proxy_below_configured"]
    assert finite_error["sharper_chernoff_exponent"]["decimal"] > 1024 / 9
    assert report["checks"]["fresh_component_meets_requested_tail"]
    assert report["checks"]["proof_aligned_parity_key_law_implemented"]
    assert report["checks"]["proof_aligned_paired_error_law_implemented"]
    assert report["checks"]["exact_finite_error_decoder_implemented"]
    assert report["checks"]["exact_finite_error_mgf_certificate_formalized"]
    assert report["checks"]["finite_error_proxy_matches_parameter"]
    assert report["checks"]["exact_ordinary_ternary_rlwe_reduction_formalized"]
    assert not report["checks"]["modified_format_implemented"]


def _print_human(report: dict[str, Any]) -> None:
    candidate = report["candidate"]
    print("Correlated source-aligned lvl02 candidate")
    print(
        "  prefix / suffix / ring degree: "
        f"{candidate['prefix_dimension']} / "
        f"{candidate['suffix_dimension']} / {candidate['ring_degree']}"
    )
    print(
        "  gadget levels / basebit / coverage: "
        f"{candidate['gadget_levels']} / {candidate['gadget_basebit']} / "
        f"{candidate['decomposition_bits']} bits"
    )
    print(
        "  BRK rows / suffix RLWE samples: "
        f"{candidate['brk_rows']} / {candidate['suffix_rlwe_samples']}"
    )
    print(
        "  aligned width / factor energy: "
        f"{candidate['aligned_width']} / "
        f"{candidate['factor_energy_bound']}"
    )
    print(f"  correctness proxy sigma: {candidate['fresh_sigma']}")
    print(
        "  exact finite-law standard deviation: "
        f"{candidate['finite_error_law']['standard_deviation']:.6f}"
    )
    print(
        "  fresh-component log2 failure bound: "
        f"{candidate['log2_two_sided_tail_upper_bound']:.3f}"
    )
    print(
        "  proof-aligned key/error samplers: "
        f"{report['checks']['proof_aligned_parity_key_law_implemented']} / "
        f"{report['checks']['proof_aligned_paired_error_law_implemented']}"
    )
    print(
        "  exact ordinary ternary-RLWE reduction: "
        f"{report['checks']['exact_ordinary_ternary_rlwe_reduction_formalized']}"
    )
    print(
        "  exact finite decoder / MGF certificate: "
        f"{report['checks']['exact_finite_error_decoder_implemented']} / "
        f"{report['checks']['exact_finite_error_mgf_certificate_formalized']}"
    )
    print(
        "  native correlated format implemented: "
        f"{report['checks']['modified_format_implemented']}"
    )


def main() -> None:
    default_root = Path(__file__).resolve().parents[3] / "TFHEpp"
    default_formal_proof_root = (
        Path(__file__).resolve().parents[3] / "FormalProof4FHE"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tfhepp-root", type=Path, default=default_root)
    parser.add_argument(
        "--formal-proof-root", type=Path, default=default_formal_proof_root
    )
    parser.add_argument("--failure-bits", type=int, default=128)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test(args.tfhepp_root, args.formal_proof_root)
    report = analyse(
        args.tfhepp_root, args.failure_bits, args.formal_proof_root
    )
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        _print_human(report)
        if args.self_test:
            print("  self-test: PASS")


if __name__ == "__main__":
    main()
