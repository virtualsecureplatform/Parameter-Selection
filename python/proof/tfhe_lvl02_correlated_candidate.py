#!/usr/bin/env python3
"""Source-bound arithmetic screen for a correlated lvl02 candidate.

The candidate is intentionally a modified-format parameter profile.  It uses
the exact correlated BRK/aligned-KSK error law proved in FormalProof4FHE; it is
not executable through TFHEpp's native independent subset KSK.

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


SELECTED = Candidate(
    prefix_dimension=1024,
    ring_degree=2048,
    ring_rank=1,
    gadget_levels=18,
    gadget_basebit=2,
    torus_bits=64,
    fresh_sigma=1 << 42,
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


def analyse(root: Path, failure_bits: int = 128) -> dict[str, Any]:
    if failure_bits <= 0:
        raise ValueError("failure bits must be positive")

    header = root / "include/params/source_aligned_lvl02.hpp"
    key_header = root / "include/tfhe/key.hpp"
    header_source = header.read_text(encoding="utf-8")
    key_source = key_header.read_text(encoding="utf-8")
    candidate = SELECTED

    bindings = {
        "prefix_dimension": r"SourceAlignedLvl0Param[\s\S]*?\bn\s*=\s*1024\s*;",
        "ring_nbit": r"SourceAlignedLvl2Param[\s\S]*?\bnbit\s*=\s*11\s*;",
        "gadget_levels": r"SourceAlignedLvl2Param[\s\S]*?\bl\s*=\s*18\s*;",
        "gadget_basebit": r"SourceAlignedLvl2Param[\s\S]*?\bBgbit\s*=\s*2\s*;",
        "error_sigma": r"SourceAlignedLvl2Param[\s\S]*?error_sigma_log2\s*=\s*42\s*;",
        "binary_domain": r"SourceAlignedLvl0Param[\s\S]*?key_value_min\s*=\s*0\s*;",
        "ternary_target": r"SourceAlignedLvl2Param[\s\S]*?key_value_min\s*=\s*-1\s*;",
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
            "failure_target_bits": failure_bits,
            "source_binding_ok": True,
            "lvlh2_exclusion": (
                "The current lvlhalf key is sampled independently of lvl2, "
                "so lvlh2 is not an instance of the subset-prefix theorem."
            ),
        },
        "source_evidence": {
            "candidate_header": str(header.relative_to(root)),
            "candidate_header_sha256": _sha256(header),
            "key_header_sha256": _sha256(key_header),
            "default_subset_chain_lvl0_through_lvl1_to_lvl2": subset_chain,
            "native_correlated_aligned_ksk_detected": native_aligned_format,
            "candidate_registered_in_secret_key": False,
        },
        "candidate": {
            **asdict(candidate),
            "suffix_dimension": candidate.suffix_dimension,
            "rows_per_control": candidate.rows_per_control,
            "brk_rows": candidate.brk_rows,
            "aligned_width": candidate.aligned_width,
            "digit_radius": candidate.digit_radius,
            "factor_energy_bound": candidate.factor_energy_bound,
            "decomposition_bits": candidate.decomposition_bits,
            "fresh_variance_bound": candidate.fresh_variance_bound,
            "conservative_threshold": candidate.threshold,
            "fresh_sigma_ceiling_for_requested_tail": sigma_ceiling,
            "chernoff_exponent": _fraction(exponent),
            "log2_two_sided_tail_upper_bound": _log2_two_sided_tail(exponent),
        },
        "hardness_instances": {
            "prefix_first": {
                "kind": "binary-secret LWE",
                "n": candidate.prefix_dimension,
                "q": 1 << candidate.torus_bits,
                "m": candidate.aligned_width,
                "sigma": candidate.fresh_sigma,
                "multiplicity": 2,
            },
            "prefix_second": {
                "kind": "binary-secret LWE",
                "n": candidate.prefix_dimension,
                "q": 1 << candidate.torus_bits,
                "m": candidate.aligned_width,
                "sigma": candidate.fresh_sigma,
                "multiplicity": 2,
            },
            "suffix_proxy": {
                "kind": "known-prefix suffix-RLWE as ordinary-LWE proxy",
                "n": candidate.suffix_dimension,
                "q": 1 << candidate.torus_bits,
                "m": candidate.aligned_width,
                "sigma": candidate.fresh_sigma,
                "multiplicity": 2,
                "warning": (
                    "Scalarizing all rotations as ordinary LWE is a "
                    "conservative heuristic proxy, not a reduction theorem."
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
            "modified_format_implemented": native_aligned_format,
        },
        "remaining_obligations": [
            "finite modular-Gaussian MGF certificate",
            "complete ordinary bootstrap rounding/noise composition",
            "formal suffix-RLWE hardness assumption rather than its LWE proxy",
            "widened correlated aligned-KSK generator and evaluator",
        ],
        "conclusion": {
            "arithmetic_correctness_candidate": True,
            "native_tfhepp_parameter_certified": False,
            "reason": (
                "The isolated correction arithmetic passes, but security is "
                "heuristic and the correlated aligned format is not implemented."
            ),
        },
    }


def self_test(root: Path) -> None:
    report = analyse(root)
    candidate = report["candidate"]
    assert candidate["suffix_dimension"] == 1024
    assert candidate["rows_per_control"] == 36
    assert candidate["brk_rows"] == 36_864
    assert candidate["aligned_width"] == 75_497_472
    assert candidate["digit_radius"] == 2
    assert candidate["factor_energy_bound"] == 301_989_888
    assert candidate["chernoff_exponent"]["numerator"] == 1024
    assert candidate["chernoff_exponent"]["denominator"] == 9
    assert report["checks"]["fresh_component_meets_requested_tail"]
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
        "  BRK rows / aligned width / factor energy: "
        f"{candidate['brk_rows']} / {candidate['aligned_width']} / "
        f"{candidate['factor_energy_bound']}"
    )
    print(f"  fresh sigma: {candidate['fresh_sigma']}")
    print(
        "  fresh-component log2 failure bound: "
        f"{candidate['log2_two_sided_tail_upper_bound']:.3f}"
    )
    print(
        "  native correlated format implemented: "
        f"{report['checks']['modified_format_implemented']}"
    )


def main() -> None:
    default_root = Path(__file__).resolve().parents[3] / "TFHEpp"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tfhepp-root", type=Path, default=default_root)
    parser.add_argument("--failure-bits", type=int, default=128)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test(args.tfhepp_root)
    report = analyse(args.tfhepp_root, args.failure_bits)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        _print_human(report)
        if args.self_test:
            print("  self-test: PASS")


if __name__ == "__main__":
    main()
