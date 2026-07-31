#!/usr/bin/env python3
"""Source-bound screen for the correlated source-aligned TFHE theorem.

The screen reads TFHEpp's default parameter and implementation headers.  It
computes the exact non-bundled BRK row count, the coefficient-scalarized
aligned-KSK width, and a conservative centered-digit factor-energy bound.  It
then evaluates the *fresh correction term* left after the formal correlated
BRK-error cancellation.

This is intentionally split into three statuses:

* formal algebra/correctness component proved by the Lean theorem;
* concrete arithmetic derived from the audited source parameters; and
* still-open implementation/security premises.

In particular, the current C++ tree does not generate the widened correlated
aligned KSK.  The lattice instances printed here are inputs to the companion
Sage/lattice-estimator script; this pure-Python screen does not turn those
heuristic estimates into a theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import asdict
from fractions import Fraction
from pathlib import Path
from typing import Any

from tfhe_subset_joint_screen import (
    _integer,
    _read,
    _struct,
    read_parameters,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction_json(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": str(float(value)),
    }


def _all_fragments(source: str, fragments: tuple[str, ...]) -> bool:
    return all(fragment in source for fragment in fragments)


def _floor_correctness_sigma(
    threshold: int, factor_energy: int, failure_bits: int
) -> int:
    denominator = 2.0 * factor_energy * (failure_bits + 1) * math.log(2.0)
    return math.floor(threshold / math.sqrt(denominator))


def _log2_two_sided_tail(exponent: Fraction) -> float:
    """Return log2(2*exp(-exponent)); this is an upper-bound logarithm."""
    return 1.0 - float(exponent) / math.log(2.0)


def analyse(root: Path, failure_bits: int = 128) -> dict[str, Any]:
    if failure_bits <= 0:
        raise ValueError("failure bits must be positive")

    base, base_evidence = read_parameters(root)
    files = {
        "parameter_header": root / "include/params/128bit.hpp",
        "trgsw": root / "include/tfhe/trgsw.hpp",
        "evalkeygens": root / "include/tfhe/evalkeygens.hpp",
        "gatebootstrapping": root / "include/tfhe/gatebootstrapping.hpp",
        "cmake": root / "CMakeLists.txt",
    }
    source = {name: _read(path) for name, path in files.items()}
    lvl1 = _struct(source["parameter_header"], "lvl1param")

    ring_rank = _integer(lvl1, "k")
    levels = _integer(lvl1, "l")
    basebit = _integer(lvl1, "Bgbit")
    l_a_equals_l = bool(re.search(r"\blₐ\s*=\s*l\s*;", lvl1))
    if not l_a_equals_l:
        raise ValueError("screen currently requires lvl1 l_a = l")
    nonce_levels = levels

    # A binary level-zero control has one nonzero plaintext branch.  A rank-k
    # TGSW has k*l_a nonce rows and l body rows.
    encrypted_values_per_control = (
        base.lvl0_key_max - base.lvl0_key_min
        - (1 if base.lvl0_key_min <= 0 <= base.lvl0_key_max else 0)
        + 1
    )
    rows_per_control = ring_rank * nonce_levels + levels
    brk_rows = (
        base.lvl0_dimension * encrypted_values_per_control * rows_per_control
    )
    aligned_width = brk_rows * base.lvl1_dimension
    centered_digit_radius = 1 << (basebit - 1)
    factor_energy = aligned_width * centered_digit_radius**2
    native_subset_width = base.ksk_rows

    lvl1_sigma = 1 << (base.lvl1_torus_bits + base.lvl1_alpha_log2)
    threshold = 1 << (base.lvl1_torus_bits - 4)
    fresh_variance = lvl1_sigma**2 * factor_energy
    fresh_exponent = Fraction(threshold**2, 2 * fresh_variance)
    maximum_sigma = _floor_correctness_sigma(
        threshold, factor_energy, failure_bits
    )
    maximum_sigma_exponent = Fraction(
        threshold**2, 2 * maximum_sigma**2 * factor_energy
    )

    include_headers = list((root / "include").rglob("*.hpp"))
    aligned_tokens = ("SourceAlignedKSK", "source_aligned_ksk", "AlignedKSK")
    aligned_token_hits = []
    for path in include_headers:
        text = _read(path)
        for token in aligned_tokens:
            if token in text:
                aligned_token_hits.append(
                    {"path": str(path.relative_to(root)), "token": token}
                )

    source_evidence = {
        "sha256": {name: _sha256(path) for name, path in files.items()},
        "base_parameter_source_binding": base_evidence,
        "level_one_l_a_equals_l": l_a_equals_l,
        "nonbundled_default": bool(
            re.search(
                r'option\(USE_KEY_BUNDLE\s+"[^"]*"\s+OFF\)',
                source["cmake"],
            )
        ),
        "binary_bk_generation_has_one_nonzero_value_per_control": (
            base.lvl0_key_min == 0
            and base.lvl0_key_max == 1
            and _all_fragments(
                source["evalkeygens"],
                (
                    "P::domainP::key_value_min",
                    "P::domainP::key_value_max",
                    "if (j != 0)",
                    "trgswSymEncrypt<typename P::targetP>",
                ),
            )
        ),
        "centered_base_decomposition_matches_radius": _all_fragments(
            source["trgsw"],
            (
                "maskBg",
                "halfBg",
                "((a >> shift) & maskBg) - halfBg",
            ),
        ),
        "native_subset_ksk_is_suffix_coordinate_layout": _all_fragments(
            source["evalkeygens"],
            (
                "void subikskgen",
                "P::domainP::k * P::domainP::n - P::targetP::k * P::targetP::n",
                "ksk[i][j][k]",
                "tlweSymEncrypt<typename P::targetP>",
            ),
        ),
        "native_subset_dispatch_present": "SubsetIdentityKeySwitch<P>" in source[
            "gatebootstrapping"
        ],
        "aligned_generator_token_hits": aligned_token_hits,
    }

    binding_flags = [
        source_evidence["level_one_l_a_equals_l"],
        source_evidence["nonbundled_default"],
        source_evidence[
            "binary_bk_generation_has_one_nonzero_value_per_control"
        ],
        source_evidence["centered_base_decomposition_matches_radius"],
        source_evidence["native_subset_ksk_is_suffix_coordinate_layout"],
        source_evidence["native_subset_dispatch_present"],
    ]
    source_binding_ok = all(binding_flags)
    aligned_layout_detected = bool(aligned_token_hits)

    prefix_instance = {
        "kind": "binary-secret LWE",
        "n": base.lvl0_dimension,
        "q": 1 << base.lvl1_torus_bits,
        "q_bits": base.lvl1_torus_bits,
        "m": aligned_width,
        "secret": "binary",
        "error_sigma_baseline": lvl1_sigma,
        "error_sigma_largest_integer_passing_fresh_term_screen": maximum_sigma,
        "multiplicity_in_reduction": 2,
        "warning": (
            "The theorem needs a concrete finite error law and a hardness "
            "statement; a lattice-estimator result is heuristic evidence only."
        ),
    }
    suffix_instance = {
        "kind": "known-prefix suffix-RLWE",
        "ring_degree": base.lvl1_dimension,
        "q": 1 << base.lvl1_torus_bits,
        "q_bits": base.lvl1_torus_bits,
        "ring_samples": brk_rows,
        "scalar_equations_if_all_negacyclic_rotations_are_exposed": aligned_width,
        "known_prefix_coefficients": base.lvl0_dimension,
        "unknown_iid_ternary_suffix_coefficients": base.suffix_dimension,
        "error_sigma": lvl1_sigma,
        "multiplicity_in_reduction": 1,
        "warning": (
            "Treating this structured suffix-RLWE instance as ordinary LWE is "
            "a conservative heuristic proxy, not the formal RLWE assumption."
        ),
    }

    return {
        "scope": {
            "statement": (
                "Parameter screen for the modified correlated source-aligned "
                "BRK/KSK theorem; not an insecurity theorem."
            ),
            "tfhepp_root": str(root.resolve()),
            "source_binding_ok": source_binding_ok,
            "failure_target_bits": failure_bits,
        },
        "parameters": {
            **{
                key: value
                for key, value in asdict(base).items()
                if key != "lvl0_alpha"
            },
            "ring_rank": ring_rank,
            "brk_levels": levels,
            "brk_basebit": basebit,
            "encrypted_values_per_control": encrypted_values_per_control,
            "rows_per_control": rows_per_control,
            "brk_rows": brk_rows,
            "aligned_width": aligned_width,
            "native_subset_width": native_subset_width,
            "centered_digit_radius": centered_digit_radius,
            "factor_energy_bound": factor_energy,
            "lvl1_integer_sigma": lvl1_sigma,
            "conservative_fresh_term_threshold": threshold,
        },
        "source_evidence": source_evidence,
        "checks": {
            "formal_correlated_cancellation": {
                "reused_brk_error_after_key_switch": "exactly zero",
                "remaining_random_term": "-<factor,freshError>",
                "result": "PASS_FORMAL_EXACT",
            },
            "adaptive_factor_tail": {
                "factor_may_depend_on_complete_prior_transcript": True,
                "fresh_error_sampled_independently_after_transcript": True,
                "reachable_factor_union_loss": "none",
                "required_analytic_premise": (
                    "finite sampler MGF certificate plus a uniform factor "
                    "covariance-energy bound"
                ),
                "result": "PASS_FORMAL_CONDITIONAL_MGF",
            },
            "fresh_term_arithmetic": {
                "variance_proxy": fresh_variance,
                "chernoff_exponent": _fraction_json(fresh_exponent),
                "log2_two_sided_tail_upper_bound": _log2_two_sided_tail(
                    fresh_exponent
                ),
                "largest_integer_sigma_for_requested_tail": maximum_sigma,
                "largest_sigma_log2_two_sided_tail_upper_bound": (
                    _log2_two_sided_tail(maximum_sigma_exponent)
                ),
                "baseline_sigma_passes": (
                    _log2_two_sided_tail(fresh_exponent) <= -failure_bits
                ),
                "caveat": (
                    "This is only the fresh correction after cancellation; "
                    "ordinary blind-rotation decomposition, input, and output "
                    "rounding terms are not included."
                ),
                "result": "PASS_FRESH_COMPONENT_ONLY",
            },
            "layout": {
                "aligned_width_exceeds_native_subset_width": (
                    aligned_width > native_subset_width
                ),
                "aligned_width_minus_native_subset_width": (
                    aligned_width - native_subset_width
                ),
                "aligned_generator_detected_in_audited_headers": (
                    aligned_layout_detected
                ),
                "result": (
                    "PASS_IMPLEMENTED"
                    if aligned_layout_detected
                    else "FAIL_MODIFIED_FORMAT_NOT_IMPLEMENTED"
                ),
            },
            "finite_sampler_bridge": {
                "cxx_sampler_directly_is_proved_mgf_certificate": False,
                "required": (
                    "an exact finite MGF certificate or an explicit "
                    "implementation-to-certified-law distance"
                ),
                "result": "OPEN",
            },
            "ordinary_bootstrap_correctness_terms": {
                "included_in_this_screen": False,
                "required": (
                    "compose the fresh correction bound with the existing "
                    "input, gadget-decomposition, external-product, modulus-"
                    "switch, extraction, and decoding margins"
                ),
                "result": "OPEN_TECHNICAL_ACCOUNTING",
            },
        },
        "hardness_instances": {
            "suffix": suffix_instance,
            "prefix_first": prefix_instance,
            "prefix_second": prefix_instance,
            "formal_reduction_loss": (
                "2*Adv_suffix + 2*Adv_prefix_first + "
                "2*Adv_prefix_second + epsilon_aux"
            ),
        },
        "conclusion": {
            "current_tfhepp_implements_target_format": aligned_layout_detected,
            "fresh_correction_component_at_baseline": "PASSES_NOMINAL_SCREEN",
            "theorem_certifies_current_parameter": False,
            "reason": (
                "The exact cancellation and nominal fresh-noise margin now "
                "pass, but the widened correlated aligned KSK is absent; "
                "the two prefix/suffix hardness instances, finite sampler "
                "certificate, and full bootstrap noise composition remain."
            ),
        },
    }


def self_test(root: Path) -> None:
    report = analyse(root)
    p = report["parameters"]
    c = report["checks"]
    assert report["scope"]["source_binding_ok"]
    assert p["rows_per_control"] == 6
    assert p["brk_rows"] == 3780
    assert p["aligned_width"] == 3_870_720
    assert p["native_subset_width"] == 5516
    assert p["centered_digit_radius"] == 32
    assert p["factor_energy_bound"] == 3_963_617_280
    assert p["lvl1_integer_sigma"] == 128
    assert c["fresh_term_arithmetic"]["variance_proxy"] == 64_939_905_515_520
    exponent = c["fresh_term_arithmetic"]["chernoff_exponent"]
    assert exponent["numerator"] == "524288"
    assert exponent["denominator"] == "945"
    assert c["fresh_term_arithmetic"]["baseline_sigma_passes"]
    assert c["fresh_term_arithmetic"][
        "largest_integer_sigma_for_requested_tail"
    ] == 318
    assert c["layout"]["result"] == "FAIL_MODIFIED_FORMAT_NOT_IMPLEMENTED"
    assert not report["conclusion"]["theorem_certifies_current_parameter"]


def _print_human(report: dict[str, Any]) -> None:
    p = report["parameters"]
    c = report["checks"]
    h = report["hardness_instances"]
    print("TFHEpp correlated source-aligned theorem screen")
    print(f"  source binding: {report['scope']['source_binding_ok']}")
    print(
        "  BRK rows / aligned width / native subset width: "
        f"{p['brk_rows']} / {p['aligned_width']} / {p['native_subset_width']}"
    )
    print(f"  correlated cancellation: {c['formal_correlated_cancellation']['result']}")
    print(f"  transcript-adaptive tail: {c['adaptive_factor_tail']['result']}")
    print(
        "  baseline fresh-term log2 failure bound: "
        f"{c['fresh_term_arithmetic']['log2_two_sided_tail_upper_bound']:.3f}"
    )
    print(
        "  largest integer fresh sigma for requested component tail: "
        f"{c['fresh_term_arithmetic']['largest_integer_sigma_for_requested_tail']}"
    )
    print(f"  layout: {c['layout']['result']}")
    print(f"  finite sampler bridge: {c['finite_sampler_bridge']['result']}")
    print(
        "  prefix LWE input: "
        f"n={h['prefix_first']['n']}, q=2^{h['prefix_first']['q_bits']}, "
        f"m={h['prefix_first']['m']}, sigma="
        f"{h['prefix_first']['error_sigma_baseline']}.."
        f"{h['prefix_first']['error_sigma_largest_integer_passing_fresh_term_screen']}"
    )
    print(
        "  suffix proxy input: "
        f"unknown={h['suffix']['unknown_iid_ternary_suffix_coefficients']}, "
        f"q=2^{h['suffix']['q_bits']}, ring samples={h['suffix']['ring_samples']}"
    )
    print(
        "  theorem certifies current parameter: "
        f"{report['conclusion']['theorem_certifies_current_parameter']}"
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
