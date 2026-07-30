#!/usr/bin/env python3
"""Source-bound screen for TFHEpp's subset-key joint BRK/KSK theorem route.

The checker reads the default 128-bit TFHEpp parameters and the implementation
sites that generate the subset secret and subset key-switch key.  It then
evaluates the equal-covariance obligations and the delayed-projection fallback:

* whether the implementation error sampler is already covered by the exact
  continuous-Gaussian pair-kernel certificate;
* whether the centered-mixture interaction condition can hold on the
  covariance-compatible branch;
* whether positive-semidefinite covariance correction is compatible with a
  high-probability constrained factorization;
* whether combining at the large modulus and projecting once satisfies the
  exact word-conversion algebra and its image-aware, valuation-stratified
  two-budget count.

This is a screen of one sufficient security proof, not an attack and not an
insecurity claim.  Exact arithmetic uses ``int`` and ``Fraction`` throughout;
floating-point values are reporting aids only.
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
class Parameters:
    lvl0_dimension: int
    lvl1_dimension: int
    lvl0_torus_bits: int
    lvl1_torus_bits: int
    lvl0_key_min: int
    lvl0_key_max: int
    lvl1_key_min: int
    lvl1_key_max: int
    lvl0_alpha: Fraction
    lvl1_alpha_log2: int
    ksk_levels: int
    ksk_basebit: int

    @property
    def suffix_dimension(self) -> int:
        return self.lvl1_dimension - self.lvl0_dimension

    @property
    def ksk_digit_count(self) -> int:
        return 1 << (self.ksk_basebit - 1)

    @property
    def ksk_rows(self) -> int:
        return self.suffix_dimension * self.ksk_levels * self.ksk_digit_count

    @property
    def top_gadget(self) -> int:
        return 1 << (self.lvl0_torus_bits - self.ksk_basebit)

    @property
    def nominal_integer_variance(self) -> Fraction:
        sigma = self.lvl0_alpha * (1 << self.lvl0_torus_bits)
        return sigma * sigma


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError as error:
        raise ValueError(f"cannot read required TFHEpp source: {path}") from error


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _struct(source: str, name: str) -> str:
    match = re.search(rf"\bstruct\s+{re.escape(name)}\s*\{{", source)
    if match is None:
        raise ValueError(f"missing struct {name}")
    start = match.end()
    depth = 1
    for index in range(start, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[start:index]
    raise ValueError(f"unterminated struct {name}")


def _function(source: str, name: str) -> str:
    match = re.search(rf"\b{re.escape(name)}\s*\(", source)
    if match is None:
        raise ValueError(f"missing function {name}")
    start = source.find("{", match.end())
    if start < 0:
        raise ValueError(f"missing function body {name}")
    depth = 1
    for index in range(start + 1, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[start + 1 : index]
    raise ValueError(f"unterminated function {name}")


def _integer(body: str, name: str) -> int:
    match = re.search(rf"\b{re.escape(name)}\s*=\s*(-?\d+)\s*;", body)
    if match is None:
        raise ValueError(f"cannot parse integer field {name}")
    return int(match.group(1))


def _torus_bits(body: str) -> int:
    match = re.search(r"\busing\s+T\s*=\s*uint(16|32|64)_t\s*;", body)
    if match is None:
        raise ValueError("cannot parse torus word type")
    return int(match.group(1))


def _decimal(body: str, name: str) -> Fraction:
    match = re.search(
        rf"\b{re.escape(name)}\s*=\s*([0-9][0-9'.]*)\s*;", body
    )
    if match is None:
        raise ValueError(f"cannot parse decimal field {name}")
    return Fraction(match.group(1).replace("'", ""))


def _pow2_exponent(body: str, name: str) -> int:
    match = re.search(
        rf"\b{re.escape(name)}\s*=\s*std::pow\(2\.0,\s*(-?\d+)\)\s*;",
        body,
    )
    if match is None:
        raise ValueError(f"cannot parse power-of-two field {name}")
    return int(match.group(1))


def _fraction_json(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "decimal": str(float(value)),
    }


def _two_adic_valuation(value: int) -> int:
    """Return the exponent of two in a positive integer."""
    if value <= 0:
        raise ValueError("two-adic valuation requires a positive integer")
    valuation = 0
    while value % 2 == 0:
        value //= 2
        valuation += 1
    return valuation


def _cache_subset_state(root: Path) -> str:
    caches = sorted(root.glob("build*/CMakeCache.txt"))
    states: list[str] = []
    for cache in caches:
        text = _read(cache)
        match = re.search(r"^USE_SUBSET_KEY:BOOL=(ON|OFF)$", text, re.MULTILINE)
        if match is not None:
            states.append(f"{cache.relative_to(root)}:{match.group(1)}")
    return ",".join(states) if states else "NO_CACHE_FOUND"


def read_parameters(root: Path) -> tuple[Parameters, dict[str, Any]]:
    files = {
        "params": root / "include/params.hpp",
        "parameter_header": root / "include/params/128bit.hpp",
        "key": root / "include/tfhe/key.hpp",
        "evalkeygens": root / "include/tfhe/evalkeygens.hpp",
        "gatebootstrapping": root / "include/tfhe/gatebootstrapping.hpp",
        "keyswitch": root / "include/tfhe/keyswitch.hpp",
        "tlwe": root / "include/tfhe/tlwe.hpp",
        "utils": root / "include/utils.hpp",
        "cmake": root / "CMakeLists.txt",
    }
    text = {name: _read(path) for name, path in files.items()}

    lvl0 = _struct(text["parameter_header"], "lvl0param")
    lvl1 = _struct(text["parameter_header"], "lvl1param")
    lvl10 = _struct(text["parameter_header"], "lvl10param")
    subset_iks = _function(text["keyswitch"], "SubsetIdentityKeySwitch")
    lvl1_nbit = _integer(lvl1, "nbit")
    lvl1_dimension = 1 << lvl1_nbit

    parameters = Parameters(
        lvl0_dimension=_integer(lvl0, "n"),
        lvl1_dimension=lvl1_dimension,
        lvl0_torus_bits=_torus_bits(lvl0),
        lvl1_torus_bits=_torus_bits(lvl1),
        lvl0_key_min=_integer(lvl0, "key_value_min"),
        lvl0_key_max=_integer(lvl0, "key_value_max"),
        lvl1_key_min=_integer(lvl1, "key_value_min"),
        lvl1_key_max=_integer(lvl1, "key_value_max"),
        lvl0_alpha=_decimal(lvl0, "α"),
        lvl1_alpha_log2=_pow2_exponent(lvl1, "α"),
        ksk_levels=_integer(lvl10, "t"),
        ksk_basebit=_integer(lvl10, "basebit"),
    )

    evidence = {
        "sha256": {name: _sha256(path) for name, path in files.items()},
        "default_parameter_header_is_128bit": bool(
            re.search(r"#else\s*#include \"params/128bit\.hpp\"", text["params"])
        ),
        "subset_key_cmake_default_off": bool(
            re.search(
                r'option\(USE_SUBSET_KEY\s+"Use subset key"\s+OFF\)',
                text["cmake"],
            )
        ),
        "detected_build_cache_state": _cache_subset_state(root),
        "lvl10_domain_and_target": (
            "using domainP = lvl1param;" in lvl10
            and "using targetP = lvl0param;" in lvl10
        ),
        "lvl10_uses_lvl0_alpha": "α = lvl0param::α" in lvl10,
        "lvl0_sampled_before_lvl1": (
            text["key"].find("keyGen<lvl0param>")
            < text["key"].find("keyGen<lvl1param>")
        ),
        "subset_prefix_overwrite": all(
            fragment in text["key"]
            for fragment in (
                "#ifdef USE_SUBSET_KEY",
                "i < lvl0param::k * lvl0param::n",
                "std::get<Key<lvl1param>>(keys)[i] =",
                "std::get<Key<lvl0param>>(keys)[i]",
            )
        ),
        "general_keygen_uniform_integer_interval": all(
            fragment in text["key"]
            for fragment in (
                "uniform_int_distribution<int32_t> gen(P::key_value_min,",
                "P::key_value_max)",
                "for (typename P::T& i : key) i = gen(generator);",
            )
        ),
        "subset_ksk_suffix_loop": all(
            fragment in text["evalkeygens"]
            for fragment in (
                "P::domainP::k * P::domainP::n - P::targetP::k * P::targetP::n",
                "i < P::t",
                "1U << (P::basebit - 1)",
            )
        ),
        "subset_ksk_encrypts_under_prefix": all(
            fragment in text["evalkeygens"]
            for fragment in (
                "Key<typename P::targetP> subkey;",
                "subkey[i] = domainkey[i];",
                "tlweSymEncrypt<typename P::targetP>",
            )
        ),
        "subset_ksk_gadget_uses_target_word_bits": all(
            fragment in text["evalkeygens"]
            for fragment in (
                "domainkey[P::targetP::k * P::targetP::n + i] * (k + 1)",
                "static_cast<typename P::targetP::T>(1)",
                "numeric_limits<typename P::targetP::T>::digits -",
                "(j + 1) * P::basebit",
            )
        ),
        "subset_dispatch_guarded_by_macro": all(
            fragment in text["gatebootstrapping"]
            for fragment in (
                "#ifdef USE_SUBSET_KEY",
                "SubsetIdentityKeySwitch<P>",
                "std::is_same_v<P, lvl10param>",
            )
        ),
        "target_encryptor_uses_parameter_alpha": all(
            fragment in text["tlwe"]
            for fragment in (
                "ModularGaussian<P>(p, α)",
                "tlweSymEncrypt<P>(res, p, P::α, key)",
            )
        ),
        "subset_word_conversion_uses_half_unit_rounding": all(
            fragment in subset_iks
            for fragment in (
                "domain_digit - target_digit - 1",
                "(domain_digit - target_digit);",
            )
        ),
        "uint16_error_uses_cpp_normal": all(
            fragment in text["utils"]
            for fragment in (
                "std::is_same_v<typename P::T, uint16_t>",
                "std::normal_distribution<double> distribution(0., stdev);",
                "double err = distribution(generator);",
                "return center + dtot16(err);",
            )
        ),
        "dtot16_truncates_to_integer_torus": all(
            fragment in text["utils"]
            for fragment in (
                "inline uint16_t dtot16(double d)",
                "int64_t((d - int64_t(d)) * (1LL << 16))",
            )
        ),
    }
    return parameters, evidence


def analyse(root: Path, source_row_bound_bits: int = 127) -> dict[str, Any]:
    if source_row_bound_bits < 0:
        raise ValueError("source-row bound bits must be nonnegative")
    parameters, evidence = read_parameters(root)

    expected = {
        "lvl0_dimension": 630,
        "lvl1_dimension": 1024,
        "lvl0_torus_bits": 16,
        "lvl1_torus_bits": 32,
        "lvl0_key_min": 0,
        "lvl0_key_max": 1,
        "lvl1_key_min": -1,
        "lvl1_key_max": 1,
        "lvl1_alpha_log2": -25,
        "ksk_levels": 7,
        "ksk_basebit": 2,
    }
    observed = asdict(parameters)
    observed_without_alpha = {
        key: value for key, value in observed.items() if key != "lvl0_alpha"
    }
    expected_alpha = Fraction("0.0000925119974676756")
    parameter_match = (
        observed_without_alpha == expected and parameters.lvl0_alpha == expected_alpha
    )

    iid_ternary_suffix = (
        parameters.lvl1_key_min == -1
        and parameters.lvl1_key_max == 1
        and evidence["lvl0_sampled_before_lvl1"]
        and evidence["subset_prefix_overwrite"]
        and evidence["general_keygen_uniform_integer_interval"]
    )

    secret_variance = Fraction(2, 3)
    nominal_variance = parameters.nominal_integer_variance
    zero_row_top_residual_variance = (
        secret_variance * parameters.top_gadget * parameters.top_gadget
    )
    residual_to_target_ratio = zero_row_top_residual_variance / nominal_variance

    modulus_scale_bits = parameters.lvl1_torus_bits - parameters.lvl0_torus_bits
    modulus_scale = 1 << modulus_scale_bits
    lvl1_integer_variance = Fraction(
        1 << (2 * (parameters.lvl1_torus_bits + parameters.lvl1_alpha_log2)),
        1,
    )
    delayed_effective_source_variance = (
        lvl1_integer_variance / (modulus_scale * modulus_scale)
    )
    delayed_source_only_operator_norm_squared_ceiling = (
        nominal_variance / delayed_effective_source_variance
    )
    delayed_source_only_operator_norm_integer_radius = math.isqrt(
        delayed_source_only_operator_norm_squared_ceiling.numerator
        // delayed_source_only_operator_norm_squared_ceiling.denominator
    )
    lifted_minor_minimum_high_modulus_norm = modulus_scale
    lifted_minor_minimum_derived_variance = lvl1_integer_variance
    lifted_minor_covariance_compatible = (
        lifted_minor_minimum_derived_variance <= nominal_variance
    )
    delayed_constraint_entropy_bits = (
        parameters.lvl1_torus_bits * parameters.suffix_dimension
    )
    delayed_combined_candidate_bits_at_target = (
        delayed_constraint_entropy_bits - 128
    )
    valuation_strata = []
    for valuation in range(parameters.lvl1_torus_bits + 1):
        image_entropy_bits = (
            parameters.lvl1_torus_bits - valuation
        ) * parameters.suffix_dimension
        valuation_strata.append(
            {
                "valuation": valuation,
                "row_combination_image_entropy_bits": image_entropy_bits,
                "combined_compatible_candidate_bits_for_128_bit_bound": (
                    image_entropy_bits - 128
                    if image_entropy_bits >= 128
                    else None
                ),
            }
        )

    lifted_gadget_strata = []
    for level in range(parameters.ksk_levels):
        base_small_valuation = (
            parameters.lvl0_torus_bits
            - (level + 1) * parameters.ksk_basebit
        )
        if base_small_valuation < 0:
            raise ValueError("key-switch gadget level exceeds target word width")
        for digit in range(1, parameters.ksk_digit_count + 1):
            small_valuation = base_small_valuation + _two_adic_valuation(digit)
            lifted_valuation = small_valuation + modulus_scale_bits
            if lifted_valuation > parameters.lvl1_torus_bits:
                raise ValueError("lifted gadget coefficient vanishes modulo source word")
            image_entropy_bits = (
                parameters.lvl1_torus_bits - lifted_valuation
            ) * parameters.suffix_dimension
            lifted_gadget_strata.append(
                {
                    "level": level,
                    "digit_multiplier": digit,
                    "target_modulus_valuation": small_valuation,
                    "lifted_modulus_valuation": lifted_valuation,
                    "row_combination_image_entropy_bits": image_entropy_bits,
                    "combined_compatible_candidate_bits_for_128_bit_bound": (
                        image_entropy_bits - 128
                        if image_entropy_bits >= 128
                        else None
                    ),
                }
            )

    signed_candidate_bits = source_row_bound_bits + 1
    factorization_success_exponent = (
        parameters.lvl0_torus_bits * parameters.suffix_dimension
        - signed_candidate_bits
    )
    if factorization_success_exponent <= 0:
        raise ValueError("source-row bound is too large for a nontrivial count")

    source_binding_flags = [
        value
        for key, value in evidence.items()
        if key not in {"sha256", "detected_build_cache_state"}
        and key != "subset_key_cmake_default_off"
    ]
    source_binding_ok = parameter_match and all(source_binding_flags)
    cache_state = evidence["detected_build_cache_state"]
    subset_active_in_detected_build = ":ON" in cache_state

    return {
        "scope": {
            "statement": (
                "Necessary-condition screen for the centered-mixture joint "
                "subset-key BRK/KSK reduction; not an insecurity claim."
            ),
            "tfhepp_root": str(root.resolve()),
            "source_binding_ok": source_binding_ok,
            "subset_key_active_in_detected_build": subset_active_in_detected_build,
        },
        "parameters": {
            **observed_without_alpha,
            "lvl0_alpha": _fraction_json(parameters.lvl0_alpha),
            "suffix_dimension": parameters.suffix_dimension,
            "ksk_digit_count": parameters.ksk_digit_count,
            "ksk_rows": parameters.ksk_rows,
            "top_gadget": parameters.top_gadget,
            "source_row_bound_bits": source_row_bound_bits,
            "signed_candidate_bound_bits": signed_candidate_bits,
        },
        "source_evidence": evidence,
        "checks": {
            "secret_law": {
                "iid_uniform_ternary_suffix_when_subset_enabled": iid_ternary_suffix,
                "variance": _fraction_json(secret_variance),
                "result": "PASS" if iid_ternary_suffix else "FAIL_SOURCE_BINDING",
            },
            "implementation_error_law": {
                "nominal_integer_variance": _fraction_json(nominal_variance),
                "is_cpp_normal_then_dtot16": (
                    evidence["uint16_error_uses_cpp_normal"]
                    and evidence["dtot16_truncates_to_integer_torus"]
                ),
                "directly_is_exact_continuous_gaussian": False,
                "required_next_certificate": (
                    "exact finite pair-kernel model or an explicit "
                    "implementation-to-reference distance"
                ),
                "result": "OPEN_EXACT_IMPLEMENTATION_PAIR_KERNEL",
            },
            "interaction": {
                "covariance_compatible_integral_nonzero_branch_residual": "zero",
                "mahalanobis_gram": "zero",
                "maximum_absolute_pair_interaction": "0",
                "result": "PASS_ON_COVARIANCE_COMPATIBLE_BRANCH",
            },
            "covariance": {
                "nominal_target_variance": _fraction_json(nominal_variance),
                "zero_row_top_residual_variance": _fraction_json(
                    zero_row_top_residual_variance
                ),
                "residual_to_target_variance_ratio": _fraction_json(
                    residual_to_target_ratio
                ),
                "zero_top_row_psd_compatible": False,
                "nonzero_integral_row_shape": "one signed selector",
                "nonzero_integral_row_residual": "zero",
                "result": "FAIL_ZERO_ROW_AND_FORCE_EXACT_SIGNED_SELECTOR",
            },
            "factorization_probability": {
                "fixed_row_success_upper_bound": (
                    f"2^-{factorization_success_exponent}"
                ),
                "exponent": factorization_success_exponent,
                "complete_factorization_success_is_no_larger": True,
                "result": "FAIL_HIGH_PROBABILITY_SIGNED_SELECTOR_FACTORIZATION",
            },
            "delayed_projection": {
                "large_modulus_bits": parameters.lvl1_torus_bits,
                "small_modulus_bits": parameters.lvl0_torus_bits,
                "scale_bits": modulus_scale_bits,
                "scale": modulus_scale,
                "lvl1_nominal_integer_variance": _fraction_json(
                    lvl1_integer_variance
                ),
                "continuous_proxy_effective_source_variance_after_scaling": (
                    _fraction_json(delayed_effective_source_variance)
                ),
                "continuous_proxy_source_only_operator_norm_squared_ceiling": _fraction_json(
                    delayed_source_only_operator_norm_squared_ceiling
                ),
                "continuous_proxy_source_only_operator_norm_ceiling_decimal": str(
                    float(delayed_source_only_operator_norm_squared_ceiling) ** 0.5
                ),
                "exact_source_only_operator_norm_integer_radius": (
                    delayed_source_only_operator_norm_integer_radius
                ),
                "next_integer_radius_exceeds_budget": (
                    Fraction(
                        lvl1_integer_variance
                        * (delayed_source_only_operator_norm_integer_radius + 1) ** 2,
                        modulus_scale * modulus_scale,
                    )
                    > nominal_variance
                ),
                "one_row_primitive_constraint_entropy_bits": (
                    delayed_constraint_entropy_bits
                ),
                "primitive_combined_candidate_bits_giving_128_bit_upper_bound": (
                    delayed_combined_candidate_bits_at_target
                ),
                "valuation_strata": valuation_strata,
                "implementation_lifted_gadget_strata": lifted_gadget_strata,
                "algebra": "PASS_FORMAL_DELAYED_ROUNDED_HIGH_WORD_IDENTITY",
                "counting": "PASS_FORMAL_PRIMITIVE_TWO_BUDGET_BOUND",
                "image_aware_counting": (
                    "PASS_FORMAL_EXACT_ROW_IMAGE_AND_COMPATIBLE_RESIDUAL_BOUND"
                ),
                "nonprimitive_2adic_strata": (
                    "PASS_FORMAL_FIXED_AND_MIXED_POWER_OF_TWO_STRATUM_BOUNDS"
                ),
                "valuation_stratum_compatible_residual_counts": (
                    "CLOSED_FOR_EXACT_SOLVER_BY_ZERO_RESIDUAL"
                ),
                "candidate_row_stratum_witnesses": (
                    "PASS_FORMAL_AUTOMATIC_LEAST_POWER_OF_TWO_STRATUM"
                ),
                "invertible_minor_solver": (
                    "PASS_FORMAL_EXACT_DISJOINT_BLOCK_SOLVER"
                ),
                "invertible_minor_minimum_high_modulus_norm": (
                    lifted_minor_minimum_high_modulus_norm
                ),
                "invertible_minor_minimum_derived_variance": _fraction_json(
                    lifted_minor_minimum_derived_variance
                ),
                "invertible_minor_covariance_compatible": (
                    lifted_minor_covariance_compatible
                ),
                "invertible_minor_noise_budget": (
                    "PASS"
                    if lifted_minor_covariance_compatible
                    else "FAIL_SOURCE_VARIANCE_ALREADY_EXCEEDS_TARGET"
                ),
                "full_matrix_operator_norm_control": (
                    "PASS_FORMAL_FOR_INVERTIBLE_MINOR_SOLVER"
                ),
                "genuinely_short_high_modulus_preimage_solver": (
                    "OPEN_RESEARCH_INHOMOGENEOUS_SIS"
                ),
                "exact_finite_joint_error_convolution": (
                    "PASS_FORMAL_EXACT_FINITE_MASS_TABLE_AND_TV_BRIDGE"
                ),
                "concrete_cpp_sampler_table_equality_or_distance": "OPEN",
                "result": (
                    "FAIL_CURRENT_INVERTIBLE_MINOR_NOISE_BUDGET_"
                    "SHORT_HIGH_MODULUS_PREIMAGE_OPEN"
                ),
            },
        },
        "conclusion": {
            "detected_build": (
                "SUBSET_ROUTE_ACTIVE"
                if subset_active_in_detected_build
                else "SUBSET_ROUTE_NOT_ACTIVE"
            ),
            "if_compiled_with_subset_key": (
                "NOT_CERTIFIED_CURRENT_MINOR_SOLVER_FAILS_SHORT_PREIMAGE_OPEN"
            ),
            "reason": (
                "The original covariance route forces signed selectors with one-row "
                f"success at most 2^-{factorization_success_exponent}; the exact "
                "invertible-minor delayed solver exists but its minimum nonzero derived "
                "source variance already exceeds the target variance."
            ),
            "remaining_route": (
                "Find a genuinely short 32-bit preimage L with L*A = 2^16*G "
                f"and row norm at most {delayed_source_only_operator_norm_integer_radius}; "
                "the proved invertible-minor solver instead lifts a nonzero target-ring row "
                "and already exceeds the target covariance.  Then instantiate the finite "
                "mass-table equality or an explicit distance for the concrete C++ sampler."
            ),
        },
    }


def self_test(root: Path) -> None:
    report = analyse(root)
    parameters = report["parameters"]
    checks = report["checks"]
    assert report["scope"]["source_binding_ok"]
    assert parameters["suffix_dimension"] == 394
    assert parameters["ksk_digit_count"] == 2
    assert parameters["ksk_rows"] == 5516
    assert parameters["top_gadget"] == 1 << 14
    assert checks["secret_law"]["iid_uniform_ternary_suffix_when_subset_enabled"]
    assert Fraction(
        checks["covariance"]["nominal_target_variance"]["numerator"]
    ) / Fraction(
        checks["covariance"]["nominal_target_variance"]["denominator"]
    ) < 37
    assert Fraction(
        checks["covariance"]["zero_row_top_residual_variance"]["numerator"]
    ) / Fraction(
        checks["covariance"]["zero_row_top_residual_variance"]["denominator"]
    ) > 37
    assert checks["factorization_probability"]["exponent"] == 6176
    assert checks["delayed_projection"]["scale_bits"] == 16
    assert checks["delayed_projection"]["lvl1_nominal_integer_variance"] == {
        "numerator": "16384",
        "denominator": "1",
        "decimal": "16384.0",
    }
    assert (
        checks["delayed_projection"][
            "exact_source_only_operator_norm_integer_radius"
        ]
        == 3104
    )
    assert checks["delayed_projection"]["next_integer_radius_exceeds_budget"]
    assert not checks["delayed_projection"][
        "invertible_minor_covariance_compatible"
    ]
    assert (
        checks["delayed_projection"]["invertible_minor_noise_budget"]
        == "FAIL_SOURCE_VARIANCE_ALREADY_EXCEEDS_TARGET"
    )
    assert (
        checks["delayed_projection"][
            "one_row_primitive_constraint_entropy_bits"
        ]
        == 12608
    )
    assert (
        checks["delayed_projection"][
            "primitive_combined_candidate_bits_giving_128_bit_upper_bound"
        ]
        == 12480
    )
    strata = checks["delayed_projection"]["valuation_strata"]
    assert len(strata) == 33
    assert strata[0]["row_combination_image_entropy_bits"] == 12608
    assert strata[-1]["row_combination_image_entropy_bits"] == 0
    assert strata[-1]["combined_compatible_candidate_bits_for_128_bit_bound"] is None
    gadget_strata = checks["delayed_projection"][
        "implementation_lifted_gadget_strata"
    ]
    assert len(gadget_strata) == 14
    assert gadget_strata[0] == {
        "level": 0,
        "digit_multiplier": 1,
        "target_modulus_valuation": 14,
        "lifted_modulus_valuation": 30,
        "row_combination_image_entropy_bits": 788,
        "combined_compatible_candidate_bits_for_128_bit_bound": 660,
    }
    assert gadget_strata[-1]["lifted_modulus_valuation"] == 19
    assert (
        checks["delayed_projection"]["nonprimitive_2adic_strata"]
        == "PASS_FORMAL_FIXED_AND_MIXED_POWER_OF_TWO_STRATUM_BOUNDS"
    )
    assert (
        checks["delayed_projection"]["candidate_row_stratum_witnesses"]
        == "PASS_FORMAL_AUTOMATIC_LEAST_POWER_OF_TWO_STRATUM"
    )
    assert (
        checks["delayed_projection"]["exact_finite_joint_error_convolution"]
        == "PASS_FORMAL_EXACT_FINITE_MASS_TABLE_AND_TV_BRIDGE"
    )


def _print_human(report: dict[str, Any]) -> None:
    parameters = report["parameters"]
    checks = report["checks"]
    print("TFHEpp subset-key joint centered-mixture screen")
    print(f"  source binding: {report['scope']['source_binding_ok']}")
    print(
        "  detected build subset route: "
        f"{report['scope']['subset_key_active_in_detected_build']}"
    )
    print(
        "  suffix / KSK rows / top gadget: "
        f"{parameters['suffix_dimension']} / {parameters['ksk_rows']} / "
        f"2^{parameters['lvl0_torus_bits'] - parameters['ksk_basebit']}"
    )
    print(f"  implementation pair law: {checks['implementation_error_law']['result']}")
    print(f"  interaction: {checks['interaction']['result']}")
    print(f"  covariance: {checks['covariance']['result']}")
    print(
        "  signed-selector factorization: "
        f"{checks['factorization_probability']['result']} "
        f"({checks['factorization_probability']['fixed_row_success_upper_bound']})"
    )
    print(
        "  delayed-projection route: "
        f"{checks['delayed_projection']['result']}"
    )
    print(
        "  short high-modulus row radius: "
        f"{checks['delayed_projection']['exact_source_only_operator_norm_integer_radius']}"
    )
    print(
        "  invertible-minor noise budget: "
        f"{checks['delayed_projection']['invertible_minor_noise_budget']}"
    )
    print(
        "  valuation-aware count: "
        f"{checks['delayed_projection']['nonprimitive_2adic_strata']}"
    )
    print(
        "  compatible residual counts: "
        f"{checks['delayed_projection']['valuation_stratum_compatible_residual_counts']}"
    )
    print(
        "  finite rounded-error bridge: "
        f"{checks['delayed_projection']['exact_finite_joint_error_convolution']}"
    )
    print(f"  conclusion: {report['conclusion']['if_compiled_with_subset_key']}")


def main() -> None:
    default_root = Path(__file__).resolve().parents[3] / "TFHEpp"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tfhepp-root", type=Path, default=default_root)
    parser.add_argument("--source-row-bound-bits", type=int, default=127)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test(args.tfhepp_root)
    report = analyse(args.tfhepp_root, args.source_row_bound_bits)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        _print_human(report)
        if args.self_test:
            print("  self-test: PASS")


if __name__ == "__main__":
    main()
