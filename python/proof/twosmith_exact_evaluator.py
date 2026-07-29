#!/usr/bin/env python3
"""Exact certificates for the two-level Smith/Hasse overlap formula.

All probability arithmetic uses :class:`fractions.Fraction`.  JSON decimal
literals are rejected instead of being silently rounded.  The commands match
``sketch/twosmith_exact_parameter_note.tex``:

* ``error-table`` computes ``p_e``, ``beta_e`` and ``zeta_hat_e``;
* ``local`` evaluates the Lucas recursion and the two-level product bound;
* ``secret-hist`` computes selected fixed-weight ternary pair strata;
* ``ord`` computes the exact ``(1-X)``-adic order in the power-of-two ring;
* ``descriptor`` evaluates ``(t-s) * (z + f*g*(t+s))``;
* ``aggregate`` sums a full joint row-tuple histogram.

The emitted recursion traces and valuation decompositions are intended as
small, exact inputs to a Lean-side certificate checker.  This program is not
a replacement for the theorems: it only selects and serializes their finite
parameters.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence


def _reject_float(token: str) -> None:
    raise ValueError(
        f"floating-point JSON value {token!r} is forbidden; use an integer "
        "or a rational string such as \"17/32\""
    )


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(
            handle,
            parse_float=_reject_float,
            parse_constant=_reject_float,
        )


def parse_fraction(value: Any, *, field: str = "value") -> Fraction:
    """Parse an exact JSON integer or rational string.

    JSON booleans are rejected even though ``bool`` subclasses ``int`` in
    Python.  Decimal and exponential strings are rejected by ``Fraction``'s
    explicit syntax check below.
    """

    if isinstance(value, bool):
        raise ValueError(f"{field} must be an integer or rational string")
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        token = value.strip()
        if not token:
            raise ValueError(f"{field} is empty")
        pieces = token.split("/")
        if len(pieces) > 2 or any(
            not piece or piece.lstrip("+-").isdigit() is False
            for piece in pieces
        ):
            raise ValueError(
                f"{field}={value!r} is not an integer or a/b rational"
            )
        try:
            return Fraction(token)
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError(f"invalid exact rational for {field}: {value!r}") from error
    raise ValueError(
        f"{field} must be a JSON integer or rational string, not "
        f"{type(value).__name__}"
    )


def _parse_int(value: Any, *, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field} must be a JSON integer")
    return value


def _require_power_of_two(value: int, *, field: str) -> None:
    if value <= 0 or value & (value - 1):
        raise ValueError(f"{field} must be a positive power of two")


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def _json_exact(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if isinstance(value, dict):
        return {str(key): _json_exact(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_exact(item) for item in value]
    return value


def _write_json(value: Any, output: Path | None) -> None:
    text = json.dumps(_json_exact(value), indent=2, sort_keys=True) + "\n"
    if output is None:
        sys.stdout.write(text)
    else:
        output.write_text(text, encoding="utf-8")


def normalize_pmf(raw_pmf: Mapping[Any, Any]) -> dict[int, Fraction]:
    """Parse, combine, and normalize-check an exact finite coefficient PMF."""

    pmf: dict[int, Fraction] = {}
    for raw_value, raw_mass in raw_pmf.items():
        try:
            value = int(raw_value)
        except (TypeError, ValueError) as error:
            raise ValueError(f"PMF key {raw_value!r} is not an integer") from error
        if str(value) != str(raw_value) and not isinstance(raw_value, int):
            raise ValueError(f"PMF key {raw_value!r} is not a canonical integer")
        mass = parse_fraction(raw_mass, field=f"pmf[{raw_value!r}]")
        if mass < 0:
            raise ValueError("PMF masses must be nonnegative")
        pmf[value] = pmf.get(value, Fraction(0)) + mass
    if not pmf:
        raise ValueError("PMF must have finite nonempty support")
    total = sum(pmf.values(), Fraction(0))
    if total != 1:
        raise ValueError(f"PMF must sum exactly to 1, not {_fraction_text(total)}")
    return {value: mass for value, mass in sorted(pmf.items()) if mass}


def residue_collision_mass(pmf: Mapping[int, Fraction], exponent: int) -> Fraction:
    """Return ``sum_r Pr[E=r mod 2^exponent]^2`` exactly."""

    if exponent < 0:
        raise ValueError("the divisibility exponent must be nonnegative")
    modulus = 1 << exponent
    residue_mass: dict[int, Fraction] = {}
    for value, mass in pmf.items():
        residue = value % modulus
        residue_mass[residue] = residue_mass.get(residue, Fraction(0)) + mass
    return sum((mass * mass for mass in residue_mass.values()), Fraction(0))


def exact_error_table(pmf: Mapping[int, Fraction], modulus_bits: int) -> dict[str, Any]:
    """Compute every exact difference-divisibility table entry through ``K``."""

    if modulus_bits <= 0:
        raise ValueError("K must be positive")
    probabilities = [
        residue_collision_mass(pmf, exponent)
        for exponent in range(modulus_bits + 1)
    ]
    rows: list[dict[str, Any]] = []
    for exponent, probability in enumerate(probabilities):
        row: dict[str, Any] = {
            "e": exponent,
            "p": probability,
            "zeta_hat": (1 << exponent) * probability,
        }
        if exponent < modulus_bits:
            if probability == 0:
                raise ValueError("a finite normalized PMF cannot have zero collision mass")
            numerator = 2 * probabilities[exponent + 1] - probability
            beta = numerator / probability
            if not 0 <= beta <= 1:
                raise AssertionError("IID difference parity bias escaped [0,1]")
            row.update(
                {
                    "p_next": probabilities[exponent + 1],
                    "bias_numerator": numerator,
                    "beta": beta,
                }
            )
        rows.append(row)
    return {
        "K": modulus_bits,
        "pmf": {str(value): mass for value, mass in pmf.items()},
        "normalization": sum(pmf.values(), Fraction(0)),
        "rows": rows,
    }


def lucas_hasse_weight_value(
    degree: int, depth: int, point: Fraction, *, trace: list[dict[str, Any]] | None = None
) -> Fraction:
    """Evaluate ``W_(degree,depth)(point)`` by the exact Lucas recursion."""

    _require_power_of_two(degree, field="n")
    if not 0 <= depth <= degree:
        raise ValueError("v must lie in [0,n]")
    if degree == 1:
        result = Fraction(1) if depth == 0 else 1 + point
        if trace is not None:
            trace.append(
                {
                    "n": degree,
                    "v": depth,
                    "point": point,
                    "branch": "base",
                    "result": result,
                }
            )
        return result

    half = degree // 2
    if depth <= half:
        child_point = point * point
        child = lucas_hasse_weight_value(half, depth, child_point, trace=trace)
        result = child
        branch = "low"
        child_depth = depth
    else:
        denominator = 1 + point * point
        if denominator == 0:
            raise ZeroDivisionError("the high Lucas transform denominator vanished")
        child_point = 2 * point / denominator
        child_depth = depth - half
        child = lucas_hasse_weight_value(
            half, child_depth, child_point, trace=trace
        )
        result = denominator**half * child
        branch = "high"
    if trace is not None:
        trace.append(
            {
                "n": degree,
                "v": depth,
                "point": point,
                "branch": branch,
                "child_n": half,
                "child_v": child_depth,
                "child_point": child_point,
                "result": result,
            }
        )
    return result


def exact_local_factor(
    *, degree: int, exponent: int, depth: int, p: Fraction, p_next: Fraction
) -> dict[str, Any]:
    """Evaluate the exact Hasse factor and its proved product upper bound."""

    _require_power_of_two(degree, field="n")
    if exponent < 0 or not 0 <= depth <= degree:
        raise ValueError("expected e >= 0 and 0 <= v <= n")
    if p <= 0 or p_next < 0 or p_next > p:
        raise ValueError("expected 0 < p_next <= p (allowing p_next=0 only formally)")
    beta = (2 * p_next - p) / p
    if not 0 <= beta <= 1:
        raise ValueError("the supplied p and p_next do not give an IID bias in [0,1]")
    trace: list[dict[str, Any]] = []
    enumerator = lucas_hasse_weight_value(degree, depth, beta, trace=trace)
    zeta = (1 << exponent) * p
    zeta_next = (1 << (exponent + 1)) * p_next
    exact = zeta**degree * enumerator
    product_bound = zeta ** (degree - depth) * zeta_next**depth
    if exact > product_bound:
        raise AssertionError("exact Hasse factor exceeded its formal product bound")
    return {
        "n": degree,
        "e": exponent,
        "v": depth,
        "p": p,
        "p_next": p_next,
        "beta": beta,
        "zeta_hat": zeta,
        "zeta_hat_next": zeta_next,
        "weight_enumerator_value": enumerator,
        "exact_normalized_factor": exact,
        "two_level_product_bound": product_bound,
        "exact_le_product_bound": True,
        "lucas_trace_postorder": trace,
    }


def _truncated_binomial_row(power: int, cutoff: int) -> list[int]:
    return [math.comb(power, index) for index in range(min(power, cutoff) + 1)]


def hasse_kernel_weight_coefficients(
    degree: int,
    depth: int,
    cutoff: int,
    *,
    _cache: dict[tuple[int, int, int], tuple[int, ...]] | None = None,
) -> tuple[int, ...]:
    """Return coefficients of ``A_(n,v)`` through ``z^cutoff``.

    This is the dual Lucas recursion.  Truncation is exact because neither
    branch can lower Hamming weight.
    """

    _require_power_of_two(degree, field="n")
    if not 0 <= depth <= degree:
        raise ValueError("v must lie in [0,n]")
    if cutoff < 0:
        raise ValueError("the Hamming-weight cutoff must be nonnegative")
    cutoff = min(cutoff, degree)
    cache = {} if _cache is None else _cache
    key = (degree, depth, cutoff)
    if key in cache:
        return cache[key]
    if degree == 1:
        result = (1,) if cutoff == 0 or depth == 1 else (1, 1)
        cache[key] = result
        return result

    half = degree // 2
    coefficients = [0] * (cutoff + 1)
    if depth <= half:
        child = hasse_kernel_weight_coefficients(
            half, depth, cutoff, _cache=cache
        )
        for weight, multiplicity in enumerate(child):
            if multiplicity == 0 or weight > cutoff:
                continue
            residual_power = half - weight
            maximum_index = min(residual_power, (cutoff - weight) // 2)
            scale = multiplicity * (1 << weight)
            for index in range(maximum_index + 1):
                coefficients[weight + 2 * index] += (
                    scale * math.comb(residual_power, index)
                )
    else:
        child = hasse_kernel_weight_coefficients(
            half, depth - half, cutoff // 2, _cache=cache
        )
        for weight, multiplicity in enumerate(child):
            if 2 * weight <= cutoff:
                coefficients[2 * weight] = multiplicity
    result = tuple(coefficients)
    cache[key] = result
    return result


def exact_hasse_valuation_weight_count(
    degree: int,
    valuation: int,
    weight: int,
    *,
    _cache: dict[tuple[int, int, int], tuple[int, ...]] | None = None,
) -> int:
    if not 0 <= valuation < degree:
        raise ValueError("a nonzero binary word valuation must lie in [0,n)")
    if not 0 <= weight <= degree:
        return 0
    cache = {} if _cache is None else _cache
    at_level = hasse_kernel_weight_coefficients(
        degree, valuation, weight, _cache=cache
    )[weight]
    at_next = hasse_kernel_weight_coefficients(
        degree, valuation + 1, weight, _cache=cache
    )[weight]
    if at_next > at_level:
        raise AssertionError("nested Hasse kernel counts are not monotone")
    return at_level - at_next


def fixed_weight_ternary_histogram(
    degree: int, secret_weight: int, valuations: Iterable[int]
) -> dict[str, Any]:
    """Return exact signed-ternary pair counts for selected Hasse valuations."""

    _require_power_of_two(degree, field="n")
    if not 0 <= secret_weight <= degree:
        raise ValueError("w must lie in [0,n]")
    selected = sorted(set(valuations))
    if any(valuation < 0 or valuation >= degree for valuation in selected):
        raise ValueError("every requested valuation must lie in [0,n)")
    cache: dict[tuple[int, int, int], tuple[int, ...]] = {}
    support_changing: list[dict[str, Any]] = []
    same_support: list[dict[str, Any]] = []

    for a in range(1, min(secret_weight, degree - secret_weight) + 1):
        binary_weight = 2 * a
        combinatorial = (
            (1 << (2 * secret_weight))
            * math.comb(2 * a, a)
            * math.comb(degree - 2 * a, secret_weight - a)
        )
        for valuation in selected:
            primitive_count = exact_hasse_valuation_weight_count(
                degree,
                valuation,
                binary_weight,
                _cache=cache,
            )
            count = combinatorial * primitive_count
            if count:
                support_changing.append(
                    {
                        "e": 0,
                        "v": valuation,
                        "a": a,
                        "binary_weight": binary_weight,
                        "primitive_word_count": primitive_count,
                        "ordered_pair_count": count,
                    }
                )

    for flips in range(1, secret_weight + 1):
        combinatorial = (
            (1 << secret_weight)
            * math.comb(degree - flips, secret_weight - flips)
        )
        for valuation in selected:
            primitive_count = exact_hasse_valuation_weight_count(
                degree, valuation, flips, _cache=cache
            )
            count = combinatorial * primitive_count
            if count:
                same_support.append(
                    {
                        "e": 1,
                        "v": valuation,
                        "b": flips,
                        "binary_weight": flips,
                        "primitive_word_count": primitive_count,
                        "ordered_pair_count": count,
                    }
                )

    secret_cardinality = (1 << secret_weight) * math.comb(degree, secret_weight)
    covered = (
        secret_cardinality
        + sum(row["ordered_pair_count"] for row in support_changing)
        + sum(row["ordered_pair_count"] for row in same_support)
    )
    full = selected == list(range(degree))
    if full and covered != secret_cardinality * secret_cardinality:
        raise AssertionError("full ternary pair histogram failed to normalize")
    return {
        "n": degree,
        "w": secret_weight,
        "selected_valuations": selected,
        "secret_cardinality": secret_cardinality,
        "zero_difference_count": secret_cardinality,
        "support_changing": support_changing,
        "same_support": same_support,
        "covered_ordered_pair_count": covered,
        "full_histogram": full,
        "expected_full_ordered_pair_count": secret_cardinality**2,
    }


def _coerce_polynomial(raw: Any, degree: int, *, field: str) -> list[int]:
    if isinstance(raw, bool):
        raise ValueError(f"{field} must be an integer or coefficient list")
    if isinstance(raw, int):
        return [raw] + [0] * (degree - 1)
    if not isinstance(raw, list) or len(raw) != degree:
        raise ValueError(f"{field} must contain exactly n integer coefficients")
    return [_parse_int(value, field=f"{field}[{index}]") for index, value in enumerate(raw)]


def _reduce_coefficients(polynomial: Sequence[int], modulus: int) -> list[int]:
    return [coefficient % modulus for coefficient in polynomial]


def negacyclic_add(left: Sequence[int], right: Sequence[int], modulus: int) -> list[int]:
    if len(left) != len(right):
        raise ValueError("polynomial lengths differ")
    return [(a + b) % modulus for a, b in zip(left, right)]


def negacyclic_sub(left: Sequence[int], right: Sequence[int], modulus: int) -> list[int]:
    if len(left) != len(right):
        raise ValueError("polynomial lengths differ")
    return [(a - b) % modulus for a, b in zip(left, right)]


def negacyclic_mul(left: Sequence[int], right: Sequence[int], modulus: int) -> list[int]:
    """Multiply modulo ``X^n+1`` and the coefficient modulus."""

    if len(left) != len(right):
        raise ValueError("polynomial lengths differ")
    degree = len(left)
    result = [0] * degree
    for left_index, left_value in enumerate(left):
        if left_value % modulus == 0:
            continue
        for right_index, right_value in enumerate(right):
            total_index = left_index + right_index
            if total_index < degree:
                result[total_index] += left_value * right_value
            else:
                result[total_index - degree] -= left_value * right_value
    return _reduce_coefficients(result, modulus)


def _v2_mod_power_of_two(value: int, modulus_bits: int) -> int:
    residue = value % (1 << modulus_bits)
    if residue == 0:
        return modulus_bits
    return (residue & -residue).bit_length() - 1


def pi_basis_coefficients(
    polynomial: Sequence[int], modulus_bits: int
) -> list[int]:
    """Expand a degree-below-``n`` representative in ``pi=1-X``."""

    if modulus_bits <= 0:
        raise ValueError("K must be positive")
    degree = len(polynomial)
    _require_power_of_two(degree, field="polynomial length")
    modulus = 1 << modulus_bits
    result: list[int] = []
    for pi_degree in range(degree):
        coefficient = sum(
            polynomial[x_degree] * math.comb(x_degree, pi_degree)
            for x_degree in range(pi_degree, degree)
        )
        if pi_degree & 1:
            coefficient = -coefficient
        result.append(coefficient % modulus)
    return result


def pi_valuation(polynomial: Sequence[int], modulus_bits: int) -> int:
    """Exact valuation in ``(Z/2^K)[X]/(X^n+1)`` for power-of-two ``n``.

    In the ``pi=1-X`` basis, the initial degree of coefficient ``c_j`` is
    ``j + n*v2(c_j)`` because ``2`` is a unit times ``pi^n``.  Different
    ``j`` have different residues modulo ``n``, so their initial terms cannot
    cancel.  Zero is assigned the terminal value ``K*n``.
    """

    degree = len(polynomial)
    _require_power_of_two(degree, field="polynomial length")
    if modulus_bits <= 0:
        raise ValueError("K must be positive")
    terminal = modulus_bits * degree
    return min(
        terminal,
        min(
            index + degree * _v2_mod_power_of_two(coefficient, modulus_bits)
            for index, coefficient in enumerate(
                pi_basis_coefficients(polynomial, modulus_bits)
            )
        ),
    )


def valuation_record(polynomial: Sequence[int], modulus_bits: int) -> dict[str, Any]:
    degree = len(polynomial)
    modulus = 1 << modulus_bits
    reduced = _reduce_coefficients(polynomial, modulus)
    pi_coefficients = pi_basis_coefficients(reduced, modulus_bits)
    order = pi_valuation(reduced, modulus_bits)
    terminal = modulus_bits * degree
    record: dict[str, Any] = {
        "K": modulus_bits,
        "n": degree,
        "coefficients_mod_2^K": reduced,
        "pi_basis_coefficients_mod_2^K": pi_coefficients,
        "valuation": order,
        "terminal": order == terminal,
    }
    if order < terminal:
        exponent, depth = divmod(order, degree)
        record.update({"e": exponent, "v": depth})
    return record


def complete_descriptor(
    *,
    secret_left: Sequence[int],
    secret_right: Sequence[int],
    row_mask: Sequence[int],
    common_factor: Sequence[int],
    gadget: Sequence[int],
    modulus_bits: int,
) -> dict[str, Any]:
    """Compute the complete descriptor polynomial and capped valuation sum."""

    degree = len(secret_left)
    if any(
        len(polynomial) != degree
        for polynomial in (secret_right, row_mask, common_factor, gadget)
    ):
        raise ValueError("all descriptor polynomials must have the same length")
    modulus = 1 << modulus_bits
    difference = negacyclic_sub(secret_right, secret_left, modulus)
    secret_sum = negacyclic_add(secret_right, secret_left, modulus)
    multiplier = negacyclic_add(
        row_mask,
        negacyclic_mul(
            negacyclic_mul(common_factor, gadget, modulus),
            secret_sum,
            modulus,
        ),
        modulus,
    )
    delta = negacyclic_mul(difference, multiplier, modulus)
    difference_order = pi_valuation(difference, modulus_bits)
    multiplier_order = pi_valuation(multiplier, modulus_bits)
    delta_order = pi_valuation(delta, modulus_bits)
    predicted = min(modulus_bits * degree, difference_order + multiplier_order)
    if delta_order != predicted:
        raise AssertionError("descriptor valuation did not compose by capped addition")
    return {
        "K": modulus_bits,
        "n": degree,
        "secret_difference": difference,
        "secret_sum": secret_sum,
        "multiplier": multiplier,
        "descriptor_difference": delta,
        "secret_difference_valuation": difference_order,
        "multiplier_valuation": multiplier_order,
        "descriptor_valuation": delta_order,
        "capped_sum": predicted,
        "composition_matches": True,
        "valuation_record": valuation_record(delta, modulus_bits),
    }


def aggregate_joint_histogram(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Aggregate exact row factors against a full joint tuple histogram."""

    raw_rows = payload.get("rows")
    raw_histogram = payload.get("histogram")
    if not isinstance(raw_rows, list) or not raw_rows:
        raise ValueError("rows must be a nonempty list")
    if not isinstance(raw_histogram, list) or not raw_histogram:
        raise ValueError("histogram must be a nonempty list")
    row_tables: list[dict[tuple[int, int], Fraction]] = []
    for row_index, raw_row in enumerate(raw_rows):
        if not isinstance(raw_row, Mapping) or not isinstance(raw_row.get("factors"), Mapping):
            raise ValueError(f"rows[{row_index}].factors must be an object")
        table: dict[tuple[int, int], Fraction] = {}
        for raw_key, raw_factor in raw_row["factors"].items():
            pieces = str(raw_key).split(",")
            if len(pieces) != 2:
                raise ValueError("factor keys must have the form e,v")
            key = (int(pieces[0]), int(pieces[1]))
            table[key] = parse_fraction(
                raw_factor, field=f"rows[{row_index}].factors[{raw_key!r}]"
            )
        row_tables.append(table)

    total_weight = Fraction(0)
    total_overlap = Fraction(0)
    cells: list[dict[str, Any]] = []
    for cell_index, raw_cell in enumerate(raw_histogram):
        if not isinstance(raw_cell, Mapping):
            raise ValueError(f"histogram[{cell_index}] must be an object")
        weight = parse_fraction(raw_cell.get("weight"), field=f"histogram[{cell_index}].weight")
        if weight < 0:
            raise ValueError("histogram weights must be nonnegative")
        raw_tuple = raw_cell.get("tuple")
        if not isinstance(raw_tuple, list) or len(raw_tuple) != len(row_tables):
            raise ValueError("every tuple must contain one [e,v] pair per row")
        row_product = Fraction(1)
        tuple_value: list[list[int]] = []
        for row_index, raw_pair in enumerate(raw_tuple):
            if not isinstance(raw_pair, list) or len(raw_pair) != 2:
                raise ValueError("tuple entries must be [e,v] integer pairs")
            key = (
                _parse_int(raw_pair[0], field="tuple exponent"),
                _parse_int(raw_pair[1], field="tuple depth"),
            )
            try:
                factor = row_tables[row_index][key]
            except KeyError as error:
                raise ValueError(
                    f"row {row_index} has no factor for tuple cell {key}"
                ) from error
            row_product *= factor
            tuple_value.append([key[0], key[1]])
        contribution = weight * row_product
        total_weight += weight
        total_overlap += contribution
        cells.append(
            {
                "tuple": tuple_value,
                "weight": weight,
                "row_factor_product": row_product,
                "contribution": contribution,
            }
        )
    if total_weight != 1:
        raise ValueError(
            f"joint histogram must sum exactly to 1, not {_fraction_text(total_weight)}"
        )
    return {
        "row_count": len(row_tables),
        "cell_count": len(cells),
        "histogram_normalization": total_weight,
        "cells": cells,
        "exact_overlap": total_overlap,
    }


def _negacyclic_pow(base: Sequence[int], exponent: int, modulus: int) -> list[int]:
    degree = len(base)
    result = [1] + [0] * (degree - 1)
    power = list(base)
    while exponent:
        if exponent & 1:
            result = negacyclic_mul(result, power, modulus)
        power = negacyclic_mul(power, power, modulus)
        exponent >>= 1
    return result


def run_self_tests() -> dict[str, Any]:
    """Run deterministic identities that exercise every command core."""

    pmf = normalize_pmf({"-1": "1/4", "0": "1/2", "1": "1/4"})
    table = exact_error_table(pmf, 3)
    assert table["rows"][0]["p"] == 1
    assert table["rows"][0]["beta"] == 0
    local = exact_local_factor(
        degree=8,
        exponent=0,
        depth=3,
        p=Fraction(1),
        p_next=Fraction(1, 2),
    )
    assert local["exact_normalized_factor"] <= local["two_level_product_bound"]

    # Compare the recursion to direct row-space enumeration at a small size.
    point = Fraction(3, 7)
    for depth in range(9):
        direct = Fraction(0)
        for message in range(1 << depth):
            word = 0
            for row in range(depth):
                if message >> row & 1:
                    for column in range(8):
                        if column & row == row:
                            word ^= 1 << column
            direct += point ** word.bit_count()
        assert lucas_hasse_weight_value(8, depth, point) == direct

    histogram = fixed_weight_ternary_histogram(8, 2, range(8))
    assert histogram["full_histogram"]
    assert histogram["covered_ordered_pair_count"] == histogram[
        "expected_full_ordered_pair_count"
    ]

    for modulus_bits in (1, 2, 3):
        modulus = 1 << modulus_bits
        uniformizer = [1, -1, 0, 0]
        for order in range(modulus_bits * 4 + 1):
            power = _negacyclic_pow(uniformizer, order, modulus)
            assert pi_valuation(power, modulus_bits) == min(order, modulus_bits * 4)

    descriptor = complete_descriptor(
        secret_left=[1, 0, -1, 0],
        secret_right=[0, 1, 0, -1],
        row_mask=[2, 0, 0, 0],
        common_factor=[1, 0, 0, 0],
        gadget=[1, 1, 0, 0],
        modulus_bits=4,
    )
    assert descriptor["composition_matches"]

    aggregate = aggregate_joint_histogram(
        {
            "rows": [
                {"factors": {"0,0": "2", "0,1": "3"}},
                {"factors": {"0,0": "5", "0,1": "7"}},
            ],
            "histogram": [
                {"tuple": [[0, 0], [0, 1]], "weight": "1/3"},
                {"tuple": [[0, 1], [0, 0]], "weight": "2/3"},
            ],
        }
    )
    assert aggregate["exact_overlap"] == Fraction(44, 3)
    return {
        "status": "ok",
        "checks": [
            "IID difference table",
            "Lucas recursion against direct enumeration",
            "kernel recursion and full ternary histogram normalization",
            "uniformizer powers through the truncated chain",
            "complete descriptor capped valuation",
            "joint tuple aggregation",
        ],
    }


def _payload_for_command(command: str, payload: Mapping[str, Any]) -> Any:
    if command == "error-table":
        raw_pmf = payload.get("pmf")
        if not isinstance(raw_pmf, Mapping):
            raise ValueError("pmf must be a JSON object")
        return exact_error_table(
            normalize_pmf(raw_pmf),
            _parse_int(payload.get("K"), field="K"),
        )
    if command == "local":
        return exact_local_factor(
            degree=_parse_int(payload.get("n"), field="n"),
            exponent=_parse_int(payload.get("e"), field="e"),
            depth=_parse_int(payload.get("v"), field="v"),
            p=parse_fraction(payload.get("p"), field="p"),
            p_next=parse_fraction(payload.get("p_next"), field="p_next"),
        )
    if command == "secret-hist":
        degree = _parse_int(payload.get("n"), field="n")
        raw_valuations = payload.get("valuations")
        if raw_valuations is None:
            if degree > 256:
                raise ValueError(
                    "for n>256, provide an explicit valuations list to avoid an "
                    "accidental large certificate"
                )
            valuations: Iterable[int] = range(degree)
        elif isinstance(raw_valuations, list):
            valuations = [
                _parse_int(value, field="valuations entry")
                for value in raw_valuations
            ]
        else:
            raise ValueError("valuations must be a JSON integer list")
        return fixed_weight_ternary_histogram(
            degree,
            _parse_int(payload.get("w"), field="w"),
            valuations,
        )
    if command == "ord":
        modulus_bits = _parse_int(payload.get("K"), field="K")
        raw_polynomial = payload.get("polynomial")
        if not isinstance(raw_polynomial, list):
            raise ValueError("polynomial must be an integer coefficient list")
        polynomial = [
            _parse_int(value, field=f"polynomial[{index}]")
            for index, value in enumerate(raw_polynomial)
        ]
        return valuation_record(polynomial, modulus_bits)
    if command == "descriptor":
        degree = _parse_int(payload.get("n"), field="n")
        modulus_bits = _parse_int(payload.get("K"), field="K")
        return complete_descriptor(
            secret_left=_coerce_polynomial(payload.get("s"), degree, field="s"),
            secret_right=_coerce_polynomial(payload.get("t"), degree, field="t"),
            row_mask=_coerce_polynomial(payload.get("z"), degree, field="z"),
            common_factor=_coerce_polynomial(payload.get("f"), degree, field="f"),
            gadget=_coerce_polynomial(payload.get("g"), degree, field="g"),
            modulus_bits=modulus_bits,
        )
    if command == "aggregate":
        return aggregate_joint_histogram(payload)
    raise ValueError(f"unknown command {command!r}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run exact deterministic regression checks",
    )
    subparsers = parser.add_subparsers(dest="command")
    for name in (
        "error-table",
        "local",
        "secret-hist",
        "ord",
        "descriptor",
        "aggregate",
    ):
        command = subparsers.add_parser(name)
        command.add_argument("input", type=Path, help="exact JSON input")
        command.add_argument("--output", type=Path, help="write JSON here")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.self_test:
            if args.command is not None:
                parser.error("--self-test cannot be combined with a command")
            _write_json(run_self_tests(), None)
            return 0
        if args.command is None:
            parser.error("choose a command or pass --self-test")
        payload = _load_json(args.input)
        if not isinstance(payload, Mapping):
            raise ValueError("the input root must be a JSON object")
        _write_json(_payload_for_command(args.command, payload), args.output)
        return 0
    except (ValueError, ZeroDivisionError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
