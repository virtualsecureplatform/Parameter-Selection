#!/usr/bin/env python3
"""Regression tests for the exact two-level Smith certificate evaluator."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import unittest

from twosmith_exact_evaluator import (
    aggregate_joint_histogram,
    exact_error_table,
    exact_hasse_valuation_weight_count,
    exact_local_factor,
    fixed_weight_ternary_histogram,
    hasse_kernel_weight_coefficients,
    lucas_hasse_weight_value,
    negacyclic_mul,
    normalize_pmf,
    parse_fraction,
    pi_valuation,
)


class ExactRationalInputTests(unittest.TestCase):
    def test_fraction_parser_rejects_inexact_values(self) -> None:
        self.assertEqual(parse_fraction("-7/12"), Fraction(-7, 12))
        self.assertEqual(parse_fraction(5), Fraction(5))
        with self.assertRaises(ValueError):
            parse_fraction(0.5)
        with self.assertRaises(ValueError):
            parse_fraction("0.5")

    def test_iid_difference_table(self) -> None:
        pmf = normalize_pmf({"-1": "1/4", "0": "1/2", "1": "1/4"})
        report = exact_error_table(pmf, 3)
        rows = report["rows"]
        self.assertEqual(rows[0]["p"], 1)
        self.assertEqual(rows[0]["beta"], 0)
        for row in rows[:-1]:
            self.assertGreaterEqual(row["beta"], 0)
            self.assertLessEqual(row["beta"], 1)
            self.assertEqual(
                row["bias_numerator"], 2 * row["p_next"] - row["p"]
            )


class HasseRecursionTests(unittest.TestCase):
    @staticmethod
    def direct_dual_enumerator(degree: int, depth: int, point: Fraction) -> Fraction:
        total = Fraction(0)
        for message in range(1 << depth):
            word = 0
            for row in range(depth):
                if message >> row & 1:
                    for column in range(degree):
                        if column & row == row:
                            word ^= 1 << column
            total += point ** word.bit_count()
        return total

    def test_lucas_recursion_matches_native_rows(self) -> None:
        point = Fraction(2, 5)
        for degree in (1, 2, 4, 8):
            for depth in range(degree + 1):
                self.assertEqual(
                    lucas_hasse_weight_value(degree, depth, point),
                    self.direct_dual_enumerator(degree, depth, point),
                )

    def test_kernel_coefficients_partition_nonzero_words(self) -> None:
        degree = 8
        cache: dict[tuple[int, int, int], tuple[int, ...]] = {}
        for depth in range(degree + 1):
            direct = [0] * (degree + 1)
            for word in range(1 << degree):
                in_kernel = all(
                    sum(
                        ((column & row) == row) * ((word >> column) & 1)
                        for column in range(degree)
                    )
                    % 2
                    == 0
                    for row in range(depth)
                )
                if in_kernel:
                    direct[word.bit_count()] += 1
            self.assertEqual(
                hasse_kernel_weight_coefficients(
                    degree, depth, degree, _cache=cache
                ),
                tuple(direct),
            )
        for weight in range(1, degree + 1):
            exact_total = sum(
                exact_hasse_valuation_weight_count(
                    degree, valuation, weight, _cache=cache
                )
                for valuation in range(degree)
            )
            self.assertEqual(exact_total, math_comb(degree, weight))
        full_cube = hasse_kernel_weight_coefficients(degree, 0, degree)
        self.assertEqual(full_cube, tuple(math_comb(degree, k) for k in range(degree + 1)))

    def test_local_factor_is_below_product_bound(self) -> None:
        result = exact_local_factor(
            degree=8,
            exponent=2,
            depth=5,
            p=Fraction(3, 8),
            p_next=Fraction(1, 4),
        )
        self.assertLessEqual(
            result["exact_normalized_factor"], result["two_level_product_bound"]
        )


def math_comb(n: int, k: int) -> int:
    # Kept local so all imported evaluator names above are intentional API tests.
    import math

    return math.comb(n, k)


class TernaryHistogramTests(unittest.TestCase):
    def test_full_histogram_normalizes(self) -> None:
        for degree, weight in ((4, 1), (4, 2), (8, 2)):
            report = fixed_weight_ternary_histogram(degree, weight, range(degree))
            self.assertTrue(report["full_histogram"])
            self.assertEqual(
                report["covered_ordered_pair_count"],
                report["expected_full_ordered_pair_count"],
            )


class ValuationTests(unittest.TestCase):
    def test_small_rings_satisfy_capped_product_law_exhaustively(self) -> None:
        for degree, modulus_bits in ((1, 2), (2, 2), (4, 1)):
            modulus = 1 << modulus_bits
            elements = list(product(range(modulus), repeat=degree))
            terminal = modulus_bits * degree
            for left in elements:
                left_order = pi_valuation(left, modulus_bits)
                for right in elements:
                    right_order = pi_valuation(right, modulus_bits)
                    product_value = negacyclic_mul(left, right, modulus)
                    self.assertEqual(
                        pi_valuation(product_value, modulus_bits),
                        min(terminal, left_order + right_order),
                    )


class JointAggregationTests(unittest.TestCase):
    def test_correlated_tuple_cells_are_aggregated_jointly(self) -> None:
        report = aggregate_joint_histogram(
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
        self.assertEqual(report["histogram_normalization"], 1)
        self.assertEqual(report["exact_overlap"], Fraction(44, 3))


if __name__ == "__main__":
    unittest.main()
