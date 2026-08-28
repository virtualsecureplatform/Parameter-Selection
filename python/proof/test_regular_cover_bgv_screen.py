#!/usr/bin/env python3

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from regular_cover_bgv_screen import ScreenInput, minimum_public_key_rows, screen


class RegularCoverBGVScreenTest(unittest.TestCase):
    def test_leftover_threshold_is_tight(self) -> None:
        log_card = 1000
        rows = minimum_public_key_rows(log_card, 128)
        self.assertEqual(rows, 2254)
        result = screen(
            ScreenInput(
                degree=10,
                group_size=10,
                q_bits=10,
                public_key_rows=rows,
            )
        )
        self.assertTrue(result.leftover_ok)
        self.assertEqual(result.leftover_log2_bound, -128.0)

    def test_cover_costs_scale_with_group(self) -> None:
        result = screen(
            ScreenInput(
                degree=16,
                group_size=4,
                q_bits=20,
                evaluation_rows=10,
                pivot_attempts=8,
                native_failure_bits=140,
                native_input_noise_bits=30,
                native_fresh_noise_bits=20,
            )
        )
        self.assertEqual(result.cover_log2_cardinality, 1280)
        self.assertEqual(result.ciphertext_storage_bits, 2560)
        self.assertEqual(result.cover_failure_bits, 138.0)
        self.assertTrue(result.refresh_contract_ok)
        expected_rows = result.public_key_rows + 10 + 8
        self.assertEqual(result.base_rlwe_samples, 4 * expected_rows)

    def test_pivot_condition_and_failure(self) -> None:
        result = screen(
            ScreenInput(degree=8, group_size=2, q_bits=8, pivot_attempts=20)
        )
        self.assertTrue(result.pivot_condition_ok)
        self.assertLessEqual(result.pivot_failure_log2, -20.0)
        self.assertGreater(result.unit_probability, 0.5)

    def test_union_bound_can_fail(self) -> None:
        result = screen(
            ScreenInput(
                degree=8,
                group_size=8,
                q_bits=10,
                native_failure_bits=128,
            )
        )
        self.assertFalse(result.failure_ok)
        self.assertEqual(result.cover_failure_bits, 125.0)

    def test_unspecified_failure_is_not_a_certificate(self) -> None:
        result = screen(ScreenInput(degree=8, group_size=2, q_bits=8))
        self.assertTrue(math.isinf(result.cover_failure_bits))
        self.assertFalse(result.failure_ok)


if __name__ == "__main__":
    unittest.main()
