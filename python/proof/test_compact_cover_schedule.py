#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compact_cover_schedule import (
    general_baby_giant,
    index_to_sequence,
    thin_power_of_two_schedule,
)


class CompactCoverScheduleTest(unittest.TestCase):
    def test_mixed_radix(self) -> None:
        self.assertEqual(index_to_sequence(1, (2, 3)), (0, 0))
        self.assertEqual(index_to_sequence(6, (2, 3)), (1, 2))

    def test_general_baby_giant_matches_magma_heuristic(self) -> None:
        self.assertEqual(general_baby_giant((2, 16_384)), ((2, 91), (1, 181)))

    def test_target_schedule(self) -> None:
        schedule = thin_power_of_two_schedule()
        self.assertEqual(schedule.frobenius_order, 2)
        self.assertEqual(schedule.slot_count, 32_768)
        self.assertEqual(schedule.dimension_sizes, (16_384, 2))
        self.assertEqual(schedule.baby_product, 182)
        self.assertEqual(schedule.giant_product, 181)
        self.assertEqual(len(schedule.switch_key_exponents), 362)
        self.assertEqual(schedule.generated_group_size, 65_536)
        self.assertEqual(schedule.peak_live_ciphertexts, 368)


if __name__ == "__main__":
    unittest.main()
